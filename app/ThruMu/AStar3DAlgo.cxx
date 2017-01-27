#include "AStar3DAlgo.h"
#include <vector>
#include <cmath>

// larcv
#include "UBWireTool/UBWireTool.h"
#include ""

namespace larlitecv {

  
  std::vector<AStar3DNode> AStar3DAlgo::findpath( const std::vector<larcv::Image2D>& img, const int start_row, const int goal_row, 
                                                  const std::vector<int>& goal_cols, const std::vector<float>& thresh ) {
    
    const larcv::ImageMeta& meta = img.front().meta();

    if ( verbose<=1 )
      std::cout << "[[ASTAR 3D ALGO START]]"

    // turn image pixel coordinates into a 3D start and goal point
    std::vector<int> start_wids(start_cols.size());
    std::vector<int> goal_wids(goal_cols.size());    
    for (size_t p=0; p<goal_cols.size()) {
      start_wids[p] = img_v.at(p).meta().pos_x( start_cols[p] );
      goal_wids[p]  = img_v.at(p).meta().pos_x( goal_cols[p] );      
    }
    float start_tri = 0.;
    int start_crosses = 0;
    std::vector< float > poszy_start(2,0.0);
    larcv::UBWireTool::wireIntersection( start_wids, poszy_start, start_tri, start_crosses );

    float goal_tri = 0.;
    int goal_crosses = 0;
    std::vector< float > poszy_goal(2,0.0);
    larcv::UBWireTool::wireIntersection( goal_wids, poszy_goal, goal_tri, goal_crosses );

    std::vector<float> startpos(3,0);
    startpos[0] = 0;
    startpos[1] = poszy_start[1];
    startpos[2] = poszy_start[0];

    std::vector<float> goalpos(3,0);
    goalpos[0] = fabs( goal_row-start_row )*cm_per_row;
    goalpos[1] = poszy_goal[1];
    goalpos[2] = poszy_goal[0];

    if ( start_crosses==0 || goal_crosses==0 ) {
      throw std::runtime_error("AStar3DAlgo::findpath[error] start or goal point not a good 3D space point.");
    }

    // next, define the lattice
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    float cm_per_wire = 0.3;
    float cm_per_row  = cm_per_tick*meta.pixel_height();
    float cm_per_col  = cm_per_wire*meta.pixel_width();

    std::vector<float> cm_per_pixel(3);
    cm_per_pixel[0] = cm_per_row;
    cm_per_pixel[1] = cm_per_col;
    cm_per_pixel[2] = cm_per_col;

    // we set t0 to be 0 on the lattice
    float tick0 = ( meta.pos_y(start_row) > meta.pos_y(goal_row) ) ? meta.pos_y(start_row) : meta.pos_y(goal_row); // why the largest row? image2d is time-reversed

    // now we define the lattice the search will go over
    std::vector<int> lattice_widths(3);
    lattice_widths[0] = abs( goal_row - start_row ); // rows
    lattice_widths[1] = (int)fabs( startpos[1] - goalpos[1] )/cm_per_col;
    lattice_widths[2] = (int)fabs( startpos[2] - goalpos[2] )/cm_per_col;

    // add padding
    for (int i=0; i<3; i++)
      lattice_widths[i] += 2*_config.lattice_padding;

    // define the origin of the lattice in detector space
    std::vector<float> origin_lattice(3);
    origin_lattice[0] = -_config.lattice_padding*cm_per_row; 
    origin_lattice[1] = startpos[1] - _config.lattice_padding*cm_per_col;
    origin_lattice[2] = startpos[2] - _config.lattice_padding*cm_per_col;

    // finally make the lattice
    Lattice lattice( origin_lattice, lattice_widths, cm_per_pixel );

    // we now make some definitions
    
    AStar3DNodePtrList openset;
    AStar3DNodePtrList closedset;
    std::map< PixPos_t, AStar3DNode* > position_lookup;

    // make the target node (window coorindates)
    AStar3DNode* goal = lattice.getNode( goalpos );

    // make starting node (window coordinates)
    AStar3DNode* start = lattice.getNode( startpos );
    start->gscore = 0.0;
    start->fscore = 0.0;
    // start fscore gets set to the heuristic
    for (int i=0; i<3; i++ ) {
      float di= goalpos[i]-startpos[i];
      start->fscore += di*di;
    }
    start->fscore = sqrt(fscore);

    // set the start direction to heads towards the goal
    start->dir3d.resize(3,0.0);
    for (int i=0; i<3; i++) {
      start->dir3d[i] = (goalpos[i]-startpos[i])/start->fscore;
    }

    // add the start node into the openset
    openset.addnode( start );
    
    if ( verbose>0 ) {
      std::cout << "start astar algo." << std::endl;
    }

    int neighborhood_size = _config.astar_neighborhood;
    int nsteps = -1;
    AStar3DNode* current = NULL;
    while ( openset.size()>0 ) {
      // get current
      current = openset.back();
      openset.pop_back();
      nsteps++;
      if ( verbose>2 || (verbose>0 && nsteps%100==0) ) {
        std::cout << "step=" << nsteps << ": get current node [ " << current->str() << " ] "
                  << "number of remaining nodes in openset=" << openset.size() << std::endl;
      }

      if ( *current==*goal ) {
        if ( verbose>0 ) {
          std::cout << "astar goal reached. finished." << std::endl;
        }
        break;        
      }

      // scan through neighors, and ID the candidate successor node
      evaluateNeighborNodes( current, start, goal, openset, closedset, neighborhood_size, 
        min_c, min_r, win_c, win_r, img, meta, use_bad_chs, position_lookup );
      //evaluateBadChNeighbors( current, start, goal, openset, closedset, 1, 
      //  min_c, min_r, win_c, win_r, img, meta, use_bad_chs, position_lookup );


      // finished with current node, put it on the closed set
      current->closed = true;
      closedset.addnode( current );

      openset.sort();

      if ( verbose>1 ) {
        std::cout << "update on sets: " << std::endl;
        std::cout << " openset: " << std::endl;
        for ( auto node : openset ) 
          std::cout << "  * (" << meta.pos_x(node->col+min_c) << "," << meta.pos_y(node->row+min_r) << ") fscore=" << node->fscore << " gscore=" << node->gscore << std::endl;
        std::cout << " nodes in closedset: " << closedset.size() << std::endl;
      }
          

      if ( verbose>1 && nsteps%100==0 ) {
        std::cout << "[enter] to continue." << std::endl;
        //std::cin.get();
      }

    }//end of while loop

    larcv::Image2D fscoreimg = visualizeScores( "fscore", img, min_c, min_r, win_c, win_r, position_lookup );
    larcv::Image2D gscoreimg = visualizeScores( "gscore", img, min_c, min_r, win_c, win_r, position_lookup );    
    m_visualizedimgs.clear();
    m_visualizedimgs.emplace_back( std::move(fscoreimg) );
    m_visualizedimgs.emplace_back( std::move(gscoreimg) );    

    bool path_completed = false;
    std::vector<AStar3DNode> path;
    if ( *current!=*goal ) {
      closedset.sort();
      current = closedset.back();
      if ( verbose>0)
        std::cout << "could not reach goal. best node: " 
                  << " (" << meta.pos_x( current->col + min_c ) << "," << meta.pos_y( current->row+min_r) << ")"
                  << " fscore=" << current->fscore << std::endl;
    }
    path = makeRecoPath( start, current, min_r, min_c, path_completed );

    // clean up
    if ( verbose>0 )
        std::cout << "clean up" << std::endl;
    for ( auto it_node : position_lookup ) {
      delete it_node.second;
    }
    delete goal;

    return path;
  }


  void AStar3DAlgo::evaluateNeighborNodes( AStar3DNode* current, const AStar3DNode* start, const AStar3DNode* goal,
    const int neighborhood_size, const float tick0, const std::vector<larcv::Image2D>& img, 
    AStar3DNodePtrList& openset, AStar3DNodePtrList& closedset, Lattice& lattice ) {

    const A3DPixPos_t& center = current->nodeid;

    int number_updates=0;

    // make a list of lattice points we have to evaluate
    std::vector<A3DPixPos_t> neighborhood_points;
    for (int du=-neighborhood_size; du<=neighborhood_size; du++) {
      for (int dv=-neighborhood_size; du<=neighborhood_size; du++) {
        for (int dw=-neighborhood_size; du<=neighborhood_size; du++) {
          int u = center[0]+du;
          int v = center[1]+dv;
          int w = center[2]+dw;
          if ( u<0 || u>=lattice.width[0] || v<0 || v>=lattice.width[1] || w<0 || w>=lattice.width[2] )
            continue;
          if (du==0 && dv==0 && dw==0 )
            continue;

          A3DPixPos_t evalme( u, v, w);
          neighborhood_points.emplace_back( std::move(evalme) );
        }
      }
    }

    for ( auto const& latticept : neighborhood_points ) {

      // turn this lattice point into a position in the images
      std::vector<float> tyz = lattice.getPos( latticept );
      float tick = tyz[0] + tick0;
      std::vector<int> wid(3,0);
      for ( int p=0; p<3; p++ ) {
        wid[p] = larlite::Geometry
      }

        // check if out of the image
        if ( r_neigh+min_r<0 || r_neigh+min_r>=(int)meta.rows() || c_neigh+min_c<0 || c_neigh+min_c>=(int)meta.cols() ) continue; 

        // is the neighbor within the start and end pad? then its forgiven if point is below charge
        // why do I do this?
        bool within_pad = false;
        if ( abs(r_neigh-start->row)<_config.astar_start_padding 
            && abs(c_neigh-start->col)<_config.astar_start_padding ) within_pad = true;
        if ( abs(r_neigh-goal->row)<_config.astar_end_padding 
            && abs(c_neigh-goal->col)<_config.astar_end_padding ) within_pad = true;

        // is this neighbor a bad ch?
        bool a_bad_ch = false;
        if ( use_bad_chs && m_badchimg->pixel( r_neigh+min_r, c_neigh+min_c )>0.5 ) a_bad_ch = true;

        // is this neighbor the goal? we need to know so we can always include it as a possible landing point, even if in badch region.
        bool isgoal = false;
        if ( c_neigh==goal->col && r_neigh==goal->row ) isgoal = true;
          
        // we skip this pixel if its below threshold. but we keep it if its a bad channel
        if ( !within_pad && !a_bad_ch && !isgoal
            && img.pixel( r_neigh+min_r, c_neigh+min_c )<_config.astar_threshold.at((int)meta.plane()) ) 
          continue; // skip if below threshold
          
        // do we already have the neighbor in the openset?
        PixPos_t neighpos( c_neigh, r_neigh );
        auto it_lookup = position_lookup.find( neighpos );

        AStar3DNode* neighbor_node = NULL;
        if ( it_lookup!=position_lookup.end() ) {
          if ((*it_lookup).second->closed==false ) {
            // if its on the openset we use it
            neighbor_node = (*it_lookup).second;
          }
        }
        else {
          // we need to make a new node object
          neighbor_node = new AStar3DNode( c_neigh, r_neigh, std::vector<float>(2,0.0) );
          // register with lookup table
          position_lookup.insert( pos_key_t( neighpos, neighbor_node ) );
          openset.addnode( neighbor_node );
        }

        // if we didn't find a node, move on
        if ( neighbor_node==NULL)  continue;

        // define the jump cost for this node
        // first get the diff vector from here to current node
        std::vector<float> dir3d(2,0.0);
        dir3d[0] = c_neigh - current->col;
        dir3d[1] = r_neigh - current->row;
        float norm = 0.;
        for (int v=0; v<2; v++) {
          norm += dir3d[v]*dir3d[v];
        }
        norm = sqrt(norm);
        for (int v=0; v<2; v++ ) dir3d[v] /= norm;          

        float cosine = 0;
        for (int v=0; v<2; v++ ) cosine += current->dir3d[v]*dir3d[v];
        neighbor_node->jump_cost = (1.0-cosine) + norm;

        // calculate hscore: heuristic distance from neighbor to goal
        float dcol_goal = goal->col-neighbor_node->col;
        float drow_goal = goal->row-neighbor_node->row;
        float hscore = sqrt( dcol_goal*dcol_goal + drow_goal*drow_goal );

        float sscore = sqrt( (goal->col-start->col)*(goal->col-start->col) + (goal->row-start->row)*(goal->row-start->row) );

        // calculate the potential g-score if we connect current to this neighboor
        float gscore = current->gscore + norm;//neighbor_node->jump_cost*0.0;

        float fscore = gscore + hscore;

        // is this gscore better than the current one (or is it unassigned?)
        if ( neighbor_node->fscore==0 || (neighbor_node->fscore>0 && neighbor_node->fscore>fscore )) {
          // update the neighbor to follow this path
          if ( verbose>2 )
          std::cout << "  updating neighbor to with better path: from f=" << neighbor_node->fscore << " g=" << neighbor_node->gscore 
                    << " to f=" << fscore << " g=" << gscore << std::endl;
          neighbor_node->prev = current;
          neighbor_node->fscore = fscore;
          neighbor_node->gscore = gscore;
          neighbor_node->dir3d = dir3d;
          number_updates++;
        }
        else {
          if ( verbose>2 )
            std::cout << "  this neighbor already on better path. current-f=" << neighbor_node->fscore << " < proposed-f=" << fscore << std::endl;
        }

      }
    }
    if ( verbose>1 )
      std::cout << "number of node updates: " << number_updates << std::endl;
  }

  /*
  void AStar3DAlgo::evaluateBadChNeighbors( AStar3DNode* current, const AStar3DNode* start, const AStar3DNode* goal,
    AStar3DNodePtrList& openset, AStar3DNodePtrList& closedset, 
    const int neighborhood_size, const int min_c, const int min_r, const int win_c, const int win_r, 
    const larcv::Image2D& img, const larcv::ImageMeta& meta, const bool use_bad_chs, std::map< PixPos_t, AStar3DNode* >& position_lookup) {

    // here we look at possible bad channels
    int r_current = current->row;
    int c_current = current->col;

    // we go back 3 nodes to get direction
    AStar3DNode* past = current;
    for (int inode=0; inode<3; inode++) {
      if ( past->prev!=NULL)
        past = current->prev;
    }
    std::vector<float> pastdir(2,0.0);
    pastdir[0] = current->col - past->col;
    pastdir[1] = current->row - past->row;
    float normpast = sqrt( pastdir[0]*pastdir[0] + pastdir[1]*pastdir[1] );
    if (normpast>0) {
      pastdir[0] /= normpast;
      pastdir[1] /= normpast;
    }

    int number_badch_updates = 0;
    int number_badch_nodes_created = 0;

    for (int dr=-neighborhood_size; dr<=neighborhood_size; dr++) {
      for (int dc=-neighborhood_size; dc<=neighborhood_size; dc++) {
        if ( dr==0 && dc==0 ) continue; // skip self
        int r_neigh = r_current+dr;
        int c_neigh = c_current+dc;

        // check if out of bounds inside the search region
        if ( r_neigh<0 || r_neigh>=win_r || c_neigh<0 || c_neigh>=win_c ) continue;
        // check if out of the image
        if ( r_neigh+min_r<0 || r_neigh+min_r>=(int)meta.rows() || c_neigh+min_c<0 || c_neigh+min_c>=(int)meta.cols() ) continue; 

        if ( m_badchimg->pixel(r_neigh+min_r,c_neigh+min_c)<0.5 ) continue;

        // ok this is a badch in the neighborhood
        int gapch = c_neigh+min_c; // back to global coordinates
        int dcol = ( dc<0 ) ? -1 : 1;

        // find the next non-bad channel
        bool foundgoodch = false;
        bool foundthegoal = false; // track that we found the goal, because we need to propose it
        while ( !foundgoodch ) {
          if ( gapch<0 || gapch>=(int)meta.cols() ) break;
          float badchval = m_badchimg->pixel(r_neigh+min_r,gapch);
          if ( badchval>0 && (gapch-min_c)==goal->col ) foundthegoal = true; // ran into the goal sitting in a badch
          if ( badchval==0 || foundthegoal ) {
            foundgoodch = true;
            break;
          }
          gapch += dcol;
        }

        if ( !foundgoodch && !foundthegoal ) continue;

        if ( verbose>1 )
          std::cout << "stepped into badch=" <<  meta.pos_x( c_neigh+min_c ) << " and jumped to goodch=" << meta.pos_x( gapch ) << std::endl;

        // we create nodes in this channel, within the window
        for (int ncol=0; ncol<2; ncol++) {
          // we check ncol more goodchs than just the one over the gap because often that one is not super reliable
	  for ( int rgap=0; rgap<(int)meta.rows(); rgap++ ) {
  	    if ( rgap-min_r<0 || rgap-min_r>=win_r ) continue;//outside the window
            if ( gapch+ncol*dcol<0 || gapch+ncol*dcol>=(int)meta.cols()) continue;
    	    if ( img.pixel( rgap, gapch+ncol*dcol )>_config.astar_threshold.at((int)meta.plane()) 
      	      || ( (rgap-min_r)==goal->row && (gapch+ncol*dcol-min_c)==goal->col ) ) {
              // either the pixel is above threshold, or is THE GOAL. if so, evaluate this node.

              int rgap_win = rgap-min_r;
              int cgap_win = gapch+ncol*dcol-min_c;

              AStar3DNode* gap_node = NULL;

              PixPos_t pos(cgap_win,rgap_win);
              auto it = position_lookup.find(pos);
              if ( it==position_lookup.end()) {
           	// make new node
           	gap_node = new AStar3DNode( cgap_win, rgap_win, std::vector<float>(2,0.0) );
           	position_lookup.insert( std::pair<PixPos_t,AStar3DNode*>( pos, gap_node ) );
            	openset.addnode( gap_node );
              	if ( verbose>1 )
                  std::cout << "created node (" << meta.pos_x( gapch+ncol*dcol ) << "," << meta.pos_y( rgap ) << ") "
                    << "local=(" << cgap_win << "," << rgap_win << ") from badch crossing." << std::endl;
              	  number_badch_nodes_created++;
            	}
            	else {
              	  gap_node = (*it).second;
            	}

            	if ( gap_node==NULL || gap_node->closed==true ) continue;

            	// evaluate this node
            	std::vector<float> dir2node(2,0.0);
            	dir2node[0] = gap_node->col-current->col;
            	dir2node[1] = gap_node->row-current->row;
            	float norm_dist2node = sqrt( dir2node[0]*dir2node[0] + dir2node[1]*dir2node[1] );
            	dir2node[0] /= norm_dist2node;
            	dir2node[1] /= norm_dist2node;

            	float cosine = dir2node[0]*pastdir[0] + dir2node[1]*pastdir[1];
            	if ( normpast==0 )
              	cosine = 0.0;

            	// we add the distance*(1-cosine)*0.5 score to the gscore
            	//float penalty = norm_dist2node*0.5*(1.0-cosine)*10.0;
	    	float penalty = norm_dist2node*norm_dist2node*0.5*(1.0-cosine);
            	//float penalty = norm_dist2node;
            	float gscore = current->gscore + norm_dist2node + penalty;
            	float hscore = sqrt( (goal->col-gap_node->col)*(goal->col-gap_node->col) + (goal->row-gap_node->row)*(goal->row-gap_node->row) );
            	float fscore = gscore + hscore;

            	if ( gap_node->fscore==0 || gap_node->fscore>fscore ) {
              	// we update this node
              	if ( verbose>1 )
                	std::cout << "updating node (" << meta.pos_x( gapch ) << "," << meta.pos_y( rgap ) << ") from badch crossing: " 
                  	<< " current-f=" << gap_node->fscore << " f-update=" << fscore << " gap-penalty=" << penalty << " cosine=" << cosine
                  	<< std::endl;
              	gap_node->fscore = fscore;
              	gap_node->gscore = gscore;
              	gap_node->prev = current;    
              	gap_node->dir3d = dir2node;
              	number_badch_updates++;            
              }

            }
          }
        }
      }
    }
    if ( verbose>1 ) {
      std::cout << "number of bad channel updates: " << number_badch_updates << std::endl;
      std::cout << "number of bad channel nodes created: " << number_badch_nodes_created << std::endl;    
      if ( number_badch_nodes_created>0 && verbose>2)
        std::cin.get();
    }
  }
  */

  std::vector<AStar3DNode> AStar3DAlgo::makeRecoPath( AStar3DNode* start, AStar3DNode* goal, int origin_row, int origin_col, bool& path_completed ) {
                                                                                                
    path_completed = true;
    std::vector<AStar3DNode> path;
    AStar3DNode* current = goal;
    while ( *current!=*start ) {
      AStar3DNode translated( *current );
      translated.col += origin_col;
      translated.row += origin_row;
      path.emplace_back(std::move(translated));
      current = current->prev;
      if ( current==NULL )
        break;
    }
    if ( current==NULL || *current!=*start ) {
      path_completed = false;
      return path;
    }
    AStar3DNode startout( *start );
    startout.col += origin_col;
    startout.row += origin_row;
    path.emplace_back( std::move(startout) );
    return path;
  }

  larcv::Image2D AStar3DAlgo::visualizeScores( std::string score_name, const larcv::Image2D& orig_img, 
    const int min_c, const int min_r, const int win_c, const int win_r, 
    const std::map<PixPos_t,AStar3DNode*>  position_lookup) {

    // create the blank image
    const larcv::ImageMeta& orig_meta = orig_img.meta();
    // create the new meta
    double width = win_c*orig_meta.pixel_width();
    double height = win_r*orig_meta.pixel_height();
    double origin_x = orig_meta.pos_x( min_c );
    double origin_y = orig_meta.pos_y( min_r );
    larcv::ImageMeta meta( width, height, win_r, win_c, origin_x, origin_y, orig_img.meta().plane() );
    larcv::Image2D img( meta );
    img.paint(0.0);
    for ( int r=0; r<win_r; r++ ) {
      for (int c=0; c<win_c; c++ ) {
        PixPos_t pos(c,r);
        auto it = position_lookup.find( pos );
        if ( it!=position_lookup.end() ) {
          if ( score_name=="fscore" )
            img.set_pixel(r,c,(*it).second->fscore);
          else if ( score_name=="gscore" )
            img.set_pixel(r,c,(*it).second->gscore);
          else
            throw std::runtime_error("unrecognized score to visualize");
        }
      }
    }
    return img;
  }

  // =========================================================================================================
  // LATTICE METHODS

  AStar3DNode* Lattice::getNode( const A3DPixPos_t& nodeid ) {

    // check bounds
    for (int i=0; i<3; i++) {
      if ( nodeid[i]<0 || nodeid[i]>=m_widths[i]) return nullptr;
    }

    // search map
    auto it_node = find( mode );

    AStar3DNode* node = nullptr;

    if ( it_node==end() ) {
      // create a node
      node = new AStar3DNode( u, v, w, std::vector<float>(3,0) );
      insert( a3dpos_pair_t(key,node) );
    }
    else {
      node = (*it_node).second;
    }

    return node;
  }

  AStar3DNode* Lattice::getNode( const int u, const int v, const int w ) {

    // check for node in the map
    A3DPixPos_t nodeid(u,v,w);
    return getNode( nodeid );
  }



  A3DPixPos_t Lattice::getNodePos( std::vector<float>& pos ) {
    A3DPixPos_t nodeid(0,0,0);
    for (int i=0; i<3; i++) {
      nodeid[i] = (int) ( pos[i]-m_origin[i])/m_cm_per_pixel[i];
      if ( nodeid[i]<0 || nodeid[i]>=m_widths[i]) 
        return A3DPixPos_t(); // outside lattice, so send and empty one
    }
    return nodeid;
  }


  AStar3DNode* Lattice::getNode( const std::vector<float>& pos ) {
    return getNode( getNodePos(pos) );
  }

  std::vector<float> Lattice::getPos( const int u, const int v, const int w ) {
    A3DPixPos_t nodeid( u, v, w );
    return getPos( nodeid );
  }

  std::vector<float> Lattice::getPos( const A3DPixPos_t& nodeid ) {
    std::vector<float> pos(3,0.0); 
    for (int i=0; i<3; i++) {
      // check if within lattice?
      pos[i] = m_origin[i] + nodeid[i]*m_cm_per_pixel[i];
    }
    return pos;
  }
  
}

