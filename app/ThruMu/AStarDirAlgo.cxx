#include "AStarDirAlgo.h"
#include <vector>
#include <cmath>

namespace larlitecv {

  
  std::vector<AStarDirNode> AStarDirAlgo::findpath( const larcv::Image2D& img, 
                                                    int start_row, int start_col, int goal_row, int goal_col, 
                                                    float thresh, bool use_bad_chs ) {
    
    if ( verbose<=1 )
      std::cout << "[[ASTAR DIR ALGO START: going from (c,r)=(" << start_col << "," << start_row << ") to (" << goal_col << "," << goal_row << ")" << std::endl;

    const larcv::ImageMeta& meta = img.meta();

    // first we define the limits of the search by making a bounding box
    int min_r = ( start_row<goal_row ) ? start_row : goal_row;
    int max_r = ( start_row>goal_row ) ? start_row : goal_row;
    int min_c = ( start_col<goal_col ) ? start_col : goal_col;
    int max_c = ( start_col>goal_col ) ? start_col : goal_col;
    float pixel_w = meta.pixel_width();
    float pixel_h = meta.pixel_height();
    
    // extend the the bounding box
    min_r = ( min_r-_config.image_padding>0 ) ? min_r - _config.image_padding : 0;
    max_r = ( max_r+_config.image_padding<=(int)meta.rows() ) ? max_r + _config.image_padding : meta.rows()-1;
    min_c = ( min_c-_config.image_padding>0 ) ? min_c - _config.image_padding : 0;
    max_c = ( max_c+_config.image_padding<=(int)meta.cols() ) ? max_c + _config.image_padding : meta.cols()-1;
    int win_r = std::abs( max_r-min_r )+1;
    int win_c = std::abs( max_c-min_c )+1;

    // we now make some definitions
    // index of a pixel with position (r,c) is: idx=wincol*r + c where r=row index, c=col index
    
    AStarDirNodePtrList openset;
    AStarDirNodePtrList closedset;
    std::map< PixPos_t, AStarDirNode* > position_lookup;

    // make the target node (window coorindates)
    AStarDirNode* goal = new AStarDirNode( goal_col-min_c, goal_row-min_r, std::vector<float>(2,0.0) );

    // make starting node (window coordinates)
    AStarDirNode* start = new AStarDirNode( start_col-min_c, start_row-min_r, std::vector<float>(2,0.0) );
    start->gscore = 0.0;
    start->fscore = sqrt(pixel_w*(start_row-goal_row)*pixel_w*(start_row-goal_row) + pixel_h*(start_col-goal_col)*pixel_h*(start_col-goal_col)); // starting heuristic

    // set the start direction to heads towards the goal
    start->dir2d.resize(2,0.0);
    start->dir2d[0] = goal->col-start->col;
    start->dir2d[1] = goal->row-start->row;
    float dist = sqrt( start->dir2d[0]*start->dir2d[0] + start->dir2d[1]*start->dir2d[1] );
    start->dir2d[0] /= dist;
    start->dir2d[1] /= dist;

    position_lookup.insert( pos_key_t( PixPos_t(start->col, start->row), start ) );

    openset.addnode( start );
    
    if ( verbose>0 ) {
      std::cout << "start astar algo." << std::endl;
      std::cout << " in window coordinates (c,r): start(" << start->col << "," << start->row << ") -> (" << goal->col << "," << goal->row << ")" << std::endl;
    }

    int neighborhood_size = _config.astar_neighborhood.at(meta.plane());
    int nsteps = -1;
    AStarDirNode* current = NULL;
    while ( openset.size()>0 ) {
      // get current
      current = openset.back();
      openset.pop_back();
      nsteps++;
      if ( verbose>2 || (verbose>0 && nsteps%100==0) ) {
        std::cout << "step=" << nsteps << ": get current node (c,r) (" << current->col << "," << current->row << ") "
                  << " f=" << current->fscore << " g=" << current->gscore << " h=" << current->fscore-current->gscore << ". "
                  << "number of remaining nodes in openset=" << openset.size() << std::endl;
      }

      if ( *current==*goal ) {
        if ( verbose>0 ) {
          std::cout << "astar goal reached. finished." << std::endl;
        }
        break;        
      }

      // fill the candidate neighbors list
      std::vector< AStarDirNode* > neighborhood_list;
      struct neighboor_compare_t {
        bool operator()( AStarDirNode* lhs, AStarDirNode* rhs ) {
          if ( lhs->jump_cost>rhs->jump_cost ) return true;
          return false;
        };
      } neighbor_compare;

      int r_current = current->row;
      int c_current = current->col;
      bool already_jumped_gap_forward = false; // marks when we have added neighbors across a bad ch gap, so we don't do it again
      bool already_jumped_gap_backward = false; // marks when we have added neighbors across a bad ch gap, so we don't do it again

      // scan through neighors, and ID the candidate successor node
      evaluateNeighborNodes( current, start, goal, openset, closedset, neighborhood_size, 
        min_c, min_r, win_c, win_r, img, meta, use_bad_chs, position_lookup );
      evaluateBadChNeighbors( current, start, goal, openset, closedset, 1, 
        min_c, min_r, win_c, win_r, img, meta, use_bad_chs, position_lookup );


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
    std::vector<AStarDirNode> path;
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
  
  void AStarDirAlgo::evaluateNeighborNodes( AStarDirNode* current, const AStarDirNode* start, const AStarDirNode* goal,
    AStarDirNodePtrList& openset, AStarDirNodePtrList& closedset, 
    const int neighborhood_size, const int min_c, const int min_r, const int win_c, const int win_r, 
    const larcv::Image2D& img, const larcv::ImageMeta& meta, const bool use_bad_chs, std::map< PixPos_t, AStarDirNode* >& position_lookup ) {

    int r_current = current->row;
    int c_current = current->col;

    int number_updates=0;
    for (int dr=-neighborhood_size; dr<=neighborhood_size; dr++) {
      for (int dc=-neighborhood_size; dc<=neighborhood_size; dc++) {
        if ( dr==0 && dc==0 ) continue; // skip self
        int r_neigh = r_current+dr;
        int c_neigh = c_current+dc;

        // check if out of bounds inside the search region
        if ( r_neigh<0 || r_neigh>=win_r || c_neigh<0 || c_neigh>=win_c ) continue;
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

        AStarDirNode* neighbor_node = NULL;
        if ( it_lookup!=position_lookup.end() ) {
          if ((*it_lookup).second->closed==false ) {
            // if its on the openset we use it
            neighbor_node = (*it_lookup).second;
          }
        }
        else {
          // we need to make a new node object
          neighbor_node = new AStarDirNode( c_neigh, r_neigh, std::vector<float>(2,0.0) );
          // register with lookup table
          position_lookup.insert( pos_key_t( neighpos, neighbor_node ) );
          openset.addnode( neighbor_node );
        }

        // if we didn't find a node, move on
        if ( neighbor_node==NULL)  continue;

        // define the jump cost for this node
        // first get the diff vector from here to current node
        std::vector<float> dir2d(2,0.0);
        dir2d[0] = c_neigh - current->col;
        dir2d[1] = r_neigh - current->row;
        float norm = 0.;
        for (int v=0; v<2; v++) {
          norm += dir2d[v]*dir2d[v];
        }
        norm = sqrt(norm);
        for (int v=0; v<2; v++ ) dir2d[v] /= norm;          

        float cosine = 0;
        for (int v=0; v<2; v++ ) cosine += current->dir2d[v]*dir2d[v];
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
          neighbor_node->dir2d = dir2d;
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

  void AStarDirAlgo::evaluateBadChNeighbors( AStarDirNode* current, const AStarDirNode* start, const AStarDirNode* goal,
    AStarDirNodePtrList& openset, AStarDirNodePtrList& closedset, 
    const int neighborhood_size, const int min_c, const int min_r, const int win_c, const int win_r, 
    const larcv::Image2D& img, const larcv::ImageMeta& meta, const bool use_bad_chs, std::map< PixPos_t, AStarDirNode* >& position_lookup) {

    // here we look at possible bad channels
    int r_current = current->row;
    int c_current = current->col;

    // we go back 3 nodes to get direction
    AStarDirNode* past = current;
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

              AStarDirNode* gap_node = NULL;

              PixPos_t pos(cgap_win,rgap_win);
              auto it = position_lookup.find(pos);
              if ( it==position_lookup.end()) {
                // make new node
                gap_node = new AStarDirNode( cgap_win, rgap_win, std::vector<float>(2,0.0) );
                position_lookup.insert( std::pair<PixPos_t,AStarDirNode*>( pos, gap_node ) );
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
                gap_node->dir2d = dir2node;
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

  std::vector<AStarDirNode> AStarDirAlgo::makeRecoPath( AStarDirNode* start, AStarDirNode* goal, int origin_row, int origin_col, bool& path_completed ) {
                                                                                                
    path_completed = true;
    std::vector<AStarDirNode> path;
    AStarDirNode* current = goal;
    while ( *current!=*start ) {
      AStarDirNode translated( *current );
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
    AStarDirNode startout( *start );
    startout.col += origin_col;
    startout.row += origin_row;
    path.emplace_back( std::move(startout) );
    return path;
  }

  larcv::Image2D AStarDirAlgo::visualizeScores( std::string score_name, const larcv::Image2D& orig_img, 
    const int min_c, const int min_r, const int win_c, const int win_r, 
    const std::map<PixPos_t,AStarDirNode*>  position_lookup) {

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
  
}

