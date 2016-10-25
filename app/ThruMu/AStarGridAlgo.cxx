#include "AStarGridAlgo.h"
#include <vector>
#include <cmath>

namespace larlitecv {

  
  std::vector<AStarNode> AStarGridAlgo::findpath( const larcv::Image2D& img, 
						  int start_row, int start_col, int goal_row, int goal_col, 
						  float thresh, bool use_bad_chs ) {
    
    if ( verbose<=1 )
      std::cout << "[[ASTAR ALGO START: going from (c,r)=(" << start_col << "," << start_row << ") -> (" << goal_col << "," << goal_row << ")" << std::endl;

    const larcv::ImageMeta& meta = img.meta();

    // first we define the limits of the search by making a bounding box
    int min_r = ( start_row<goal_row ) ? start_row : goal_row;
    int max_r = ( start_row>goal_row ) ? start_row : goal_row;
    int min_c = ( start_col<goal_col ) ? start_col : goal_col;
    int max_c = ( start_col>goal_col ) ? start_col : goal_col;
    float pixel_w = meta.pixel_width();
    float pixel_h = meta.pixel_height();
    
    // extend the the bounding box
    min_r = ( min_r-10>0 ) ? min_r - 10 : 0;
    max_r = ( max_r+10<=meta.rows() ) ? max_r + 10 : meta.rows()-1;
    min_c = ( min_c-10>0 ) ? min_c - 10 : 0;
    max_c = ( max_c+10<=meta.cols() ) ? max_c + 10 : meta.cols()-1;
    int win_r = std::abs( max_r-min_r )+1;
    int win_c = std::abs( max_c-min_c )+1;

    // we now make some definitions
    // index of a pixel with position (r,c) is: idx=wincol*r + c where r=row index, c=col index
    
    AStarSet openset; // a stack of pixel indices (as defined above)
    AStarSet closedset;  // a set of pixel indices (as defined above)

    std::map< int, float > gscore;
    std::map< int, float > fscore;

    // make the target node (window coorindates)
    AStarNode goal( goal_col-min_c, goal_row-min_r );
    goal.gscore = 0.0;
    goal.fscore = 0.0;

    // make starting node (window coordinates)
    AStarNode start( start_col-min_c, start_row-min_r );
    start.gscore = 0.0;
    start.fscore = sqrt(pixel_w*(start_row-goal_row)*pixel_w*(start_row-goal_row) + pixel_h*(start_col-goal_col)*pixel_h*(start_col-goal_col)); // starting heuristic

    openset.add( start );
    
    std::map< AStarNode, AStarNode, node_pos_compare > camefrom;

    //std::cout << "start astar algo." << std::endl;
    //std::cout << " in window coordinates (c,r): start(" << start.col << "," << start.row << ") -> (" << goal.col << "," << goal.row << ")" << std::endl;
    int neighborhood_size = _config.astar_neighborhood.at(meta.plane());
    int nsteps = -1;
    while ( openset.nentries()>0 ) {
      // get current
      AStarNode current = openset.gettop();
      // move to closed se
      closedset.add( current );
      nsteps++;
      if ( verbose==0 || (verbose==1 && nsteps%100==0) ) {
 	std::cout << "step=" << nsteps << ": get current node (c,r) (" << current.col << "," << current.row << ") "
 		  << " f=" << current.fscore << " g=" << current.gscore << " h=" << current.fscore-current.gscore << ". "
		  << "number of remaining nodes in openset=" << openset.nentries() << std::endl;
      }

      if ( current==goal ) { // make reco path
	bool completed = false;
	std::vector< AStarNode > path = makeRecoPath( start, current, camefrom, min_r, min_c, completed );
	if (!completed)
	  path.clear();
	if ( verbose==0 ){
	  std::cout << "returning completed path" << std::endl;
	  std::cin.get();
	}
	return path;
      }
      int r_current = current.row;
      int c_current = current.col;
      bool already_jumped_gap_forward = false; // marks when we have added neighbors across a bad ch gap, so we don't do it again
      bool already_jumped_gap_backward = false; // marks when we have added neighbors across a bad ch gap, so we don't do it again

      // scan through neighors, make open set
      for (int dr=-neighborhood_size; dr<=neighborhood_size; dr++) {
	for (int dc=-neighborhood_size; dc<=neighborhood_size; dc++) {
	  if ( dr==0 && dc==0 ) continue; // skip self
	  int r_neigh = r_current+dr;
	  int c_neigh = c_current+dc;

	  // check if out of bounds inside the search region
	  if ( r_neigh<0 || r_neigh>=win_r || c_neigh<0 || c_neigh>=win_c ) continue;
	  // check if out of the image
	  if ( r_neigh+min_r<0 || r_neigh+min_r>=meta.rows() || c_neigh+min_c<0 || c_neigh+min_c>=meta.cols() ) continue; 

	  // is the neighbor within the start and end pad? then its forgiven if point is below charge
	  bool within_pad = false;
	  if ( abs(r_neigh-start.row)<_config.astar_start_padding 
	       && abs(c_neigh-start.col)<_config.astar_start_padding ) within_pad = true;
	  if ( abs(r_neigh-goal.row)<_config.astar_end_padding 
	       && abs(c_neigh-goal.col)<_config.astar_end_padding ) within_pad = true;


	  // is this neighbor a bad ch?
	  bool a_bad_ch = false;
	  if ( use_bad_chs && m_badchimg->pixel( r_neigh+min_r, c_neigh+min_c )>0.5 ) a_bad_ch = true;
	  
	  // we skip this pixel if its below threshold. but we keep it if its a bad channel
	  if ( !within_pad && !a_bad_ch && img.pixel( r_neigh+min_r, c_neigh+min_c )<_config.astar_threshold.at((int)meta.plane()) ) continue; // skip if below threshold
	  
	  // make the neighbor
	  AStarNode neighbor( c_neigh, r_neigh );

	  if ( closedset.contains(neighbor) ) continue; // already searched through here

	  // here we look to continue the path, what we do, however depends on if its a good or bad channel
	  if ( !a_bad_ch ) {
	    
	    //if ( closedset.contains(neighbor) ) continue; // already searched through here
	    
	    float dw = dc*pixel_w;
	    float dh = dr*pixel_h;
	    float tentative_gscore = current.gscore + sqrt( dw*dw + dh*dh );
	    
	    float goaldist = 0.;
	    goaldist += pixel_h*pixel_h*(neighbor.row - goal.row)*(neighbor.row - goal.row); // (must add start because nodes are rel.)
	    goaldist += pixel_w*pixel_w*(neighbor.col - goal.col)*(neighbor.col - goal.col); // (must add start because nodes are rel.)
	    goaldist = sqrt(goaldist);
	    float fscore = tentative_gscore + goaldist;
	    neighbor.gscore = tentative_gscore;
	    neighbor.fscore = fscore;
	    
	    if ( !openset.contains(neighbor) ) {
	      // not in the openset, add the neighbor
	      openset.add( neighbor );
	      //if ( verbose==0 || (verbose==1 && nsteps%100==0) )
	      //std::cout << "adding neighbor(" << neighbor.col << "," << neighbor.row << ") f=" << neighbor.fscore << std::endl;
	    }
	    else if ( tentative_gscore>=openset.get_gscore(neighbor) ) {
	      // else if exists, we check if this neighbor is a good move.
	      //if ( verbose==0 || (verbose==1 && nsteps%100==0) )
	      //std::cout << "   neighbor not a better path: tentative_gscore=" << tentative_gscore << " neigh=" << openset.get_gscore(neighbor) << std::endl;
	      continue; // not a better path, we can come back to it later
	    }
	  
	    // path is best until now
	    camefrom[neighbor] = current;
// 	    if ( verbose==0 || (verbose==1 && nsteps%100==0) ) {
// 	      std::cout << "appending neighbor(" << neighbor.col << "," << neighbor.row << ") f=" << neighbor.fscore << " to path: " << std::endl;
// 	      auto itcurrent = camefrom.find( neighbor );
// 	      if ( itcurrent!=camefrom.end() ) {
// 		AStarNode prevnode( 0, 0 );
// 		prevnode = itcurrent->second;
// 		std::cout << "  from (" << prevnode.row << ", " << prevnode.col << ")" << std::endl;
// 	      }
// 	    }
	  }//end of if a good channel
	  else if ( a_bad_ch ) {
	    // a bad channel

	    // get the next good channel. but first we need to know direction
	    float dir[2] = { 0.0, 0.0};
	    dir[0] = neighbor.col-current.col;
	    dir[1] = neighbor.row-current.row;
	   
	    int stepdir = 1;
	    if ( dir[0]<0 )
	      stepdir = -1;
	    if (dir[0]==0)
	      continue;

	    int cnext = neighbor.col + stepdir;
	    
	    while ( cnext+min_c>=0 && cnext+min_c<meta.cols() ) {
	      if ( m_badchimg->pixel( r_neigh+min_r, cnext+min_c )<1.0 )
		break;
	      cnext += stepdir;
	    }

	    if ( cnext+min_c<0 || cnext+min_c>=meta.cols() ) {
	      continue;
	    }

	    if ( verbose==0 || (verbose==1 && nsteps%100==0) ) {
	      std::cout << "hit neighbor that is a bad ch= current.col @ (c,r)=(" << current.col << "," << current.row << ")!" << std::endl;
	      std::cout << "  was going in direction=(" << dir[0] << "," << dir[1] << ")" << std::endl; 
	      std::cout << "  to jump gap, go to col=" << cnext << " in abs ch=" << cnext+min_c << std::endl;
	    }
	    if ( stepdir>0 )
	      already_jumped_gap_forward = true;
	    else
	      already_jumped_gap_backward = true;

	    // we load in neighbors along the jumping col (and one more in case of bad wires)
	    for (int ic=0; ic<2; ic++) {
	      int cskip = cnext + ic*stepdir;
	      if ( cskip+min_c<0 || cskip+min_c>=meta.cols() ) continue;
	      for (int jrow=0; jrow<win_r; jrow++ ) {
		if ( jrow+min_r<meta.rows() && img.pixel( jrow+min_r, cskip+min_c )>_config.astar_threshold.at((int)meta.plane()) ) {
		  
		  AStarNode jneighbor( cskip, jrow );
		  
		  if ( closedset.contains(jneighbor) ) continue; // already searched through here
		  
		  float dw = (jneighbor.col-current.col)*pixel_w;
		  float dh = (jneighbor.row-current.row)*pixel_h;
		  float tentative_gscore = current.gscore + sqrt( dw*dw + dh*dh );
		  float goaldist = 0.;
		  //goaldist += pixel_h*pixel_h*(jneighbor.row - goal_row)*(jneighbor.row - goal_row);
		  //goaldist += pixel_w*pixel_w*(jneighbor.col - goal_col)*(jneighbor.col - goal_col);
		  goaldist += pixel_h*pixel_h*(jneighbor.row - goal.row)*(jneighbor.row - goal.row);
		  goaldist += pixel_w*pixel_w*(jneighbor.col - goal.col)*(jneighbor.col - goal.col);
		  goaldist = sqrt(goaldist);
		  float fscore = tentative_gscore + goaldist;
		  jneighbor.gscore = tentative_gscore;
		  jneighbor.fscore = fscore;
		
		  if ( !openset.contains(jneighbor) ) {
		    // not in the openset, add the neighbor
		    openset.add( jneighbor );
		    //if ( verbose==0 || (verbose==1 && nsteps%100==0) )
		    //std::cout << "adding jump-gap neighbor(" << jneighbor.col << "," << jneighbor.row << ") f=" << jneighbor.fscore << " h(n)=" << goaldist << std::endl;
		  }
		  else if ( tentative_gscore>=openset.get_gscore(jneighbor) ) {
		    // else if exists, we check if move to this neighbor from our position is a good move. it's not.
		    //if ( verbose==0 || (verbose==1 && nsteps%100==0) )
		    //std::cout << "   jump-gap neighbor not a better path: tentative_gscore=" << tentative_gscore << " neigh=" << openset.get_gscore(jneighbor) << std::endl;
		    continue;
		  }
		
		  // path to jneighbor is best until now
		  camefrom[jneighbor] = current;
		  // 		if ( verbose==0 || (verbose==1 && nsteps%100==0) ) {
		  // 		  std::cout << "jumping gap to jneighbor(" << jneighbor.col << "," << jneighbor.row << ") f=" << jneighbor.fscore << " to path: " << std::endl;
		  // 		  auto itcurrent = camefrom.find( jneighbor );
		  // 		  if ( itcurrent!=camefrom.end() ) {
		  // 		    AStarNode prevnode( 0, 0 );
		  // 		    prevnode = itcurrent->second;
		  // 		    //std::cout << "  from (" << prevnode.col << ", " << prevnode.row << ")" << std::endl;
		  // 		  }
		  // 		}
		}//end of if gap crossing pixel is above threshold
	      }// end of loop over jumping rows
	    }//end of loop over jumping columns
	    //std::cin.get();
	  }//end of if a bad channel
	}//end of dc loop
      }//end of dr loop
      
      if ( verbose==0 && nsteps%100==0 ) {
	std::cout << "[enter] to continue." << std::endl;
	std::cin.get();
      }
    }//end of while loop

    if ( verbose<=1 ){
      std::cout << "returning UN-completed path" << std::endl;
      AStarNode closest = closedset.gettop();
      bool completed = false;
      std::vector< AStarNode > path = makeRecoPath( start, closest, camefrom, min_r, min_c, completed );
      if ( completed ) {
	std::cout << "returning closest path for diagnostics" << std::endl;
	return path;
      }
      std::cin.get();
    }

    // return empty path
    std::vector<AStarNode> empty;
    return empty;
  }
  
  std::vector<AStarNode> AStarGridAlgo::makeRecoPath( const AStarNode& start, const AStarNode& goal, const std::map< AStarNode, AStarNode, node_pos_compare >& camefrom,
						      int origin_row, int origin_col, bool& path_completed ) {
    path_completed = true;
    std::vector<AStarNode> path;
    AStarNode current(goal);
    while ( current!=start ) {
      AStarNode translated( current );
      translated.col += origin_col;
      translated.row += origin_row;
      path.push_back(translated);
      auto it = camefrom.find(current);
      if ( it!=camefrom.end() )
	current = it->second;
      else {
	std::cout << "MAKE RECO PATH ERROR: could not find origin of one of the nodes" << std::endl;
	path_completed = false;
	break;
      }
    }
    if ( path_completed ) {
      AStarNode translated(start);
      translated.col += origin_col;
      translated.row += origin_row;
      path.push_back(translated);
    }      
    return path;
  }
  
}

