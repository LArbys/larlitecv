#include "AStarGridAlgo.h"
#include <vector>
#include <cmath>

namespace larlitecv {

  
  std::vector<AStarNode> AStarGridAlgo::findpath( const larcv::Image2D& img, int start_row, int start_col, int goal_row, int goal_col, float thresh ) {
    
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
    max_r = ( max_r+10<meta.rows() ) ? max_r + 10 : meta.rows();
    min_c = ( min_c-10>0 ) ? min_c - 10 : 0;
    max_c = ( max_c+10<=meta.cols() ) ? max_c + 10 : meta.cols()-1;
    int win_r = std::abs( max_r-min_r );
    int win_c = std::abs( max_c-min_c );

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

    int nsteps = -1;
    while ( openset.nentries()>0 ) {
      // get current
      AStarNode current = openset.gettop();
      // move to closed se
      closedset.add( current );
      nsteps++;
      if ( verbose==0 || (verbose==1 && nsteps%100==0) ) {
 	std::cout << "step=" << nsteps << ": get current node (" << current.row << "," << current.col << ") "
 		  << " f=" << current.fscore << " g=" << current.gscore << "."
		  << "number of remaining nodes in openset=" << openset.nentries() << std::endl;
      }

      if ( current==goal ) // make reco path
	return makeRecoPath( start, current, camefrom, min_r, min_c );
      int r_current = current.row;
      int c_current = current.col;

      // scan through neighors, make open set
      for (int dr=-1; dr<=1; dr++) {
	for (int dc=-1; dc<=1; dc++) {
	  if ( dr==0 && dc==0 ) continue; // skip self
	  int r_neigh = r_current+dr;
	  int c_neigh = c_current+dc;
	  if ( r_neigh<0 || r_neigh>=win_r || c_neigh<0 || c_neigh>=win_c ) continue;
	  if ( r_neigh+min_r<0 || r_neigh+min_r>=meta.rows() || c_neigh+min_c<0 || c_neigh+min_c>=meta.cols() ) continue; // skip if outside the image of course
	  if ( img.pixel( r_neigh+min_r, c_neigh+min_c )<_config.astar_threshold.at((int)meta.plane()) ) continue; // skip if below threshold

	  // make the neighbor
	  AStarNode neighbor( c_neigh, r_neigh );
	  
	  if ( closedset.contains(neighbor) ) continue; // already searched through here
	  
	  float dw = dc*pixel_w;
	  float dh = dr*pixel_h;
	  float tentative_gscore = current.gscore + sqrt( std::fabs( dw*dw ) + std::fabs( dh*dh ) );

	  float goaldist = 0.;
	  goaldist += pixel_h*pixel_h*(neighbor.row+start_row - goal_row)*(neighbor.row+start_row - goal_row);
	  goaldist += pixel_w*pixel_w*(neighbor.col+start_col - goal_col)*(neighbor.col+start_col - goal_col);
	  goaldist = sqrt(goaldist);
	  float fscore = tentative_gscore + goaldist;
	  neighbor.gscore = tentative_gscore;
	  neighbor.fscore = fscore;
	  
	  if ( !openset.contains(neighbor) ) {
	    // not in the openset, add the neighbor
	    openset.add( neighbor );
	    if ( verbose==0 || (verbose==1 && nsteps%100==0) )
	      std::cout << "adding neighbor(" << neighbor.row << "," << neighbor.col << ") f=" << neighbor.fscore << std::endl;
	  }
	  else if ( tentative_gscore>=openset.get_gscore(neighbor) ) {
	    // else if exists, we check if this neighbor is a good move.
	    if ( verbose==0 || (verbose==1 && nsteps%100==0) )
	      std::cout << "   neighbor not a better path: tentative_gscore=" << tentative_gscore << " neigh=" << openset.get_gscore(neighbor) << std::endl;
	    continue; // not a better path, we can come back to it later
	  }
	  
	  // path is best until now
	  camefrom[neighbor] = current;
	  if ( verbose==0 || (verbose==1 && nsteps%100==0) ) {
	    std::cout << "appending neighbor(" << neighbor.row << "," << neighbor.col << ") f=" << neighbor.fscore << " to path: " << std::endl;
	    auto itcurrent = camefrom.find( neighbor );
	    if ( itcurrent!=camefrom.end() ) {
	      AStarNode prevnode( 0, 0 );
	      prevnode = itcurrent->second;
	      std::cout << "  from (" << prevnode.row << ", " << prevnode.col << ")" << std::endl;
	    }
	  }
	  
	}//end of dc loop
      }//end of dr loop
      
      if ( verbose==0 && nsteps%100==0 ) {
	std::cout << "[enter] to continue." << std::endl;
	std::cin.get();
      }
    }//end of while loop

    // return empty path
    std::vector<AStarNode> empty;
    return empty;
  }
  
  std::vector<AStarNode> AStarGridAlgo::makeRecoPath( const AStarNode& start, const AStarNode& goal, const std::map< AStarNode, AStarNode, node_pos_compare >& camefrom,
						      int origin_row, int origin_col ) {
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
	break;
      }
    }
    return path;
  }
  
}

