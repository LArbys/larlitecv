#include "BoundaryMuonTaggerAlgo.h"
#include "dbscan/DBSCANAlgo.h"
#include "AStarGridAlgo.h"
#include "BoundaryEndPt.h"
#include <vector>
#include <cmath>
#include <assert.h>
#include <ctime>
#include "BoundaryIntersectionAlgo.h"


namespace larlitecv {

  void BoundaryMuonTaggerAlgo::run() {
    if (true)
      return;
    
    // this is still a work in progress. keeping the order of operations here as notes for now
    // eventually this will be the main function run by the user
    // searchforboundarypixels
    // clusterBoundaryPixels
  }

  int BoundaryMuonTaggerAlgo::searchforboundarypixels( const std::vector< larcv::Image2D >& imgs, std::vector< larcv::Image2D >& matchedpixels ) {
    // this checks for pixels consistent with boundary crossings at the top/bottom/upstream/downstream portions of the TPC
    // returns an image that marks the location of such pixels
    // the vector returned has N*4 elements from N planes x 4 crossing types

    int ncrossings = 4;
    int nplanes = imgs.size();
    int nresults = ncrossings*nplanes;

    if ( !_config.checkOK() )  {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Not configured." << std::endl;
      return kErr_NotConfigured;
    }
    
    if ( imgs.size()<3 ) {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Expecting 3 planes currently." << std::endl;
      return kErr_BadInput;
    }

    // clear the output container
    matchedpixels.clear();

    // get meta for input image
    const larcv::ImageMeta& meta = imgs.at(0).meta();

    // create Image2D objects to store the result of boundary matching
    // we use meta to make the same sized image as the input

    for (int i=0; i<nresults; i++) { // top, bottom, upstream, downstream x 3 planes
      larcv::Image2D matchimage( meta );
      matchimage.paint(0.0);
      matchedpixels.emplace_back( std::move(matchimage) );
    }

    // we need the wire downsampling factor
    int dsfactor = int( meta.pixel_width()+0.1 ); 
    const clock_t begin_time = clock();
    
    // now loop over over the time of the images
    for (int r=0; r<meta.rows(); r++) {
      // loop through boundary type
      for (int b=0; b<4; b++) {
	// loop through combinations (dumb)
	for (int imatch=0; imatch<m_matches.nmatches( (larlitecv::BoundaryMatchArrays::Boundary_t)b ); imatch++) {
	  
	  int triple[3];
	  m_matches.getMatch( (larlitecv::BoundaryMatchArrays::Boundary_t)b, imatch, triple[0], triple[1], triple[2] );
	  

	  // now we have a match combo
	  // look if the image has charge above some threshold, within some neighborhood
	  
	  bool hascharge[3] = { false, false, false };
	  std::vector<int> abovethresh[3]; // we save col number of pixels above threshold in each plane with this

	  for (int p=0; p<3; p++) {
	    const larcv::Image2D& img = imgs.at(p);
	    int col = triple[p]/dsfactor;
	    for (int n=-_config.neighborhoods.at(p); n<=_config.neighborhoods.at(p); n++) {
	      if (col + n <0 || col + n>=meta.cols() )  continue; // skip out of bound columns
		   
	      if ( img.pixel( r, col + n ) > _config.thresholds.at(p) ) {
		hascharge[p] = true;
		abovethresh[p].push_back(col+n);
	      }
	    }
	    // small optimization, if we see that one plane does not have charge, the match will fail. move on.
	    if ( !hascharge[p] ) break;
	  }
	  
	  if ( hascharge[0] & hascharge[1] & hascharge[2] ) {
	    // match!  we write the match to the output
	    for (int p=0; p<nplanes; p++) {
	      for ( int ipixel=0; ipixel<(int)abovethresh[p].size(); ipixel++) {
		matchedpixels.at( p*ncrossings + b ).set_pixel( r, abovethresh[p].at(ipixel), 256.0 ); // 256 is arbitrary marker
	      }
	    }
	  }
	  
	}//end of match loop
      }//end of boundary loop
    }//end of row loop
    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "boundary pixel search took " << elapsed_secs << " secs" << std::endl;

    
    return kOK;
  }

  int BoundaryMuonTaggerAlgo::searchforboundarypixels3D( const std::vector< larcv::Image2D >& imgs, std::vector< larcv::Image2D >& matchedpixels ) {
    // this checks for pixels consistent with boundary crossings at the top/bottom/upstream/downstream portions of the TPC
    // returns an image that marks the location of such pixels
    // the vector returned has N*4 elements from N planes x 4 crossing types
    // the output images are in a different space:
    //  top/bottom: wire axis -> z-axis
    //  upstream/downstream: wire-axis -> y-axis

    int ncrossings = 4;
    int nplanes = imgs.size();
    int nresults = ncrossings*nplanes;

    if ( !_config.checkOK() )  {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Not configured." << std::endl;
      return kErr_NotConfigured;
    }
    
    if ( imgs.size()<3 ) {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Expecting 3 planes currently." << std::endl;
      return kErr_BadInput;
    }

    // clear the output container
    matchedpixels.clear();

    // get meta for input image
    const larcv::ImageMeta& meta = imgs.at(0).meta();

    // create Image2D objects to store the result of boundary matching
    // we use meta to make the same sized image as the input
    enum { top=0, bot, upstream, downstream };

    for (int i=0; i<nresults; i++) { // top, bottom, upstream, downstream x 3 planes
      larcv::Image2D matchimage( meta );
      matchimage.paint(0.0);
      matchedpixels.emplace_back( std::move(matchimage) );
    }

    // we need the wire downsampling factor
    int dsfactor = int( meta.pixel_width()+0.1 ); 
    const clock_t begin_time = clock();

    // now loop over over the time of the images
    for (int r=0; r<meta.rows(); r++) {
      // loop through boundary type
      for (int b=0; b<4; b++) {
	// loop through combinations (dumb)
	for (int imatch=0; imatch<m_matches.nmatches( (larlitecv::BoundaryMatchArrays::Boundary_t)b ); imatch++) {
	  
	  int triple[3];
	  m_matches.getMatch( (larlitecv::BoundaryMatchArrays::Boundary_t)b, imatch, triple[0], triple[1], triple[2] );

	  // now we have a match combo
	  // look if the image has charge above some threshold, within some neighborhood
	  
	  bool hascharge[3] = { false, false, false };
	  std::vector<int> abovethresh[3]; // we save col number of pixels above threshold in each plane with this

	  for (int p=0; p<3; p++) {
	    const larcv::Image2D& img = imgs.at(p);
	    int col = triple[p]/dsfactor;
	    for (int n=-_config.neighborhoods.at(p); n<=_config.neighborhoods.at(p); n++) {
	      if (col + n <0 || col + n>=meta.cols() )  continue; // skip out of bound columns
		   
	      if ( img.pixel( r, col + n ) > _config.thresholds.at(p) ) {
		hascharge[p] = true;
		abovethresh[p].push_back(col+n);
	      }
	    }
	    // small optimization, if we see that one plane does not have charge, the match will fail. move on.
	    if ( !hascharge[p] ) break;
	  }
	  
	  if ( hascharge[0] & hascharge[1] & hascharge[2] ) {

	    // match!  we write the match to the output
	    for (int p=0; p<nplanes; p++) {
	      
	      if ( ( b==upstream || b==downstream) && p==2 )
		continue; // no useful info here

	      for ( int ipixel=0; ipixel<(int)abovethresh[p].size(); ipixel++) {
		int pixcol = abovethresh[p].at(ipixel);
		int wid = (int)pixcol*meta.pixel_width();
		const std::vector<float>& start = m_WireData[p].wireStart.find(wid)->second;
		const std::vector<float>& end   = m_WireData[p].wireEnd.find(wid)->second;

		if ( b==top ) {
		  // we are filling out the z-position here
		  float prev_val = matchedpixels.at( b ).pixel( r, (int)end[2] );
		  matchedpixels.at( b ).set_pixel( r, int(end[2]), prev_val+50.0 );
		}
		else if ( b==bot ) {
		  // fill z
		  float prev_val = matchedpixels.at( b ).pixel( r, (int)start[2] );
		  matchedpixels.at( b ).set_pixel( r, int(start[2]), prev_val+50.0 );
		}
		else if ( b==upstream && p<2) {
		  // fill y
		  float y = start[1];
		  if ( p==1 )
		    y = end[1];
		  float prev_val = matchedpixels.at( b ).pixel( r, (int)(y+116.0) );
		  matchedpixels.at( b ).set_pixel( r, int(y + 116.0), prev_val+50.0 );
		}
		else if ( b==downstream && p<2 ) {
		  // fill y
		  float y = end[1];
		  if ( p==1 )
		    y = start[1];
		  float prev_val = matchedpixels.at( b ).pixel( r, (int)(y+116.0) );
		  matchedpixels.at( b ).set_pixel( r, int(y + 116.0), prev_val+50.0 );
		}
	      }//end of pixel loop
	    }//end number of planes
	  }//if match
	}//end of match loop
      }//end of boundary loop
    }//end of row loop
    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "boundary pixel search took " << elapsed_secs << " secs" << std::endl;

    
    return kOK;
  }


  int BoundaryMuonTaggerAlgo::clusterBoundaryPixels( const std::vector< larcv::Image2D >& imgs, // original image
						     const std::vector< larcv::Image2D >& matchedpixels, // pixels consistent with boundary hits
						     std::vector< std::vector<BoundaryEndPt> >& end_points // clustered end points on each plane
						     ) {
    // inputs:
    //   imgs: vector of original TPC images, assuming one per plane
    //   matchedpixels: result of 'searchforboundarypixels' function. N planes x 4 crossings long
    // outputs:
    //   end_points: list of (time,wire) points.  outer list is N planes x 
    
    const int ncrossings = 4;
    int chs = (int)matchedpixels.size(); // for each plane: top/bottom/upstream/downstream
    int nplanes = chs/ncrossings;
    end_points.resize(chs); // w,t

    // loop over planes
    for (int ich=0; ich<chs; ich++) {
      // for each plane turn matched pixels into end_points via dbscan
      end_points.at(ich).clear();
      
      int plane=ich/ncrossings;
      const larcv::Image2D& img    = imgs.at(plane);
      const larcv::Image2D& hitimg = matchedpixels.at(ich);
      const larcv::ImageMeta& meta = matchedpixels.at(ich).meta();

      // we make a hit list
      dbscan::dbPoints hits;
      for ( int r=0; r<meta.rows(); r++ ) {
	for (int c=0; c<meta.cols(); c++ ) {
	  if ( hitimg.pixel( r, c ) > 1.0 ) {
	    std::vector<double> pt(2,0.0);
	    pt.at(0) = c; // x
	    pt.at(1) = r; // y
	    hits.emplace_back( pt );
	  }
	}
      }
      if ( hits.size()==0 )
	continue;
      
      // we cluster our hits
      dbscan::DBSCANAlgo dbalgo;
      dbscan::dbscanOutput clout = dbalgo.scan( hits, _config.boundary_cluster_minpixels.at(plane), _config.boundary_cluster_radius.at(plane), false, 0.0 );

      // now we get our points. we find the max and min in time
      // we choose the end where there is larger differential in the charge seen on one side versus the other
      // update: this doesn't really work.  maybe 
      for (int ic=0; ic<clout.clusters.size(); ic++) {
	int tmax = -1;
	int tmin = -1;
	int wmax = -1;
	int wmin = -1;
	for (int ichit=0; ichit<clout.clusters.at(ic).size(); ichit++) {
	  int hitidx = clout.clusters.at(ic).at(ichit);
	  int x_ = (int)hits.at(hitidx).at(0)+0.1;
	  int y_ = (int)hits.at(hitidx).at(1)+0.1;
	  if ( tmax==-1 || y_>tmax ) { tmax = y_; wmax = x_; };
	  if ( tmin==-1 || y_<tmin ) { tmin = y_; wmin = x_; };
	}

	// calculat charge above and below the tmin/tmax
	float tminq_upedge = 0.0;
	float tminq_downedge = 0.0;
	float tmaxq_upedge = 0.0;
	float tmaxq_downedge = 0.0;

	for (int dwire=-_config.edge_win_wires.at(ich); dwire<_config.edge_win_wires.at(ich); dwire++) {
	  for (int dtime=-_config.edge_win_times.at(ich); dtime<_config.edge_win_times.at(ich); dtime++) {
	    int tmin_r = tmin+dtime;
	    int tmin_w = wmin+dwire;
	    if ( tmin_r>=0 && tmin_r<meta.rows() && tmin_w>=0 && tmin_w<meta.cols() ){
	      float tmin_val = img.pixel( tmin_r, tmin_w );
	      if ( tmin_val>_config.edge_win_hitthresh.at(ich) ) {
		if ( dtime>0 ) tminq_upedge += tmin_val;
		else if ( dtime<0 )  tminq_downedge += tmin_val;
	      }
	    }

	    int tmax_r = tmax+dtime;
	    int tmax_w = wmax+dwire;
	    if ( tmax_r>=0 && tmax_r<meta.rows() && tmax_w>=0 && tmax_w<meta.cols() ){
	      float tmax_val = img.pixel( tmax_r, tmax_w );
	      if ( tmax_val>_config.edge_win_hitthresh.at(ich) ) {
		if ( dtime>0 ) tmaxq_upedge  += tmax_val;
		else if ( dtime<0 ) tmaxq_downedge += tmax_val;
	      }
	    }
	  }
	}// loop over upper and lower regions
	
	float tmin_diff = fabs( tminq_upedge-tminq_downedge );
	float tmax_diff = fabs( tmaxq_upedge-tmaxq_downedge );
// 	std::cout << "tmin=" << tmin << ", tmax=" << tmax << ": "
// 		  << "tmin_diff: " << tmin_diff << " vs. tmax diff: " << tmax_diff 
// 		  << " tmin=|" << tminq_upedge << "-" << tminq_downedge << "|" << "  vs  "
// 		  << " tmax=|" << tmaxq_upedge << "-" << tmaxq_downedge << "|" << "  vs  " << std::endl;
	
	BoundaryEndPt endpt;
	if ( tmin_diff>tmax_diff ) {
	  // set end point to tmin
	  endpt.t = tmin;
	  endpt.w = wmin;
	}
	else {
	  endpt.t = tmax;
	  endpt.w = wmax;
	}
	// set endpoint type
	endpt.type = (BoundaryEndPt::BoundaryEnd_t)(ich%ncrossings);
	// pass the endpoint into the output container
	std::vector<BoundaryEndPt>& end_point_v = end_points.at(ich);
	end_point_v.emplace_back( endpt );
      }//end of loop over clusters
    }//end of loop over channel { ich: N planes x 4 crossings }

    return 0;
  }//end of clusterBoundaryPixels

  int BoundaryMuonTaggerAlgo::clusterBoundaryPixels3D( const std::vector< larcv::Image2D >& matchedpixels, // pixels consistent with boundary hits 
						       std::vector< std::vector<BoundaryEndPt> >& end_points // list of boundary end point triples
						       ) {
    // inputs:
    //   matchedpixels: result of 'searchforboundarypixels3D' function. 4 crossings long. position of hits in edge position and time
    // outputs:
    //   end_points: list of boundary vectors.  inner vector has 3 elements, one for each plane
    
    const int ncrossings = 4;
    int chs = (int)matchedpixels.size(); // for each plane: top/bottom/upstream/downstream

    // loop over crossing types
    enum { top=0, bot, upstream, downstream };
    for (int ich=0; ich<chs; ich++) {
      
      const larcv::Image2D& hitimg = matchedpixels.at(ich);
      const larcv::ImageMeta& meta = matchedpixels.at(ich).meta();
      
      // we make a hit list
      dbscan::dbPoints hits;
      for ( int r=0; r<meta.rows(); r++ ) {
	for (int c=0; c<meta.cols(); c++ ) {
	  if ( hitimg.pixel( r, c ) > 50.0 ) {
	    std::vector<double> pt(2,0.0);
	    pt.at(0) = c; // x
	    pt.at(1) = r; // y
	    hits.emplace_back( pt );
	  }
	}
      }
      if ( hits.size()==0 )
	continue;
      
      // we cluster our hits
      dbscan::DBSCANAlgo dbalgo;
      dbscan::dbscanOutput clout = dbalgo.scan( hits, 5, 6, false, 0.0 );
      
      // now we get our points. we find the weighted mean
      for (int ic=0; ic<clout.clusters.size(); ic++) {
	float qtot = 0.0;
	float t_ave = 0.0;
	float w_ave = 0.0;
	for (int ichit=0; ichit<clout.clusters.at(ic).size(); ichit++) {
	  int hitidx = clout.clusters.at(ic).at(ichit);
	  int c_ = (int)hits.at(hitidx).at(0)+0.1;
	  int r_ = (int)hits.at(hitidx).at(1)+0.1;
	  float w_ = hitimg.pixel( r_, c_ );
	  qtot += w_;
	  t_ave += r_*w_;
	  w_ave += c_*w_;
	}
	t_ave /= qtot;
	w_ave /= qtot;

	// now we have to figure out what wire this is!
	int y_wid = -1;
	int u_wid = -1;
	int v_wid = -1;
	if ( ich==top ) {
	  // value saved as z-position on top edge of detector
	  float z = w_ave; 
	  y_wid = (z - 0.25)/0.3;
	  u_wid = (z - 0.55)/0.6 + 0;
	  v_wid = (z - 0.26)/0.6 + 672;
	}
	else if ( ich==bot ) {
	  float z = w_ave;
	  y_wid = (z - 0.25)/0.3;
	  u_wid = (z - 0.26)/0.6 + 672;
	  v_wid = (z - 0.55)/0.6 + 0;
	}
	else if ( ich==upstream  ) {
	  float y = w_ave-116.0;
	  y_wid = 0;
	  u_wid = (117.15-y)/0.346413;
	  v_wid = (y - (-115.29))/0.346413;
	}
	else if ( ich==downstream ) {
	  float y = w_ave-116.0;
	  y_wid = 3455;
	  u_wid = (117.23-y)/0.346413 + 1728;
	  v_wid = (y - (-115.29))/0.346413 + 1728;
	}
	
	// make boundary points
	std::vector< BoundaryEndPt > endpt_v;
	BoundaryEndPt endpt_Y( t_ave, (int)(y_wid/meta.pixel_width()), (BoundaryEndPt::BoundaryEnd_t)ich );
	BoundaryEndPt endpt_U( t_ave, (int)(u_wid/meta.pixel_width()), (BoundaryEndPt::BoundaryEnd_t)ich );
	BoundaryEndPt endpt_V( t_ave, (int)(v_wid/meta.pixel_width()), (BoundaryEndPt::BoundaryEnd_t)ich );
	endpt_v.emplace_back( endpt_U );
	endpt_v.emplace_back( endpt_V );
	endpt_v.emplace_back( endpt_Y );
	
	end_points.emplace_back( endpt_v );
      }//end of loop over clusters
    }//end of loop over channel { 4 crossings }
    
    return 0;
  }//end of clusterBoundaryPixels


  
  
  int BoundaryMuonTaggerAlgo::makePlaneTrackCluster( const larcv::Image2D& img, const larcv::Image2D& badchimg,
						     const std::vector< BoundaryEndPt >& top, const std::vector< BoundaryEndPt >& bot,
						     const std::vector< BoundaryEndPt >& upstream, const std::vector< BoundaryEndPt >& downstream,
						     const std::vector< BoundaryEndPt >& anode, const std::vector< BoundaryEndPt >& cathode,
						     const std::vector< BoundaryEndPt >& imgends,
						     std::vector< larlitecv::BMTrackCluster2D >& trackclusters ) {
    // inputs:
    //   img: TPC image
    //   badchimg: image where bad channels are marked (not used just yet)
    //   top,bot, upstream, downstream, anode, cathode, imgends: these are all lists of endpoints found by various routines
    //     top, bot, upstream, downstream come from BoundaryMuonTaggerAlgo::clusterBoundaryPixels
    //     anode, cathode come from FlashMuonTaggerAlgo::findTrackEnds
    //     imgends come from FlashMuonTaggerAlgo::findImageBoundaryEnds
    // output:
    //   container of candidate 2D tracks
    
    // wrap references into a struct so i can treat them more like an array
    // the order matches the order in BoundaryEnd_t enum
    struct endpt_s {
      const std::vector< BoundaryEndPt >& top;
      const std::vector< BoundaryEndPt >& bot;
      const std::vector< BoundaryEndPt >& upstream;
      const std::vector< BoundaryEndPt >& downstream;
      const std::vector< BoundaryEndPt >& anode;
      const std::vector< BoundaryEndPt >& cathode;
      const std::vector< BoundaryEndPt >& imgends;
      const std::vector< BoundaryEndPt >& operator[](int i) {
	if (i==0) return top;
	else if (i==1) return bot;
	else if (i==2) return upstream;
	else if (i==3) return downstream;
	else if (i==4) return anode;
	else if (i==5) return cathode;
	else if (i==6) return imgends;
      };
    };
    endpt_s endpts = { top, bot, upstream, downstream, anode, cathode, imgends };
    int nendpts = 7;

    // instantiate astar algo
    larlitecv::AStarAlgoConfig astar_config;
    astar_config.astar_threshold    = _config.astar_thresholds;
    astar_config.astar_neighborhood = _config.astar_neighborhood;
    larlitecv::AStarGridAlgo algo( astar_config );
    algo.setVerbose(2);

    // pair up containers
    int ntotsearched = 0;
    // poor-mans profiling
    const clock_t begin_time = clock();
    for (int i=0; i<nendpts; i++) {
      for (int j=i+1; j<nendpts; j++) {
	const std::vector< BoundaryEndPt >& pts_a = endpts[i];
	const std::vector< BoundaryEndPt >& pts_b = endpts[j];
	// combinations from a and b
	int ncombinations = pts_a.size()*pts_b.size();
	std::cout << "endpoints " << i << "->" << j << ": search combinations=" << ncombinations << std::endl;
	int icombo = 1;
	for (int ia=0; ia<(int)pts_a.size(); ia++) {
	  const BoundaryEndPt& pta = pts_a.at(ia);
	  for (int ib=0; ib<(int)pts_b.size(); ib++) {
	    const BoundaryEndPt& ptb = pts_b.at(ib);

	    //ok, we have two points, we do an A* star path search to between the two points
	    larlitecv::AStarNode start( pta.w, pta.t );
	    larlitecv::AStarNode goal( ptb.w, ptb.t );
	    // for debug
//  	    if ( j==6 && goal.row>=1006)
//  	      algo.setVerbose( 0 );
//  	    else
//  	      algo.setVerbose( 2 );
	    
	    if ( abs(start.row-goal.row)>785 ) // make this a parameter
	      continue; // physically impossible to connect

	    ntotsearched++;

	    //std::cout << "track path-finding: start=(" << start.row << "," << start.col << ") end=(" << goal.row << "," << goal.col << ")" << std::endl;
	    std::vector< larlitecv::AStarNode > path = algo.findpath( img, start.row, start.col, goal.row, goal.col, 20.0, false );
// 	    std::cout << "path returned: " << path.size() << " nodes. ";
// 	    if ( path.size()==0 )
// 	      std::cout << " FAILED" << std::endl;
// 	    else
// 	      std::cout << " COMPLETED" << std::endl;
// 	    std::cout << "[enter] to continue." << std::endl;
// 	    std::cout << "endpoint combination: " << icombo << " of " << ncombinations << std::endl;
	    //std::cin.get();

	    if ( path.size()>10 ) {
	      //larcv::Pixel2DCluster pixelpath;
	      BMTrackCluster2D track2d;
	      // set image ends
	      track2d.start.t = pta.t;
	      track2d.start.w = pta.w;
	      track2d.start.type = pta.type;
	      track2d.end.t = ptb.t;
	      track2d.end.w = ptb.w;
	      track2d.end.type = ptb.type;
	      track2d.pixelpath.clear();
	      for ( auto& node : path ) {
		larcv::Pixel2D pixel( node.col, node.row );
		pixel.Intensity( img.pixel( node.row, node.col ) );
		track2d.pixelpath += pixel; // uniary operator!
	      }
	      std::cout << "icombo=" << icombo << ". Storing Track with length=" << path.size() << std::endl;
	      trackclusters.emplace_back( track2d );
	    }// if successful path
	    icombo++;
	  }//loop over pt b
	}//loop over pt a
      }//endpoint container j
    }//end point container i
    
    
    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "total paths searched: " << ntotsearched << " in " << elapsed_secs << " secs" << std::endl;
    
    return 0;
  }

  int BoundaryMuonTaggerAlgo::markImageWithTrackClusters( const std::vector<larcv::Image2D>& imgs, const std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters,
							  std::vector<larcv::Image2D>& markedimgs ) {
    for ( auto &img : imgs ) {
      const larcv::ImageMeta& meta = img.meta();
      larcv::Image2D markedimg( meta );
      markedimg.paint(0.0);

      int plane = (int)meta.plane();
      
      const std::vector< larlitecv::BMTrackCluster2D >& planetracks = trackclusters.at(plane);

      for ( auto &track : planetracks ) {
	for ( auto &pixel : track.pixelpath ) {
	  int col = pixel.X();
	  int row = pixel.Y();

	  for ( int dc=-_config.astar_neighborhood.at(plane); dc<=_config.astar_neighborhood.at(plane); dc++ ) {
	    for ( int dr=-_config.astar_neighborhood.at(plane); dr<=_config.astar_neighborhood.at(plane); dr++ ) {
	      int r = row+dr;
	      int c = col+dc;
	      if ( r<0 || r>=meta.rows() ) continue;
	      if ( c<0 || c>=meta.cols() ) continue;
	      float val = img.pixel( row, col );
	      if ( val>_config.thresholds.at(plane) ) {
		markedimg.set_pixel( r, c, 100.0 );
	      }
	    }
	  }
	}//end of pixel list
      }//end of planetracks
      markedimgs.emplace_back( markedimg );
    }
    
    return 0;
  }
  /*
  int BoundaryMonTaggerAlgo::runAstar( const BoundaryEndPt& start, const BoudaryEndPt& end, const larcv::Image2D& img, int start_pad, int end_pad ) {
    // This wraps/interfaces with the AStar algorithm on our image. runs one start to end check.
    // inputs
    // ------
    //  start: starting boundary point
    //  end: ending boundary point
    //  img: image with charge deposition we will follow
    //  start_pad: radius around starting point where we will include all hits in case not start point directly on track
    //  end_pad:   radius around ending point where we will include all hits in case end point not directly on track

    
    larlitecv::AStarAlgoConfig astar_config;
    astar_config.astar_threshold    = _config.astar_thresholds;
    astar_config.astar_neighborhood = _config.astar_neighborhood;
    larlitecv::AStarGridAlgo algo( astar_config );
    algo.setVerbose(2);

    larlitecv::AStarNode start( start.w, start.t );
    larlitecv::AStarNode goal( end.w, end.t );
    
    
    
  }
  */

  /*
  void BoundaryMuonTaggerAlgo::filterShortAndKinkedTracks( int mintracklength, std::vector< larlitecv::BMTrackCluster2D >& tracks, 
							   std::vector< larlitecv::BMTrackCluster2D >& filtered, bool savefiltered ) {
    // this removes tracks with large kinks
    // input:
    //   mintracklength: filter out tracks shorter than this
    //   input: a list of bmtrackclusters.  we will remove elements
    //   output: filtered tracks. will be filled based on 'savefiltered' values
    //   savefiltered (optional): if true, passes filtered tracks into output container

    std::vector<int> remove_index_list;
    
    // loop over input tracks
    for ( auto &track : tracks ) {
      continue;
    }
    
  }
  */
  
  void BoundaryMuonTaggerAlgo::matchTracksStage1( const std::vector< larcv::Image2D >& imgs,
						  const std::vector< std::vector< larlitecv::BMTrackCluster2D >* >& plane2dtracks, 
						  std::vector< larlitecv::BMTrackCluster3D >& output  ) {
    // we do dumb O(N3) matching search.  
    // in the future we could structure start and end points into a tree for more efficient search: to make it Nlog(N)
    // also data products are starting to become heavy, might want to start using storage container and refer to indicices
    //     we'll see when things start to get bogged down.

    std::set< larlitecv::BMTrackCluster3D > track_set;

    int nplanes = (int)plane2dtracks.size();
    for (int p1=0; p1<1; p1++) {
      for (int p2=p1+1; p2<2; p2++) {
	
	const std::vector< larlitecv::BMTrackCluster2D >& p1tracks = *(plane2dtracks.at(p1));
	const std::vector< larlitecv::BMTrackCluster2D >& p2tracks = *(plane2dtracks.at(p2));

	for ( int idx1=0; idx1<p1tracks.size(); idx1++) {
	  const larlitecv::BMTrackCluster2D& track1 = p1tracks.at(idx1);
	  BoundaryEndPt::BoundaryEnd_t endpt_type = track1.start.type;
	  std::vector< int > match_indices;
	  std::vector< float > tot_start_diff;
	  std::vector< float > tot_end_diff;
	  for ( int idx2=0; idx2<p2tracks.size(); idx2++) {
	    const larlitecv::BMTrackCluster2D& track2 = p2tracks.at(idx2);
	    float sdiff,ediff;
	    bool start2start;
	    if ( !doTracksMatch( track1, track2, sdiff, ediff, start2start ) ) 
	      continue; // move on
	    
	    // check the third plane. what is it?
	    int p3 = -1;
	    if ( p1==0 && p2==1 ) p3 = 2;
	    else if ( p1==0 && p2==2 ) p3 = 1;
	    else if ( p1==1 && p2==2 ) p3 = 0;
	    if ( p3==-1 ) {
	      std::cout << "[BoundaryMuonTagger::matchTracksStage2] Invalid plane combination." << std::endl;
	      assert( false );
	    }

	    // look for third plane match
	    const std::vector< larlitecv::BMTrackCluster2D >& p3tracks = *(plane2dtracks.at(p3));
	    for ( int idx3=0; idx3<(int)p3tracks.size(); idx3++ ) {
	      const larlitecv::BMTrackCluster2D& track3 = p3tracks.at(idx3);
	      
	      float p1sdiff, p2sdiff, p1ediff, p2ediff;
	      bool p1start2start, p2start2start;
	      
	      if ( doTracksMatch( track1, track3, p1sdiff, p1ediff, p1start2start ) && doTracksMatch( track2, track3, p2sdiff, p2ediff, p2start2start ) ) {
		
		std::cout << " match on planes=(" << p1 << "," << p2 << "," << p3 << ") tracks=(" << idx1 << "," << idx2 << "," << idx3 << ")" << std::endl;

		// wow, made it all this way
		larlitecv::BMTrackCluster3D track3d;
		//
		track3d.trackidx.resize(3);
		track3d.trackidx.at(0) = idx1;
		track3d.trackidx.at(1) = idx2;
		track3d.trackidx.at(2) = idx3;

		// start/end time is average of all three
		int t_start  = track1.start.t;
		t_start += ( (start2start) ?   track2.start.t : track2.end.t );
		t_start += ( (p1start2start) ? track3.start.t : track3.end.t );
		t_start /= 3;

		int t_end = track1.end.t;
		t_end += ( (start2start) ?   track2.end.t : track2.start.t );
		t_end += ( (p1start2start) ? track3.end.t : track3.start.t );
		t_end /= 3;
		
		track3d.row_start = t_start;
		track3d.row_end   = t_end;
		track3d.tick_start = imgs.at(0).meta().pos_y( track3d.row_start );
		track3d.tick_end   = imgs.at(0).meta().pos_y( track3d.row_end );
		
		// set type
		track3d.start_type = track1.start.type;
		track3d.end_type   = track1.end.type;

		// set start wires
		track3d.start_wire.resize(3);
		track3d.start_wire.at(0) = track1.start.w;
		track3d.start_wire.at(1) = ( (start2start) ?   track2.start.w : track2.end.w );
		track3d.start_wire.at(2) = ( (p1start2start) ? track3.start.w : track3.end.w );

		// set end wires
		track3d.end_wire.resize(3);
		track3d.end_wire.at(0) = track1.end.w;
		track3d.end_wire.at(1) = ( (start2start) ?   track2.end.w : track2.start.w );
		track3d.end_wire.at(2) = ( (p1start2start) ? track3.end.w : track3.start.w );
		
		// boundary intersection algo
		larlitecv::BoundaryIntersectionAlgo algo;
		algo.determine3Dpoint( track3d.start_wire, track3d.start3D, track3d.start_type );
		algo.determine3Dpoint( track3d.end_wire, track3d.end3D, track3d.end_type );
		
		track_set.insert( track3d );
		output.emplace_back( track3d );
		
	      }
	    }
	    
	  }//end of loop over plane 2
	}
      }//end of loop over plane2
    }//end of loop over plane3

    std::cout << "Number of tracks in the set: " << track_set.size() << std::endl;

  }


  bool BoundaryMuonTaggerAlgo::doTracksMatch( const larlitecv::BMTrackCluster2D& track1, const larlitecv::BMTrackCluster2D& track2, 
					      float& start_t_diff, float& end_t_diff, bool& start2start ) {
    // inputs:
    //  track1: a TrackCluster2D
    //  track2: a TrackCluster2D
    // outputs:
    //  start_t_diff: difference in start point times
    //  end_t_diff: difference in end point times
    //  start2start: true if start point to start point matches. false if start to end points match.

    // does the first end point match in type?
    start2start = true; // 0=start; 1=end
    if ( track1.start.type == track2.start.type )
      start2start = true;
    else if ( track1.start.type!=track2.end.type )
      start2start = false;
    else
      return false; // does not match either
    
    // check if start_matches in time
    start_t_diff = 1e9;
    if ( start2start )
      start_t_diff = abs( track1.start.t - track2.start.t );	      
    else if ( start2start==1 )
      start_t_diff = abs( track1.start.t - track2.end.t );
    
    if ( start_t_diff > 10 )
      return false; // move on to the next one
    
    // does the end point match in type and time?
    end_t_diff = 1e9;
    bool end_matches_in_type = false;
    if ( start2start ) {
      end_t_diff = abs( track1.end.t - track2.end.t );
      if ( track1.end.type==track2.end.type ) end_matches_in_type = true;
    }
    else {
      end_t_diff = abs( track1.end.t - track2.start.t );
      if ( track1.end.type==track2.start.type ) end_matches_in_type = true;
    }
    
    if ( end_t_diff > 10 || !end_matches_in_type ) 
      return false;
    
    // profit!!!
    return true;
  }

  void BoundaryMuonTaggerAlgo::loadGeoInfo() {

    TFile fGeoFile( Form("%s/app/PMTWeights/dat/geoinfo.root",getenv("LARCV_BASEDIR")), "OPEN" );

    // Get the PMT Info
    fNPMTs = 32;
    TTree* fPMTTree  = (TTree*)fGeoFile.Get( "imagedivider/pmtInfo" );
    int femch;
    float pos[3];
    fPMTTree->SetBranchAddress( "femch", &femch );
    fPMTTree->SetBranchAddress( "pos", pos );
    for (int n=0; n<fNPMTs; n++) {
      fPMTTree->GetEntry(n);
      for (int i=0; i<3; i++) {
	pmtpos[femch][i] = pos[i];
      }
      //std::cout << "[POS " << femch << "] " << " (" << pmtpos[femch][0] << "," << pmtpos[femch][1] << "," << pmtpos[femch][2] << ")" << std::endl;
    }

    // Get the Wire Info
    TTree* fWireTree = (TTree*)fGeoFile.Get( "imagedivider/wireInfo" );
    int wireID;
    int planeID;
    float start[3];
    float end[3];
    fWireTree->SetBranchAddress( "wireID", &wireID );
    fWireTree->SetBranchAddress( "plane",  &planeID );
    fWireTree->SetBranchAddress( "start", start );
    fWireTree->SetBranchAddress( "end", end );
      
    int nentries = fWireTree->GetEntries();
    for ( int ientry=0; ientry<nentries; ientry++ ) {
      fWireTree->GetEntry(ientry);
      if ( m_WireData.find( planeID )==m_WireData.end() ) {
	// cannot find instance of wire data for plane id. make one.
	m_WireData[planeID] = larcv::pmtweights::WireData( planeID );
      }
      // start is always the one with the lowest y-value
      if ( start[1]>end[1] ) {
	// swap start end end
	for (int i=0; i<3; i++) {
	  float temp = start[i];
	  start[i] = end[i];
	  end[i] = temp;
	}
      }
      m_WireData[planeID].addWire( wireID, start, end );
    }
    
    fGeoFile.Close();

  }


  
}
