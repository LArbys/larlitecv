#include "BoundaryMuonTaggerAlgo.h"
#include "dbscan/DBSCANAlgo.h"
#include <vector>

namespace larlitecv {

  int BoundaryMuonTaggerAlgo::searchforboundarypixels( const std::vector< larcv::Image2D >& imgs, std::vector< larcv::Image2D >& matchedpixels ) {

    if ( !_config.checkOK() ) 
      return kErr_NotConfigured;
    
    if ( imgs.size()<3 ) {
      return kErr_BadInput;
    }

    // clear the output container
    matchedpixels.clear();

    // get meta for input image
    const larcv::ImageMeta& meta = imgs.at(0).meta();

    // create Image2D objects to store the result of boundary matching
    // we use meta to make the same sized image as the input

    for (int i=0; i<12; i++) { // top, bottom, upstream, downstream x 3 planes
      larcv::Image2D matchimage( meta );
      matchimage.paint(0.0);
      matchedpixels.emplace_back( std::move(matchimage) );
    }

    // we need the wire downsampling factor
    int dsfactor = int( meta.pixel_width()+0.1 ); 

    // now loop over over the time of the images
    for (int r=0; r<meta.rows(); r++) {
      // loop through boundary type
      for (int b=0; b<4; b++) {
	// loop through combinations
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
	    for (int p=0; p<3; p++) {
	      for ( int ipixel=0; ipixel<(int)abovethresh[p].size(); ipixel++) {
		//matchedpixels.at( p*4 + b ).set_pixel( r, abovethresh[p].at(ipixel), imatch ); // give it the match index
		matchedpixels.at( p*4 + b ).set_pixel( r, abovethresh[p].at(ipixel), 256.0 ); // give it the match index
	      }
	    }
	  }
	}//end of match loop
      }//end of boundary loop
    }//end of row loop
    
    return kOK;
  }

  int BoundaryMuonTaggerAlgo::clusterBoundaryPixels( const std::vector< larcv::Image2D >& imgs, // original image
						     const std::vector< larcv::Image2D >& matchedpixels, // pixels consistent with boundary hits
						     std::vector< std::vector<BoundaryEndPt> >& end_points // clustered end points on each plane
						     ) {
    int chs = (int)matchedpixels.size(); // for each plane: top/bottom/upstream/downstream
    int nplanes = chs/4;
    end_points.resize(chs); // w,t
    for (int p=0; p<chs; p++) {
      // for each plane turn matched pixels into end_points via dbscan
      end_points.at(p).clear();

      const larcv::Image2D& img    = imgs.at(p/4);
      const larcv::Image2D& hitimg = matchedpixels.at(p);
      const larcv::ImageMeta& meta = matchedpixels.at(p).meta();

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
      dbscan::dbscanOutput clout = dbalgo.scan( hits, 5, 3.0, false, 0.0 );

      // now we get our points. we find the max and min in time
      // we choose the end where there is large differential in the charge seen on one side versus the other
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

	for (int dwire=-_config.edge_win_wires.at(p); dwire<_config.edge_win_wires.at(p); dwire++) {
	  for (int dtime=-_config.edge_win_times.at(p); dtime<_config.edge_win_times.at(p); dtime++) {
	    int tmin_r = tmin+dtime;
	    int tmin_w = wmin+dwire;
	    if ( tmin_r>=0 && tmin_r<meta.rows() && tmin_w>=0 && tmin_w<meta.cols() ){
	      float tmin_val = img.pixel( tmin_r, tmin_w );
	      if ( tmin_val>_config.edge_win_hitthresh.at(p) ) {
		if ( dtime>0 ) tminq_upedge += tmin_val;
		else if ( dtime<0 )  tminq_downedge += tmin_val;
	      }
	    }

	    int tmax_r = tmax+dtime;
	    int tmax_w = tmin+dwire;
	    if ( tmax_r>=0 && tmax_r<meta.rows() && tmax_w>=0 && tmax_w<meta.cols() ){
	      float tmax_val = img.pixel( tmax_r, tmax_w );
	      if ( tmax_val>_config.edge_win_hitthresh.at(p) ) {
		if ( dtime>0 ) tmaxq_upedge  += tmax_val;
		else if ( dtime<0 ) tmaxq_downedge += tmax_val;
	      }
	    }
	  }
	}// loop over upper and lower regions
	
	float tmin_diff = fabs( tminq_upedge-tminq_downedge );
	float tmax_diff = fabs( tmaxq_upedge-tmaxq_downedge );
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
	std::vector<BoundaryEndPt>& end_point_v =end_points.at(p);
	end_point_v.emplace_back( endpt );
      }//end of loop over clusters
    }//end of loop over planes

    return 0;
  }//end of clusterBoundaryPixels


}
