#include "FlashMuonTaggerAlgo.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

namespace larlitecv {
  
  
  
  bool FlashMuonTaggerAlgo::findTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const larcv::Image2D& tpc_img,
					   std::vector< BoundaryEndPt >& trackendpts, larcv::Image2D& markedimg ) {
    
    markedimg = std::move( larcv::Image2D(tpc_img.meta()) );
    
    for ( auto& ptr_event_flash : opflashsets ) {
      for ( auto& opflash : *ptr_event_flash ) {
	
	const larcv::ImageMeta& meta = tpc_img.meta();
	int plane = meta.plane();
	
	float tick_target = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
	if ( fSearchMode==kCathode ) {
	  tick_target += fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick;
	}

	// check if we can search for this opflash
	if ( tick_target<meta.min_y() || tick_target>=meta.max_y() )
	  continue;
	
	int row_target = meta.row( tick_target );

	if ( fConfig.verbosity<1 ) {
	  std::cout << "============================================================================================" << std::endl;
	  std::cout << "[Opflash search]" << std::endl;
	  std::cout << "  opflash: p= " << plane 
		    << " tick_target=" << tick_target
		    << " row_target=" << row_target
		    << " drift_t=" << fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick << " ticks"
		    << std::endl;
	}
	
	// scan across wires, looking for a hit
	for (int iwire=0; iwire<meta.cols(); iwire++) {
	  if ( markedimg.pixel( row_target, iwire )<= 1.0 // not visited
	       && tpc_img.pixel( row_target, iwire )>fConfig.pixel_value_threshold[plane] ) {// above threshold  
	    
	    // cluster and return cluster that contains the query point
	    dbscan::dbscanOutput clout;
	    dbscan::dbPoints winpoints;
	    int containing_cluster = -1;
	    bool ok = getClusters( tpc_img, row_target, iwire, clout, winpoints, containing_cluster );
	    
	    if ( !ok ) continue;

	    // find extremities of cluster (in time)
	    int tmax = -1;
	    int tmin = -1;
	    int wmax = -1;
	    int wmin = -1;
	    for ( int ichit=0; ichit<clout.clusters.at(containing_cluster).size(); ichit++ ) {
	      int hitidx = clout.clusters.at(containing_cluster).at(ichit);
	      int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
	      int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
	      //markedimg.set_pixel( y_, x_, 10.0 );
	      if ( tmax==-1 || y_>tmax ) { tmax = y_; wmax = x_; };
	      if ( tmin==-1 || y_<tmin ) { tmin = y_; wmin = x_; };
	    }
	    if ( fConfig.verbosity<2 ) {
	      std::cout << "end points: max (r,c)=(" << tmax << ", " << wmax << ")"
			<< " tw=(" << meta.pos_y(tmax) << "," << meta.pos_x(wmax) << ")"
			<< "  min (r,c)=(" << tmin << "," << wmin << ")" 
			<< " tw=(" << meta.pos_y(tmin) << "," << meta.pos_x(wmin) << ")"
			<< "  versus: query_row=" << row_target
			<< std::endl;
	    }

	    // is this a marked flash-end?
	    // if extrema matches the annode flash hypothesis row. mark as interesting (Score 200)
	    bool success = false;
	    if ( abs(tmin-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) {
	      if ( fConfig.verbosity<1 ) std::cout << "MATCHED MIN END" << std::endl;
	      markedimg.set_pixel( tmin, wmin, 200.0 );
	      success = true;
	      BoundaryEndPt endpt;
	      endpt.t = tmin;
	      endpt.w = wmin;
	      trackendpts.emplace_back( std::move(endpt) );
	    }
	    else if ( abs(tmax-row_target)<=fConfig.endpoint_time_neighborhood.at(plane) ) { 
	      if ( fConfig.verbosity<1 ) std::cout << "MATCHED MAX_END" << std::endl;
	      success =true;
	      markedimg.set_pixel( tmax, wmax, 200.0 );
	      BoundaryEndPt endpt;
	      endpt.t = tmax;
	      endpt.w = wmax;
	      trackendpts.emplace_back( std::move(endpt) );
	    }

	    if ( success ) {
	      // mark good cluster, so we don't use it again
	      for ( int ichit=0; ichit<clout.clusters.at(containing_cluster).size(); ichit++ ) {
		int hitidx = clout.clusters.at(containing_cluster).at(ichit);
		int x_ = (int)winpoints.at(hitidx).at(0)+0.1;
		int y_ = (int)winpoints.at(hitidx).at(1)+0.1;
		if ( markedimg.pixel( y_, x_ )<100 )
		  markedimg.set_pixel( y_, x_, 10.0 );      
	      }	      
	    }

	  }//end of if point is interesting
	}//end of loop over wires
      }//end of loop over flashes
    }//end of loop over flash containers
    
  }//end of marked output
  
  // subalgos
  bool FlashMuonTaggerAlgo::getClusters( const larcv::Image2D& tpc_img, int query_row, int query_col, 
					 dbscan::dbscanOutput& cluster_info,  dbscan::dbPoints& winpoints, int& containing_cluster) {
    winpoints.clear();

    const larcv::ImageMeta& meta = tpc_img.meta();
    int plane = meta.plane();
    
    // new, unexplored region!
    // we define a window around this point: (query_row, c)
    float t1 = meta.pos_y(query_row);
    float w  = meta.pos_x(query_col);

    if ( fConfig.verbosity<2) 
      std::cout << "[flashmuontagger: clustering region, return cluster and cluster endpoints]" << std::endl;

    int centeridx = -1;
    for (int dwire=-fConfig.clustering_wire_neighborhood.at(plane); dwire<=fConfig.clustering_wire_neighborhood.at(plane); dwire++) {
      for (int dtwin=-fConfig.clustering_time_neighborhood.at(plane); dtwin<=fConfig.clustering_time_neighborhood.at(plane); dtwin++) {
	int r = query_row+dtwin;
	int c = query_col+dwire;
	// check validity
	if ( r>=0 && r<meta.rows() && c>=0 && c<meta.cols() && tpc_img.pixel(r,c)>fConfig.pixel_value_threshold[plane] ) {
	  // mark as visited
	  //stage1_annode_hits.at(plane).set_pixel( r, c, 10.0 );
	  std::vector<double> pt(2,0.0);
	  pt[0] = c;
	  pt[1] = r;
	  winpoints.emplace_back( pt );
	  if ( dtwin==0 && dwire==0 ) centeridx = winpoints.size()-1;
	}// if valid point
      }//end of time window
    }//end of wire window
    if ( fConfig.verbosity<1 ) {
      std::cout << "  exploring around: plane=" << plane
		<< " (c,r)=" << query_col << ", " << query_row << " (w,t)=(" << w << "," << t1 << ")  npoints=" << winpoints.size() << std::endl;
    }

    // clustering
    dbscan::DBSCANAlgo dbalgo;
    cluster_info = dbalgo.scan( winpoints, 3, 3.0, false, 0.0 );

    // which cluster is the center point in?
    containing_cluster = cluster_info.clusterid.at(centeridx);
    if ( fConfig.verbosity<2 ) {
      std::cout << "  connected clusterid=" << containing_cluster  << std::endl;
      for (int ic=0; ic<cluster_info.clusters.size(); ic++) {
	std::cout << "    clusterid=" << ic << ", size=" << cluster_info.clusters.at(ic).size() << std::endl;
      }
    }
		  
    // cluster that hit is a part of is too small
    if ( cluster_info.clusters.at(containing_cluster).size()<5 )
      return false;

    if ( fConfig.verbosity<1 )
      std::cout << "  valid connected cluster: " << containing_cluster << " cluster size=" << cluster_info.clusters.at(containing_cluster).size() << std::endl;
    
    return true;
    
      

  }
  
}

