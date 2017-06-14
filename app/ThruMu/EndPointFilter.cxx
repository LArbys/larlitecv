#include "EndPointFilter.h"

// larlite
#include "LArUtil/LArProperties.h"

namespace larlitecv {

	void EndPointFilter::filterEndPts( const std::vector< const BoundarySpacePoint* >& source,
		const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {
		int nendpts = (int)source.size();

		// set the end points to pass
		endpoint_passes.resize( nendpts, 1 );

		// our job will be to filter them out
		removeBoundaryAndFlashDuplicates( source, img_v, badch_v, endpoint_passes);
    removeSameBoundaryDuplicates( source, img_v, badch_v, endpoint_passes);
	}

  void EndPointFilter::removeBoundaryAndFlashDuplicates( const std::vector< const BoundarySpacePoint* >& source,
  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {
  	// simply look for boundary (top,bottom,upstream,downstream) and flash (anode,cathode,imageends) that are close and on the same cluster
  	// for anode/cathode remove them when near (top,bottom,upstream,downstream)
  	// for imageends, we favor these over top/bottom/upstream/downstream
    // replace these with parameters at some point

  	// set passing flag to 'yes' if the pass vector has not been setup
  	if ( endpoint_passes.size()!=source.size() )
	  	endpoint_passes.resize( source.size(), 1 );

  	// get the indices of flash end points and boundary end points
  	std::vector<int> flash_indices;
  	std::vector<int> boundary_indices;
  	for ( size_t idx=0; idx<source.size(); idx++ ) {
      if ( source.at(idx)==NULL || source.at(idx)->size()!=3 ) {
        std::stringstream ss;
        ss << __FILE__ << ":" << __LINE__
            << ": BoundarySpacePoints are not well formed."
            << " size()=" << source.at(idx)->size()
            << " idx=" << idx
            << std::endl;
        throw std::runtime_error(ss.str());
      }
  		if ( source.at(idx)->type()==larlitecv::kAnode || source.at(idx)->type()==larlitecv::kCathode ) {
  			flash_indices.push_back(idx);
  		}
  		else if ( source.at(idx)->type()!=larlitecv::kImageEnd )
  			boundary_indices.push_back(idx);
  	}

  	// now check combinations of boundary and flashes
  	int number_flashpts_filtered = 0;
  	for ( auto& boundary_idx : boundary_indices ) {
  		for ( auto& flash_idx : flash_indices ) {
  			if ( endpoint_passes.at(boundary_idx)==0 || endpoint_passes.at(flash_idx)==0 )
  				continue; // this has already been rejected

  			const BoundarySpacePoint& sp_boundary = *(source.at(boundary_idx));
  			const BoundarySpacePoint& sp_flash    = *(source.at(flash_idx));

        // data well formed?
        if ( source.at(boundary_idx)==NULL || source.at(flash_idx)==NULL || sp_boundary.size()!=3 && sp_flash.size()!=3 ) {
          std::stringstream ss;
          ss << __FILE__ << ":" << __LINE__
              << ": BoundarySpacePoints are not well formed."
              << " sp_boundary.size()=" << sp_boundary.size() << " sp_flash.size()=" << sp_flash.size()
              << " boundary_idx=" << boundary_idx << "/" << source.size()
              << " flash_idx=" << flash_idx << "/" << source.size()
              << std::endl;
          throw std::runtime_error(ss.str());
        }

  			// check if endpts are on same cluster
  			int planes_on_same_endpt = 0;
  			std::vector<dbscan::dbPoints> opt_clusters;
  			for (size_t p=0; p<sp_flash.size(); p++ ) {
  				const BoundaryEndPt& boundary_endpt = sp_boundary.at(p);
  				const BoundaryEndPt& flash_endpt    = sp_flash.at(p);
  				dbscan::dbPoints cluster_hits;
  				if ( areEndPointsNearbyAndOnSameCluster( boundary_endpt, flash_endpt, img_v.at(p), badch_v.at(p), 15.0, false, cluster_hits ) ) {
  					opt_clusters.emplace_back( std::move(cluster_hits) );
  					planes_on_same_endpt++;
  				}
  			}

  			if ( planes_on_same_endpt>=2 ) {
  				// we filter out the flash point
  				// we use the direction to get direction of each space point
  				// 3D dir
  				/*
  				std::vector<float> dir(3,0.0);
  				if ( fabs(sp_boundary.dwall()) < fabs(sp_flash.dwall()) ) {
  					// boundary is closer to the wall
  					dir[1] = sp_flash.pos()[1] - sp_boundary.pos()[1];
  					dir[2] = sp_flash.pos()[2] - sp_boundary.pos()[2];
  				}
  				else {
  					// flash is closer to the wall.
  					dir[1] = sp_boundary.pos()[1] - sp_flash.pos()[1];
  					dir[2] = sp_boundary.pos()[2] - sp_flash.pos()[2];
  				}
  				float norm = sqrt( dir[1]*dir[1] + dir[2]*dir[2] );
  				sp_boundary.dir( dir );
  				*/
  				endpoint_passes.at( flash_idx ) = 0;
  				number_flashpts_filtered++;
  			}

  		}
  	}
  	std::cout << "Number of Flash Endpoints filtered: " << number_flashpts_filtered << std::endl;
  }

  void EndPointFilter::removeSameBoundaryDuplicates( const std::vector< const BoundarySpacePoint* >& source,
  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {
  	// for same boundary duplicates: those close to one another and on same cluster
 		// set passing flag to 'yes' if the pass vector has not been setup
  	if ( endpoint_passes.size()!=source.size() )
	  	endpoint_passes.resize( source.size(), 1 );

  	// get the indices of flash end points and boundary end points
  	std::vector< std::vector<int> > boundary_indices( (int)larlitecv::kNumEndTypes );
  	for ( size_t idx=0; idx<source.size(); idx++ ) {
  		int endpt_idx = (int)source.at(idx)->type();
  		if ( endpt_idx>=0 && endpt_idx<(int)larlitecv::kNumEndTypes )
  			boundary_indices.at(endpt_idx).push_back(idx);
  	}

  	// now check combinations of same boundary types
  	int number_endpts_filtered = 0;
  	for ( auto& boundary_idx : boundary_indices ) {
  		for (int ia=0; ia<(int)boundary_idx.size(); ia++) {
  			int idx_a = boundary_idx.at(ia);
  			if ( endpoint_passes.at(idx_a)==0 ) continue; // already filtered out

  			for (int ib=ia+1; ib<(int)boundary_idx.size(); ib++ ) {
  				int idx_b = boundary_idx.at(ib);
  				if ( endpoint_passes.at(idx_b)==0 ) continue;

	  			const BoundarySpacePoint& sp_a = *(source.at(idx_a));
  				const BoundarySpacePoint& sp_b = *(source.at(idx_b));

	  			// check if endpts are on same cluster
	  			int planes_on_same_endpt = 0;
  				std::vector<dbscan::dbPoints> opt_clusters;
  				for (size_t p=0; p<sp_a.size(); p++ ) {
	  				const BoundaryEndPt& endpt_a = sp_a.at(p);
  					const BoundaryEndPt& endpt_b = sp_b.at(p);
	  				dbscan::dbPoints cluster_hits;
  					if ( areEndPointsNearbyAndOnSameCluster( endpt_a, endpt_b, img_v.at(p), badch_v.at(p), 5.0, false, cluster_hits ) ) {
	  					opt_clusters.emplace_back( std::move(cluster_hits) );
  						planes_on_same_endpt++;
  					}
  				}

  				if ( planes_on_same_endpt>=2 ) {
  					// we filter out the end point with the bigger dwall
  					if ( sp_a.dwall()>sp_b.dwall()) {
  						// sp A is furtherest from the wall
  						endpoint_passes.at(idx_a) = 0;
  						break; // break b-loop as a is excluded now
  					}
  					else {
  						endpoint_passes.at(idx_b) = 0;
  					}
	  				number_endpts_filtered++;
	  			}
  			}
  		}
  	}
  	std::cout << "Number of Same-type Endpoints filtered: " << number_endpts_filtered << std::endl;
  }

  void EndPointFilter::removeDiffBoundaryDuplicates( const std::vector< const BoundarySpacePoint* >& source,
  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {

  	// for different boundary duplicates: those close to one another and on same cluster
 		// set passing flag to 'yes' if the pass vector has not been setup
  	if ( endpoint_passes.size()!=source.size() )
	  	endpoint_passes.resize( source.size(), 1 );

  	// get the indices of flash end points and boundary end points
  	std::vector< std::vector<int> > boundary_indices( (int)larlitecv::kDownstream+1 );
  	for ( size_t idx=0; idx<source.size(); idx++ ) {
  		if ( endpoint_passes.at(idx)==0 ) continue; // don't waste time with this one
  		int endpt_idx = (int)source.at(idx)->type();
  		if ( endpt_idx>=0 && endpt_idx<(int)boundary_indices.size() )
  			boundary_indices.at(endpt_idx).push_back(idx);
  	}

  	// now check combinations of boundary and flashes
  	int number_endpts_filtered = 0;
  	for (int itypea=0; itypea<(int)boundary_indices.size();  itypea++ ) {
	  	for (int itypeb=itypea+1; itypeb<(int)boundary_indices.size();  itypeb++ ) {

	  		std::vector< int >& indices_a = boundary_indices.at(itypea);
	  		std::vector< int >& indices_b = boundary_indices.at(itypeb);

  			for (int ia=0; ia<(int)indices_a.size(); ia++) {
		  		int idx_a = indices_a.at(ia);

	  			for (int ib=0; ib<(int)indices_b.size(); ib++ ) {

	  				int idx_b = indices_b.at(ib);

		  			const BoundarySpacePoint& sp_a = *(source.at(idx_a));
  					const BoundarySpacePoint& sp_b = *(source.at(idx_b));

	  				// check if endpts are on same cluster
	  				int planes_on_same_endpt = 0;
  					std::vector<dbscan::dbPoints> opt_clusters;
  					for (size_t p=0; p<sp_a.size(); p++ ) {
	  					const BoundaryEndPt& endpt_a = sp_a.at(p);
  						const BoundaryEndPt& endpt_b = sp_b.at(p);
	  					dbscan::dbPoints cluster_hits;
  						if ( areEndPointsNearbyAndOnSameCluster( endpt_a, endpt_b, img_v.at(p), badch_v.at(p), 20.0, false, cluster_hits ) ) {
	  						opt_clusters.emplace_back( std::move(cluster_hits) );
  							planes_on_same_endpt++;
  						}
  					}
	  				if ( planes_on_same_endpt>=2 ) {
  						// we filter out the end point with the bigger dwall
  						if ( sp_a.dwall()>sp_b.dwall()) {
  							// sp A is furtherest from the wall
  							endpoint_passes.at(idx_a) = 0;
  							break; // break b-loop as a is excluded now
  						}
  						else {
  							endpoint_passes.at(idx_b) = 0;
  						}
	  					number_endpts_filtered++;
	  				}
  				}// end of b loop
  			}//end of a loop
  		}//end of container b loop
  	}//end of container a loop

  	std::cout << "Number of Different-type Endpoints filtered: " << number_endpts_filtered << std::endl;
  }

  bool EndPointFilter::areEndPointsNearbyAndOnSameCluster( const larlitecv::BoundaryEndPt& pta, const larlitecv::BoundaryEndPt& ptb,
  	const larcv::Image2D& img, const larcv::Image2D& badch, const float radius_cm, bool return_cluster, dbscan::dbPoints& opt_cluster ) {
  	// static utility function
  	// we check if points are within radius in the projected (x,dwire) 2D-plane
  	// if they are, we form a window around the two points and cluster, checking if they are on the same cluster.
  	// if not, then we use AStar* to try and bridge cluster.
  	// we return the cluster the end points sit on

  	const larcv::ImageMeta& meta = img.meta();
  	float dist = 0.;
  	float dcol = pta.col-ptb.col;
  	float drow = pta.row-ptb.row;
  	float cm_per_wire = 0.3;
  	float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
  	dcol *= meta.pixel_width()*cm_per_wire;
  	drow *= meta.pixel_height()*cm_per_tick;

  	dist = sqrt( dcol*dcol + drow*drow );

  	if ( dist>radius_cm ) {
  		return false;
  	}

  	int row_min = ( pta.row<ptb.row ) ? pta.row : ptb.row;
  	int row_max = ( pta.row>ptb.row ) ? pta.row : ptb.row;
  	int col_min = ( pta.col<ptb.col ) ? pta.col : ptb.col;
  	int col_max = ( pta.col>ptb.col ) ? pta.col : ptb.col;

  	// get hits within window
  	dbscan::dbPoints winhits;
  	dbscan::extractPointsFromImageBounds( img, 10.0, row_min-10, row_max+10, col_min-10, col_max+10);

  	dbscan::DBSCANAlgo algo;
  	dbscan::dbscanOutput cluster_output = algo.scan( winhits, 5, 3.0, false, 0.0 );
  	std::vector<double> pos_a(2,0.0);
  	pos_a[0] = pta.col;
  	pos_a[1] = pta.row;
  	std::vector<double> pos_b(2,0.0);
  	pos_b[0] = ptb.col;
  	pos_b[1] = ptb.row;

  	int clusterid_a = cluster_output.findMatchingCluster( pos_a, winhits, 1.0 );
  	int clusterid_b = cluster_output.findMatchingCluster( pos_b, winhits, 1.0 );

  	if ( clusterid_a==clusterid_b ) {
  		// same cluster
  		if ( return_cluster ) {
	  		opt_cluster.clear();
  			for ( int ihit=0; ihit<(int)cluster_output.clusters.at(clusterid_a).size(); ihit++ ) {
  				int hitidx = cluster_output.clusters.at(clusterid_a).at(ihit);
  				std::vector<double> hit = winhits.at(hitidx);
  				opt_cluster.emplace_back( std::move(hit) );
  			}
  		}
  		return true;
  	}
  	else {
  		std::cout << "nearby clusters (" << meta.pos_x(pta.col) << "," << meta.pos_y(pta.row) << ")" << " and"
	  		<< " (" << meta.pos_x(ptb.col) << "," << meta.pos_y(ptb.row) << ") were close (" << dist << " cm) but not on same cluster."
	  		<< std::endl;
  	}

  	// can we connect via dead channel jumping. (maybe later)

  	return false;
  }


}
