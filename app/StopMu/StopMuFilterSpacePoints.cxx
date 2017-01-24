#include "StopMuFilterSpacePoints.h"
#include <cmath>

#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv/app
#include "dbscan/DBSCANAlgo.h"

namespace larlitecv {

  StopMuFilterSpacePointsConfig MakeStopMuFilterSpacePointsConfig( larcv::PSet pset ) {
    StopMuFilterSpacePointsConfig cfg;
    stopmu_filter_cfg.filter_thrumu_by_tag = pset.get<bool>("FilterThruByTag");
    stopmu_filter_cfg.duplicate_radius_cm  = pset.get<float>("DuplicateRadiuscm");
    stopmu_filter_cfg.pixel_threshold      = pset.get<float>("PixelThreshold");
    stopmu_filter_cfg.track_width_pixels   = pset.get<int>("TrackLabelingPixelWidth");
    stopmu_filter_cfg.row_tag_neighborhood = pset.get<int>("TagRowNeighborhood");
    stopmu_filter_cfg.col_tag_neighborhood = pset.get<int>("TagColNeighbothood");
    return cfg;
  }
  
  StopMuFilterSpacePoints::StopMuFilterSpacePoints( const StopMuFilterSpacePointsConfig& config  ) {
    m_config = &config;
  }

  std::vector< std::vector< const larcv::Pixel2D* > > StopMuFilterSpacePoints::filterSpacePoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
												  const std::vector<larcv::Image2D>& thrumu_pixels,
												  const std::vector<larcv::Image2D>& badch_imgs) {
    // this function selects spacepoints we want to pass to the stop-mu tracker
    // the goal is to
    // (1) filter out space points that below any through-going muon points
    // (2) merge points that are on the same track (choosing the one closest to the boundary)

    std::cout << "StopMuFilterSpacePoints::StopMuFilterSpacePoints" << std::endl;


    // [ Remove Duplicates ]
    std::vector< std::vector<const larcv::Pixel2D*> > not_duplicagted;
    removeDuplicateEndPoints( spacepoints_list, thrumu_pixels, not_duplicagted );
    std::cout << "number of non-duplicated endpoints: " << not_duplicagted.size() << std::endl;

    // [ Remove Through-going Muons ]
    std::vector< std::vector<const larcv::Pixel2D*> > not_thrumu;
    if ( !m_config->filter_thrumu_by_tag ) {
      std::cout << "filter thrumu end points by proximinty to thrumu-tagged pixels" << std::endl;
      removeThroughGoingEndPointsFromPixVectors( not_duplicagted, thrumu_pixels, not_thrumu );
    }
    else {
      std::cout << "filter thrumu end points using tags" << std::endl;
      removeThroughGoingEndPointsFromTags( not_duplicagted, thrumu_pixels, not_thrumu );
    }
    
    std::cout << "number of non-duplicated, non-thru-mu endpoints: " << not_thrumu.size() << std::endl;

    return not_thrumu;
  }

  // ----------------------------------------------------------------------------------------------------
  // Through-going end points

  void StopMuFilterSpacePoints::removeThroughGoingEndPoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
							     const std::vector<larcv::Image2D>& thrumu_pixels,
							     std::vector< std::vector<const larcv::Pixel2D*> >& passlist ) {
    
    for ( auto &ev_ptr : spacepoints_list ) {
      size_t nendpts = (*ev_ptr).Pixel2DArray(0).size();
      for (size_t i=0; i<nendpts; i++) {
	// collect the pixels across the planes
	std::vector<const larcv::Pixel2D*> pix_v(3,0);
	for (int p=0; p<3; p++) {
	  pix_v[p] = &(ev_ptr->Pixel2DArray(p).at(i));
	}
	
	bool isthrumu = isEndPtNearThruMuTag( pix_v, thrumu_pixels );
	if ( !isthrumu )
	  passlist.push_back( pix_v );
	
      }//end of endpoint loop
    }//end of list loop
    
  }

  void StopMuFilterSpacePoints::removeThroughGoingEndPointsFromPixVectors( std::vector< std::vector<const larcv::Pixel2D*> >& spacepoints_list,
									   const std::vector<larcv::Image2D>& thrumu_pixels,
									   std::vector< std::vector<const larcv::Pixel2D*> >& passlist ) {
    
    for ( auto &pix_v : spacepoints_list ) {      
      bool isthrumu = isEndPtNearThruMuTag( pix_v, thrumu_pixels );
      if ( !isthrumu )
	passlist.push_back( pix_v );
    }//end of list loop

  }
  
  bool StopMuFilterSpacePoints::isEndPtNearThruMuTag( const std::vector<const larcv::Pixel2D*>& pix_v, const std::vector<larcv::Image2D>& thrumu_pixels ) {
    
    // we check the thrumu_pixel, is there a thru-mu tag nearby?
    int nplanes_tagged = 0;
    for (int p=0; p<3; p++) {
      int pix_row=pix_v[p]->Y();
      int pix_col=pix_v[p]->X();
      for (int r=-m_config->row_tag_neighborhood; r<=m_config->row_tag_neighborhood; r++) {
	for (int c=-m_config->col_tag_neighborhood; c<=m_config->col_tag_neighborhood; c++) {
	  int row = pix_row+r;
	  int col = pix_col+c;
	  if ( row<0 || row>=thrumu_pixels.at(p).meta().rows() || col<0 || col>=thrumu_pixels.at(p).meta().cols() )
	    continue;
	  if ( thrumu_pixels.at(p).pixel(row,col)>0 )
	    nplanes_tagged++;
	}
      }//end of r loop
    }//end of plane loop

    if ( nplanes_tagged>=2 )
      return true;
    return false;
  }//end of removeThroughGoingEndPoints

  void StopMuFilterSpacePoints::removeThroughGoingEndPointsFromTags( std::vector< std::vector<const larcv::Pixel2D*> >& spacepoints_list,
								     const std::vector<larcv::Image2D>& thrumu_pixels,
								     std::vector< std::vector<const larcv::Pixel2D*> >& passlist ) {
    for ( auto &pix_v : spacepoints_list ) {      
      if ( pix_v.at(0)->Width()<=0.1 )
	passlist.push_back( pix_v );
    }//end of list loop
    
  }

							     
  // ----------------------------------------------------------------------------------------------------

  void StopMuFilterSpacePoints::removeDuplicateEndPoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
							  const std::vector< larcv::Image2D >& img_v,
							  std::vector< std::vector<const larcv::Pixel2D*> >& outputlist ) {

    // replace these with parameters at some point
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    float cm_per_wire = 0.3;
    float trig_tick = 3200.0;
    int duplicate_cluster_row_border = 5;
    int duplicate_cluster_col_border = 5;
    float duplicate_cluster_maxpixdist = 3.0;
    float duplicate_cluster_threshold = 10.0;

    int nendpts_considered = 0;

    // are pixels of the same time next to one another?
    for ( auto &ev_ptr : spacepoints_list ) {
      size_t nendpts = (*ev_ptr).Pixel2DArray(0).size();
      std::vector<int> rejected(nendpts,0);

      for (size_t i=0; i<nendpts; i++) {
	
	nendpts_considered++;

	if ( rejected[i]==1 ) continue;

	// get the 3D points for i
	int tick_i = img_v.at(0).meta().pos_y( (*ev_ptr).Pixel2DArray(0).at(i).Y() ); // row -> tick
	std::vector<int> wids_i(3,0);
	for (size_t p=0; p<3; p++) {
	  wids_i[p] = img_v.at(p).meta().pos_x( (*ev_ptr).Pixel2DArray(p).at(i).X() ); // col -> wire ID
	}
	int crosses_i = 0;
	double triangle_area_i = -1;
	std::vector<float> intersection_i;
	larcv::UBWireTool::wireIntersection( wids_i, intersection_i, triangle_area_i, crosses_i );

	std::vector<float> pos3d_i(3,0);
	pos3d_i[0] = (tick_i-trig_tick)*cm_per_tick; // x position
	//if ( pos3d_i[0] < 0 || pos3d_i[0]>270.0 ) {
	//rejected[i]=1;
	//continue;
	//}
	pos3d_i[1] = intersection_i[1];
	pos3d_i[2] = intersection_i[0];
	
	
	for (size_t j=i+1; j<nendpts; j++) {
	  if ( rejected[j]==1 ) continue;

	  // get the 3D points for j
	  int tick_j = img_v.at(0).meta().pos_y( (*ev_ptr).Pixel2DArray(0).at(j).Y() ); // row -> tick
	  std::vector<int> wids_j(3,0);
	  for (size_t p=0; p<3; p++) {
	    wids_j[p] = img_v.at(p).meta().pos_x( (*ev_ptr).Pixel2DArray(p).at(j).X() ); // col -> wire ID
	  }
	  int crosses_j = 0;
	  double triangle_area_j = -1;
	  std::vector<float> intersection_j;
	  larcv::UBWireTool::wireIntersection( wids_j, intersection_j, triangle_area_j, crosses_j );

	  std::vector<float> pos3d_j(3,0);
	  pos3d_j[0] = (tick_j-trig_tick)*cm_per_tick; // x position
	  //if ( pos3d_j[0] < 0 || pos3d_j[0]>270.0 ) {
	  //rejected[j]=1;
	  //continue;
	  //}
	  pos3d_j[1] = intersection_j[1];
	  pos3d_j[2] = intersection_j[0];

	  // get the distance
	  float separation = 0.;
	  for (int v=0; v<3; v++) {
	    float dv = pos3d_i[v]-pos3d_j[v];
	    separation += dv*dv;
	  }
	  separation = sqrt(separation);

	  
	  // if close, are the pixels connected by a line of pixels
	  if ( separation<m_config->duplicate_radius_cm ) {

	    // are they on the same cluster on each plane?
	    int same_cluster_on_nplanes = 0;
	    for (size_t p=0; p<3; p++) {
	      dbscan::dbPoints data;
	      int row_i = (*ev_ptr).Pixel2DArray(p).at(i).Y();
	      int row_j = (*ev_ptr).Pixel2DArray(p).at(j).Y();
	      int col_i = (*ev_ptr).Pixel2DArray(p).at(i).X();
	      int col_j = (*ev_ptr).Pixel2DArray(p).at(j).X();
	      int row_min = (row_i < row_j ) ? row_i  : row_j;
	      int row_max = (row_i < row_j ) ? row_j  : row_i;
	      int col_min = (col_i < col_j ) ? col_i : col_j;
	      int col_max = (col_i < col_j ) ? col_i : col_j;
	      row_min = ( row_min>=duplicate_cluster_row_border ) ? row_min-duplicate_cluster_row_border : 0;
	      row_max = ( row_max+duplicate_cluster_row_border<img_v.at(p).meta().rows() ) ? row_max+duplicate_cluster_row_border : img_v.at(p).meta().rows()-1;
	      col_min = ( col_min>=duplicate_cluster_col_border ) ? col_min-duplicate_cluster_col_border : 0;
	      col_max = ( col_max+duplicate_cluster_col_border<img_v.at(p).meta().cols() ) ? col_max+duplicate_cluster_col_border : img_v.at(p).meta().cols()-1;
	      for (int r=row_min; r<=row_max; r++) {
		for (int c=col_min; c<=col_max; c++) {
		  if ( img_v.at(p).pixel(r,c)>duplicate_cluster_threshold ) {
		    std::vector<double> hitpix(2,0.0);
		    hitpix[0] = c;
		    hitpix[1] = r;
		    data.push_back( hitpix );
		  }
		}
	      }
	      dbscan::DBSCANAlgo dbalgo;
	      dbscan::dbscanOutput cluster = dbalgo.scan( data, 3, duplicate_cluster_maxpixdist );
	      std::vector<double> test_i(2,0.0);
	      std::vector<double> test_j(2,0.0);
	      test_i[0] = col_i;
	      test_i[1] = row_i;
	      test_j[0] = col_j;
	      test_j[1] = row_j;
	      int match_i = cluster.findMatchingCluster( test_i, data, duplicate_cluster_maxpixdist );
	      int match_j = cluster.findMatchingCluster( test_j, data, duplicate_cluster_maxpixdist );
	      if ( match_i==match_j )
		same_cluster_on_nplanes++;
	    }//end of loop over planes for clustering

	    std::cout << "clusters i=" << i << " and j=" << j << " are close in real space. dist=" << separation << " cm" 
		      << " and on the same cluster in " << same_cluster_on_nplanes << " planes." << std::endl;
	      
	    if ( same_cluster_on_nplanes>=2 ) {
	      // which has the smallest dwall
	      float dwall_i = dwall( pos3d_i );
	      float dwall_j = dwall( pos3d_j );
	      std::cout << " -- reject either cluster i=" << i << " and j=" << j << " dwalls are: i=" << dwall_i << " cm and j=" << dwall_j << std::endl;
	      if ( dwall_i < dwall_j )
		rejected[j] = 1;
	      else if ( dwall_j < dwall_j )
		rejected[i] = 1;
	      else
		rejected[j] = 1;
	    }
	  }//end of if end points are close in real space
	  
	}//end of j-loop
      }//end of i-loop
      
      // review the starting points, pass the good ones to the output container
      for (size_t i=0; i<nendpts; i++) {
	if ( rejected[i]==1 ) continue;
	std::vector< const larcv::Pixel2D* > passed;
	for (int p=0; p<3; p++)  {
	  const larcv::Pixel2D* ptr_pix2d = &( (*ev_ptr).Pixel2DArray(p).at(i) );
	  passed.push_back( ptr_pix2d );
	}
	outputlist.push_back( passed );
      }//end of spacepoint loop
    }//end of spacepoint list loop

    std::cout << "Filter saved " << outputlist.size() << " out of " << nendpts_considered << " end points." << std::endl;

  }//end of removeDuplicateEndPoints

  float StopMuFilterSpacePoints::dwall( std::vector<float>& pos3d ) {
    // z
    std::vector<float> dwall(3,1.0e9);
    dwall[0] = ( fabs(pos3d[0]) < fabs(pos3d[0]-260.0) ) ? fabs(pos3d[0]) : fabs(pos3d[0]-260.0); 
    dwall[1] = ( fabs(118.0-pos3d[1]) < fabs(-118.0-pos3d[1]) ) ? fabs(118.0-pos3d[1]) : fabs(-118.0-pos3d[1]);
    dwall[2] = ( fabs(pos3d[2]) < fabs(1050.0-pos3d[2]) ) ? fabs(pos3d[2]) : fabs(1050.0-pos3d[2]);
    
    size_t idx_closest = 0;
    float closest = dwall[0];
    for (size_t i=1; i<dwall.size(); i++) {
      if ( dwall[i]<closest ) {
	closest = dwall[i];
	idx_closest++;
      }
    }

    return dwall[idx_closest];
  }

}
