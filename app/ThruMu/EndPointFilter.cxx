#include "EndPointFilter.h"

namespace larlitecv {

	void EndPointFilter::filterEndPts( const std::vector< const BoundarySpacePoint* >& source, 
		const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {
		int nendpts = (int)source.size();

		// set the end points to pass
		endpoint_passes.resize( nendpts, 1 );

		// our job will be to filter them out
		removeBoundaryAndFlashDuplicates( source, img_v, badch_v, endpoint_passes);

	}

  void EndPointFilter::removeBoundaryAndFlashDuplicates( const std::vector< const BoundarySpacePoint* >& source, 
  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {
  	// simply look for boundary (top,bottom,upstream,downstream) and flash (anode,cathode,imageends) that are close and on the same cluster
  	// for anode/cathode remove them when near (top,bottom,upstream,downstream)
  	// for imageends, we favor these over top/bottom/upstream/downstream
    // replace these with parameters at some point

    // we expect the input vector to be in the following order according to type:
    // (top,bottom,upstream,downstream,anode,cathode,imageend)

  	// check the end points
  	//int type_index = 0;
  	//for ( auto endpt_v : source ) {
  	//}

    /*
    // old code borrowed from stopmufilter
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
    */
  }

  void EndPointFilter::removeSameBoundaryDuplicates( const std::vector< const BoundarySpacePoint* >& source, 
  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) { 
  	// for same boundary duplicates: those close to one another and on same cluster 	
  }

  void EndPointFilter::removeDiffBoundaryDuplicates( const std::vector< const BoundarySpacePoint* >& source, 
  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes ) {  	
  }

  bool EndPointFilter::areEndPointsNearbyAndOnSameCluster( const larlitecv::BoundaryEndPt& pta, const larlitecv::BoundaryEndPt& ptb, 
  	const larcv::Image2D& img, const larcv::Image2D& badch, const float radius_cm, dbscan::dbPoints& opt_cluster ) {
  	// static utility function
  	// we check if points are within radius in the projected (x,dwire) 2D-plane
  	// if they are, we form a window around the two points and cluster, checking if they are on the same cluster.
  	// if not, then we use AStar* to try and bridge cluster. 
  	// we return the cluster the end points sit on
  }


}
