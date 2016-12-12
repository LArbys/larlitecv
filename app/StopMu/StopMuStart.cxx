#include "StopMuStart.h"
#include <assert.h>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "UBWireTool/UBWireTool.h"
#include "dbscan/DBSCANAlgo.h"


namespace larlitecv {

  StopMuStart::StopMuStart() : verbose(0) {}
  StopMuStart::~StopMuStart() {}
  
  
  void StopMuStart::getStartDirection( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
				       const std::vector<larcv::Pixel2D>& start, const int rneighbor,const  int cneighbor, const float fThreshold,
				       std::array<float,3>& start_spacepoint, std::vector< std::array<float,2> >& start_dir2d, 
				       std::array<float,3>& start_dir3d  ) {
    
    float ftrigger_tick  = 3200.0;
    float step_size = 3.0;
    float max_time_step_cm = 1.0; // cm    
    
    const int nplanes = img_v.size();
    std::vector<dbscan::dbPoints> combo_points(nplanes);
    std::vector<dbscan::dbscanOutput> clustering_out;
    int endpoint_hitindex[nplanes];
      
    for (int p=0; p<nplanes; p++) {
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      const larcv::Pixel2D& endpoint = start.at(p);
      //std::cout << "plane=" << p << " start point: (" << endpoint.X() << "," << endpoint.Y() << ")" << std::endl;
      // gather pixels above threshold that are around endpoint
      for (int dr=-rneighbor;dr<=rneighbor; dr++) {
	int row = endpoint.Y() + dr;
	if ( row<0 || row >= meta.rows() ) continue;
	for ( int dc=-cneighbor;dc<=cneighbor; dc++ ) {
	  int col = endpoint.X() + dc;
	  if ( col<0 || col>= meta.cols() ) continue;
	  if ( img_v.at(p).pixel( row, col ) > fThreshold || (dc==0 && dr==0) ) {
	    // if pixel above threshold
	    std::vector<double> point(2);
	    point[0] = (double)col;
	    point[1] = (double)row;
	    combo_points[p].emplace_back( point );
	    if ( dc==0 && dr==0 ) {
	      endpoint_hitindex[p] = combo_points[p].size()-1;
	    }
	  }
	}
      }
      // cluster points
      dbscan::DBSCANAlgo algo;
      dbscan::dbscanOutput clout = algo.scan( combo_points[p], 5, 3 );
      //std::cout << "neighborhood on plane=" << p <<  " has " <<  clout.clusters.size() << " clusters" << std::endl;
      // store clustering results
      clustering_out.emplace_back( clout );
    }//end of plane loop
    

    //  (b) we get the cluster that the point is attached to
    int pointclusters[nplanes];
    for (int p=0; p<nplanes; p++) {
      pointclusters[p] = -1;
      for (int c=0; c<clustering_out.at(p).clusters.size(); c++) {
	if ( pointclusters[p]!=-1 ) break; // already found the cluster, move on
	const std::vector<int>& cluster = clustering_out.at(p).clusters.at(c);
	for (int ihit=0; ihit<cluster.size(); ihit++) {
	  int hitidx = cluster.at(ihit);
	  if ( hitidx==endpoint_hitindex[p] ) {
	    pointclusters[p] = c;
	    break;
	  }
	}
      }
      //std::cout << "end point for plane=" << p << " is on cluster=" << pointclusters[p]
      //	<< " with " << clustering_out.at(p).clusters.at(pointclusters[p]).size() << " hits" << std::endl;
    }//end of loop over planes
      
      
    //  (c) we find the furthest pixel and get the direction
    float plane_dir[nplanes][2];
    float maxdist[nplanes];
    for (int p=0; p<nplanes; p++) {
      maxdist[p] = 0.0;
      const larcv::ImageMeta& meta = img_v.at(p).meta();
      plane_dir[p][0] = plane_dir[p][1] = 0;
      const std::vector<int>& cluster = clustering_out.at(p).clusters.at(pointclusters[p]);
      float hitpos[2] = {0};
      for (int ihit=0; ihit<cluster.size(); ihit++) {
	float dx = (combo_points[p].at(ihit)[0]-start.at(p).X())*meta.pixel_width()*0.3; // cm
	float dy = -(combo_points[p].at(ihit)[1]-start.at(p).Y())*meta.pixel_height()*0.5*::larutil::LArProperties::GetME()->DriftVelocity(); // cm 
	// (note, above negative, because of reverse time order used in image 2D)
	float dist = sqrt( dx*dx + dy*dy );
	if ( dist>maxdist[p] ) {
	  maxdist[p] = dist;
	  hitpos[0] = combo_points[p].at(ihit)[0];
	  hitpos[1] = combo_points[p].at(ihit)[1];
	  plane_dir[p][0] = dx/dist;
	  plane_dir[p][1] = dy/dist;
	}
      }
      //std::cout << "plane=" << p << " cluster dir=(" << plane_dir[p][0] << "," << plane_dir[p][1] << ") maxdist=" << maxdist[p] << " cm" << std::endl;
      std::array<float,2> sdir;
      sdir[0] = plane_dir[p][0];
      sdir[1] = plane_dir[p][1];
      start_dir2d.emplace_back( sdir );
    }//end of loop over planes
      
      
    //  (d) we use the 2D directions to provide a space point
    // find initial distance to move one tick, use that as guide to find target tick, picking the shortest one.
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    float init_dist[nplanes];
    float target_tick[nplanes];
    float closest_tick = 1e6;
    for (int p=0; p<nplanes; p++) {
      if ( plane_dir[p][1]!=0.0 )
	init_dist[p] = fabs((step_size*img_v.at(p).meta().pixel_height()*cm_per_tick)/plane_dir[p][1]);// how many unit vectors required to go the step size
      else
	init_dist[p] = max_time_step_cm/cm_per_tick; 
      //std::cout << "init dist[" << p << "] " << init_dist[p] << " vs. " << max_time_step_cm/cm_per_tick << std::endl;
      if ( init_dist[p]>max_time_step_cm/cm_per_tick )
	init_dist[p] = max_time_step_cm/cm_per_tick;
      target_tick[p] = img_v.at(p).meta().pos_y( start.at(p).Y() ) + init_dist[p]*(plane_dir[p][1]/cm_per_tick);
      if ( target_tick[p]<closest_tick )
	closest_tick = target_tick[p];
      //std::cout << "initial direction determined by going from tick=" << img_v.at(p).meta().pos_y(start.at(p).Y()) << " --> " << closest_tick << std::endl;
    }// end of plane loop
      
    // get the wire coordinate of the target tick
    std::vector<int> target_col(nplanes,-1);
    std::vector<int> target_wires(nplanes,-1);
    std::vector<int> start_wires(nplanes,-1);
    for (int p=0; p<nplanes; p++) {
      if ( plane_dir[p][1]!=0.0 ) {
	float start_tick = img_v.at(p).meta().pos_y( start.at(p).Y() );
	//float s = (target_tick[p]-start_tick)*cm_per_tick/plane_dir[p][1];
	float s = (closest_tick-start_tick)*cm_per_tick/plane_dir[p][1];
	float dcm = plane_dir[p][0]*s; // change in cm's in wire direction after 's' units
	float dwire  = dcm/0.3; // [cm] / [cm/wire]
	float dcol   = dwire/img_v.at(p).meta().pixel_width();  // [dwire] / [dwire/pixel col]
	target_col[p] = (int)(start.at(p).X() + dcol);
      }
      else {
	target_col[p] = start.at(p).X();
      }
      start_wires[p]  = img_v.at(p).meta().pos_x( start.at(p).X() );
      target_wires[p] = img_v.at(p).meta().pos_x( target_col[p] );
    }
    //     std::cout << "initial column directions: " 
    // 	      << " u: " << start.at(0).X() << " --> " << target_col[0]
    // 	      << " v: " << start.at(1).X() << " --> " << target_col[1]
    // 	      << " y: " << start.at(2).X() << " --> " << target_col[2]
    // 	      << std::endl;
      
    // get the 3D position
    std::vector<float> spacepoint[2];
    double triarea[2] = {0};
    int crosses[2] = {0};
    // note these return (z,y,x), because, of course they do.
    larcv::UBWireTool::wireIntersection( start_wires,  spacepoint[0], triarea[0], crosses[0] );
    larcv::UBWireTool::wireIntersection( target_wires, spacepoint[1], triarea[1], crosses[1] );
      
    // use 3d points to get initial direction and position
    float start_tick = img_v.at(0).meta().pos_y( start.at(0).Y() );
    std::vector<float> init_pos(3,0.0);
    std::vector<float> end_pos(3,0.0);
    std::vector<float> init_dir(3,0.0);
    init_pos[0] = (start_tick-ftrigger_tick)*cm_per_tick;
    init_pos[1] = spacepoint[0][1];
    init_pos[2] = spacepoint[0][0];
    end_pos[0]  = (closest_tick-ftrigger_tick)*cm_per_tick;
    end_pos[1]  = spacepoint[1][1];
    end_pos[2]  = spacepoint[1][0];
    //std::cout << "start 3D point: (" << init_pos[0] << "," << init_pos[1] << "," << init_pos[2] << ") area=" << triarea[0] << "(" << crosses[0] << ")" << std::endl;
    //std::cout << "end 3D point:   (" << end_pos[0] << "," << end_pos[1] << "," << end_pos[2] << ") area=" << triarea[1] << "(" << crosses[1] << ")" << std::endl;
      
    for (int i=0; i<3; i++) {
      start_spacepoint[i] = init_pos[i];
    }
      
    float start_norm = 0.0;
    for (int i=0; i<3; i++) {
      init_dir[i] = end_pos[i]-init_pos[i];
      start_norm += init_dir[i]*init_dir[i];
    }
    start_norm = sqrt(start_norm);
    for (int i=0; i<3; i++) { 
      init_dir[i] /= start_norm;
      start_dir3d[i] = init_dir[i];
    }
      
  }//end of getStartDirection
    
  void StopMuStart::getStartDirectionV( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
					const std::vector<larcv::Pixel2D>& start, const int rneighbor, const int cneighbor, const float fThreshold,
					std::vector<float>& start_spacepoint, 
					std::vector< std::vector<float> >& start_dir2d, std::vector<float>& start_dir3d ) {
    // wraps the above function, because python bindings can't handle array yet
    Vec3D_t start_sp;
    PlaneVec2D_t start2d;
    Vec3D_t start3d;
    getStartDirection( img_v, badch_v, start, rneighbor, cneighbor, fThreshold, start_sp, start2d, start3d );
    start_spacepoint.resize(3,0.0);
    start_dir3d.resize(3,0.0);
    for (int i=0; i<3; i++) {
      start_spacepoint[i] = start_sp[i];
      start_dir3d[i] = start3d[i];
    }

    // translate
    for (int p=0; p<3; p++) {
      std::vector<float> dir2d(2,0.0);
      for (int i=0; i<2; i++) {
	dir2d[i] = start2d.at(p)[i];
      }
      start_dir2d.emplace_back( dir2d );
    }
  }

}
