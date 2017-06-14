#include <iostream>
#include <cmath>
#include <cstdlib>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventChStatus.h"
#include "dbscan/DBSCANAlgo.h"

// larlitecv
#include "ThruMu/EmptyChannelAlgo.h"

#include "KFStopMu.h"

int main( int nargs, char** argv ) {

  std::string cfg_file = "kfstop.cfg";

  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( "output_larcv_testextbnb_001.root", "larcv" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "KFStopMu" );
  dataco.initialize();

  dataco.goto_entry(0, "larcv");

  // PARAMETERS
  float fThreshold = 10.0;
  int rneighbor = 10;
  int cneighbor = 10;

  // for the test, we target a top-passing, stop muon
  // at time around tick 3880
  larcv::EventImage2D* imgs            = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
  larcv::EventPixel2D* top_spacepoints = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "topspacepts" );
  
  // make the bad channel image
  larlitecv::EmptyChannelAlgo emptyalgo;
  larcv::EventChStatus* ev_status      = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
  std::vector< larcv::Image2D > badchimgs = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status ); //< fill this out!
  std::cout << "number of bad ch imgs: " << badchimgs.size() << std::endl;
  
  const std::vector<larcv::Image2D>& img_v = imgs->Image2DArray();
  const larcv::ImageMeta& meta = img_v.at(0).meta();

  int nendpts = top_spacepoints->Pixel2DArray(0).size();
  std::vector<larcv::Pixel2D> start;
  for (int i=0; i<nendpts; i++) {
    const larcv::Pixel2D& pix = top_spacepoints->Pixel2DArray(0).at(i);
    float tick = meta.pos_y( pix.Y() );
    std::cout << "top endpoint #" << i << ": tick=" << tick << std::endl;
    if ( tick>3850 && tick<3950 ) {
      std::cout << "Found test start point:tick= " << tick << std::endl;
      for (int p=0; p<3; p++) {
	larcv::Pixel2D copy( top_spacepoints->Pixel2DArray(p).at(i) );
	start.emplace_back( copy );
      }
    }
  }
  
  std::cout << "start point: " << start.size() << std::endl;

  // ==========================================================================================
  // Kalman Filter for Stopping Muons
  
  // (1) determine initial direction
  //  (a) we cluster around the neighborhood of the point


  dbscan::dbPoints combo_points[3];
  std::vector<dbscan::dbscanOutput> clustering_out;
  int endpoint_hitindex[3] = {0};

  for (int p=0; p<3; p++) {
    const larcv::Pixel2D& endpoint = start.at(p);
    // gather pixels around endpoint
    for (int dr=-rneighbor;dr<=rneighbor; dr++) {
      int row = endpoint.Y() + dr;
      if ( row<0 || row >= img_v.at(p).meta().rows() ) continue;
      for ( int dc=-cneighbor;dc<=cneighbor; dc++ ) {
	int col = endpoint.X() + dc;
	if ( col<0 || col>= img_v.at(p).meta().cols() ) continue;
	if ( img_v.at(p).pixel( row, col ) > fThreshold ) {
	  // if pixel above threshold
	  std::vector<double> point(2);
	  point[0] = (double)col;
	  point[1] = (double)row;
	  combo_points[p].emplace_back( point );
	  if ( dc==0 && dr==0 )
	    endpoint_hitindex[p] = combo_points[p].size()-1;
	}
      }
    }
    // cluster points
    dbscan::DBSCANAlgo algo;
    dbscan::dbscanOutput clout = algo.scan( combo_points[p], 5, 3 );
    std::cout << "neighborhood on plane=" << p <<  ": " <<  clout.clusters.size() << " clusters" << std::endl;
    clustering_out.emplace_back( clout );
  }//end of plane loop

  //  (b) we get the cluster that the point is attached to
  int pointclusters[3] = {-1};
  for (int p=0; p<3; p++) {
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
    std::cout << "end point for plane=" << p << " is on cluster=" << pointclusters[p]
	      << " with " << clustering_out.at(p).clusters.at(pointclusters[p]).size() << " hits" << std::endl;
  }


  //  (c) we find the furthest pixel and get the direction
  float plane_dir[3][2];
  float maxdist[3] = {-1.0};
  for (int p=0; p<3; p++) {
    plane_dir[p][0] = plane_dir[p][1] = 0;
    const std::vector<int>& cluster = clustering_out.at(p).clusters.at(pointclusters[p]);
    float hitpos[2] = {0};
    for (int ihit=0; ihit<cluster.size(); ihit++) {
      float dx = (combo_points[p].at(ihit)[0]-start.at(p).X())*img_v.at(p).meta().pixel_width()*0.3; // cm
      float dy = -(combo_points[p].at(ihit)[1]-start.at(p).Y())*img_v.at(p).meta().pixel_height()*0.5*::larutil::LArProperties::GetME()->DriftVelocity(); // cm 
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
    std::cout << "plane=" << p << " cluster dir=(" << plane_dir[p][0] << "," << plane_dir[p][1] << ") maxdist=" << maxdist[p] << " cm" << std::endl;
  }


  //  (d) we use the 2D directions to provide a space point
  // find initial distance to move one tick, use that as guide to find target tick, picking the shortest one.
  float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
  std::cout << "cm per tick: " << cm_per_tick << std::endl;
  float init_dist[3] = {1.0};
  float target_tick[3] = {0.0};
  float max_time_step_cm = 1.0; // cm
  float closest_tick = 1e6;
  for (int p=0; p<3; p++) {
    if ( plane_dir[p][1]!=0.0 )
      init_dist[p] = fabs((3.0*img_v.at(p).meta().pixel_height()*cm_per_tick)/plane_dir[p][1]);// how many unit vectors required to go one tick
    else
      init_dist[p] = max_time_step_cm/cm_per_tick; 
    std::cout << "init dist[" << p << "] " << init_dist[p] << " vs. " << max_time_step_cm/cm_per_tick << std::endl;
    if ( init_dist[p]>max_time_step_cm/cm_per_tick )
      init_dist[p] = max_time_step_cm/cm_per_tick;
    target_tick[p] = img_v.at(p).meta().pos_y( start.at(p).Y() ) + init_dist[p]*(plane_dir[p][1]/cm_per_tick);
    if ( target_tick[p]<closest_tick )
      closest_tick = target_tick[p];
  }
  std::cout << "initial direction determined by going from tick=" << img_v.at(0).meta().pos_y(start.at(0).Y()) << " --> " << closest_tick << std::endl;

  // get the wire coordinate of the target tick
  int target_col[3] = { -1 };
  for (int p=0; p<3; p++) {
    if ( plane_dir[p][1]!=0.0 ) {
      float start_tick = img_v.at(p).meta().pos_y( start.at(p).Y() );
      float s = (target_tick[p]-start_tick)*cm_per_tick/plane_dir[p][1];
      float dcm = plane_dir[p][0]*s; // change in cm's in wire direction after 's' units
      float dwire  = dcm/0.3; // [cm] / [cm/wire]
      float dcol   = dwire/img_v.at(p).meta().pixel_width();  // [dwire] / [dwire/pixel col]
      target_col[p] = (int)(start.at(p).X() + dcol);
    }
    else {
      target_col[3] = start.at(p).X();
    }
  }
  std::cout << "initial wire directions: " 
	    << " u: " << start.at(0).X() << " --> " << target_col[0]
	    << " v: " << start.at(1).X() << " --> " << target_col[1]
	    << " y: " << start.at(2).X() << " --> " << target_col[2]
	    << std::endl;
  
  // get the 3D position
  
  
  // (1') determine initial direction: local PCA?

  // (2) start kalman filter

  // (3) kalman filter loop

  // (4) loop: make predicted measurement
  
  // (5) loop: search for true measurement by finding largest hit in the neighborhood

  // (6) loop: check stopping condition

  return 0;
}
