#include "StopMuTracker.h"

#include <stdio.h>
#include <time.h>
#include <sstream>

#ifndef __CINT__
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "StopMuSkeleton.h"

namespace larlitecv {

  StopMuTracker::StopMuTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v ) {
				
    
    // skeletonize image
    StopMuSkeleton skeleton_op;

    for (int p=0; p<3; p++) {
      larcv::Image2D skel = skeleton_op.skeletonize( img_v.at(p), 10.0, 3 );
      skel_v.emplace_back( std::move(skel) );
    }
    
    // cluster skeleton pixel, but mask thru-mu pixels
    time_t cluster_start = time(NULL);
    
    for (int p=0; p<3; p++) {
      dbscan::dbPoints data;
      const larcv::Image2D& skelimg = skel_v.at(p);
      for (int r=0; r<skelimg.meta().rows(); r++) {
	for (int c=0; c<skelimg.meta().cols(); c++) {
	  if ( skelimg.pixel(r,c)==0 ) continue; // not skeleton
	  if ( thrumu_v.at(p).pixel(r,c)>0 ) continue; // mask thrumu
	  std::vector<double> point(2);
	  point[0] = c; // X 
	  point[1] = r; // Y
	  data.emplace_back( point );
	}
      }
      dbscan::DBSCANAlgo dbalgo;
      dbscan::dbscanOutput cluster = dbalgo.scan( data, 5, 5.0 );
      imghits.emplace_back( data );
      clusters.emplace_back( cluster );
      std::cout << "number of clusters on plane " << p << ": " << clusters.at(p).clusters.size() << std::endl;
    }//loop over clusters
    
    time_t cluster_finished = time(NULL);
    //double dt_clusters = difftime(cluster_start,cluster_finished);
    double dt_clusters = cluster_finished-cluster_start;
    std::cout << "clustered in " << dt_clusters << " seconds." << std::endl;
    
    // initialize hit list
    for (int p=0; p<3; p++) {
      hitlists[p] = NULL;
      current_hit[p] = -1;
    }

  }

  void StopMuTracker::trackStopMu( const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
				   const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d, Step3D& trackstart ) {

    // parameters: move these to configuration file later
    float fCloseHitThreshold_cm = 1.0; // 3 wires
    float fStepSize_cm = 0.1;

    // find matching cluster
    std::vector<int> clusterid;

    const larcv::ImageMeta& meta = skel_v.at(0).meta();
    larcv::Image2D img_cluster( meta );
    img_cluster.paint(0.0);

    for (int p=0; p<3; p++) {
      std::vector<double> testpoint(2);
      testpoint[0] = start2d.at(p)[0]; // X
      testpoint[1] = start2d.at(p)[1]; // Y
      std::cout << "plane " << p << " testpoint: (" << testpoint[0] << "," << testpoint[1] << ")" << std::endl;
      int match = clusters.at(p).findMatchingCluster( testpoint, imghits.at(p), 5.0 );
      clusterid.push_back(match);
      std::cout << "plane=" << p << " matching cluster index=" << match << std::endl;

      // we make a sorted list of pixels by distance
      // we also need an initial direction
      if ( hitlists[p]!=NULL )
	delete hitlists[p];
      hitlists[p] = new Hit2DList;

      for (int ihit=0; ihit<clusters.at(p).clusters.at(match).size(); ihit++) {
	int hitidx = clusters.at(p).clusters.at(match).at(ihit);
	int x = imghits.at(p).at(hitidx)[0];
	int y = imghits.at(p).at(hitidx)[1];

	// dir from start to point
	std::vector<float> dir(2);
	dir[0] = x-start2d.at(p)[0];
	dir[1] = y-start2d.at(p)[1];
	float norm = sqrt( dir[0]*dir[0] + dir[1]*dir[1] );
	for (int i=0; i<2; i++) dir[i] /= norm;
	
	float cosine = 0.;
	for (int i=0; i<2; i++) cosine += dir[i]*start_dir2d.at(p)[i];

	Hit2D hit;
	hit[0] = x;
	hit[1] = y;
	if ( cosine>0 )
	  hit.distance = norm;
	else
	  hit.distance = -norm;
	if ( cosine>0 ) {
	  hitlists[p]->emplace(std::move(hit));
	  hitlists[p]->sort();
	}
	
	img_cluster.set_pixel((int)y,(int)x,250.0);
      }

      
      std::cout << "plane " << p << " number of hits=" << hitlists[p]->size() << std::endl;
      for (int ihit=0; ihit<(int)hitlists[p]->size(); ihit++) {
	std::cout << " [#" << ihit << "] (" << hitlists[p]->at(ihit)[0] << "," << hitlists[p]->at(ihit)[1] << ")"
		  << " " << hitlists[p]->at(ihit).distance << std::endl; 
      }
      
    }//end of loop over planes for sorting hits
    
    cv::Mat imgmat = larcv::as_mat( img_cluster );
    std::stringstream ss;
    ss << "baka.jpg";
    cv::imwrite( ss.str().c_str(), imgmat );

    int istep = 0;
    bool isfinished = false;

    // current step, direction
    std::vector< float > current_pos =  start_pos3d;
    std::vector< float > current_dir =  start_dir3d;
    std::vector< float > proposed_pos(3,0.0);
    std::vector< float > proposed_dir(3,0.0);
    _norm( current_dir );

    // make the first Step
    trackstart.pos = current_pos;
    trackstart.dir = current_dir;
    
    Step3D* current_step = &trackstart;

    while ( !isfinished && istep<100000 ) {
      
      std::cout << "===============================================================" << std::endl;
      std::cout << " Step " << istep << std::endl;

      // update 3D step
      makeProposedPos( current_pos, current_dir, proposed_pos, fStepSize_cm );

      std::cout << " from current pos=(" << current_pos[0] << "," << current_pos[1] << "," << current_pos[2] << ")"
		<< " and current dir=(" << current_dir[0] << "," << current_dir[1] << "," << current_dir[2] << ")" << std::endl;
      std::cout << " stepsize=" << fStepSize_cm << " cm" << std::endl;
      std::cout << " made proposal: " << proposed_pos[0] << "," << proposed_pos[1] << "," << proposed_pos[2] << ")" << std::endl;

      // calculate 2D positions
      int tick;
      std::vector<int> wid;
      imagePositions( proposed_pos, tick, wid );
      std::cout << "proposal location on the 2D planes: tick=" << tick << " planes U=" << wid[0] << " V=" << wid[1] << " Y=" << wid[2] << std::endl;

      std::vector<int> pixel_cols;
      int pixel_row;
      _wire2pixel( tick, wid, meta, pixel_cols, pixel_row );
      std::cout << "proposal location on the image: row=" << pixel_row << " plane col U=" << pixel_cols[0] << " V=" << pixel_cols[1] << " Y=" << pixel_cols[2] << std::endl;

      // find closest hits in cluster on each plane
      std::vector< std::vector< std::pair<int,double> > > closest_pixels(3);
      for (int p=0; p<3; p++) {
	std::vector<int> test_pos(2,0);
	test_pos[0] = pixel_cols[p]; // X
	test_pos[1] = tick;         // Y
	getClosestHitsInPlane( clusterid[p], test_pos, imghits.at(p), clusters.at(p), meta, closest_pixels.at(p) );
	std::cout << "plane " << p << " hitlist size: " << closest_pixels.at(p).size();
	if ( closest_pixels.at(p).size()> 0 )
	  std::cout << " closest dist=" << closest_pixels.at(p).at(0).second << " furthest=" << closest_pixels.at(p).back().second << std::endl;
	else
	  std::cout << std::endl;
      }

      
      // is there a close hit in each plane?
      std::cout << "closest distance to hits on each plane: ";
      bool closehits[3] = {false};
      int ngoodplanes = 0;
      for (int p=0; p<3; p++) {
	if ( closest_pixels.at(p).size()>0 && closest_pixels.at(p).at(0).second<fCloseHitThreshold_cm ) {
	  ngoodplanes++;
	  closehits[p] = true;
	  std::cout << " [p=" << p << ", " << closest_pixels.at(p).at(0).second << " cm] ";
	}
      }
      std::cout << std::endl;

      bool continue_search = true;
      if ( ngoodplanes>=2 ) {
	//  everything is fine continue
	current_pos = proposed_pos;
	std::cout << "number of good planes=" << ngoodplanes << ". continue tracker" << std::endl;
      }
      else {
	//  current path is no good. need to re-orient or stop
	std::cout << "need to re-direct track" << std::endl;
	continue_search = findNewDirection();
      }

      if ( !continue_search ) {
	std::cout << "Stopping condition met." << std::endl;
	break;
      }
      else {

	// is hit 3D consitent?
      
	// update tracker
	// prepare wire info
	std::vector<Point2D_t> closest_hits_list(3);
	std::vector<Point2D_t> plane_pos_list(3);
	for (int p=0; p<3; p++) {
	  std::vector<int> theclosesthit(2);
	  if ( closest_pixels.at(p).size()>0 ) {
	    int hitidx = closest_pixels.at(p).at(0).first; 
	    for (size_t i=0; i<2; i++)
	      theclosesthit[i] = imghits.at(p).at(hitidx)[i];
	  }
	  else {
	    theclosesthit[0] = theclosesthit[1] = -1;
	  }
	  closest_hits_list[p] = theclosesthit;
	  
	  std::vector<int> theplanehit(2);
	  theplanehit[1] = pixel_row;
	  theplanehit[0] = pixel_cols[p];
	  plane_pos_list[p] = theplanehit;
	}
	
	Step3D* next_step = new Step3D( current_pos, current_dir, closest_hits_list, plane_pos_list, *current_step );
	current_step = next_step;
	std::cout << "updated Step3D list." << std::endl;
      }
      
      istep++;
    }
    std::cout << "[End of stopmu tracker. numebr of steps=" << istep << ".]" << std::endl;
    std::cout << "===============================================================" << std::endl;    
  }

  // ------------------------------------------------------------------------------------------------
  // *** Primary Tracker Loop Functions ****
  // ------------------------------------------------------------------------------------------------

  void StopMuTracker::makeProposedPos( const std::vector<float>& currentpos, const std::vector<float>& currentdir, std::vector<float>& proposedpos, const float stepsize ) {
    if ( currentpos.size()!=3 && currentdir.size()!=3 ) {
      throw std::runtime_error("StopMuTracker::makeProposedPos length of input pos and dir vectors not 3");
    }
    float dirnorm = 0.;
    proposedpos.resize(3,0.0);
    for (int i=0; i<3; i++) {
      dirnorm += currentdir[i]*currentdir[i];
      proposedpos[i] = currentpos[i] + currentdir[i]*stepsize;
    }
    dirnorm = sqrt(dirnorm);
    if ( std::fabs(dirnorm-1.0)>1.0e-4 ) {
      throw std::runtime_error("StopMuTracker::makeProposedPos direction vector was not normalized to 1");
    }
  }

  void StopMuTracker::imagePositions( const std::vector<float>& currentpos, int& tick, std::vector<int>& wid ) {
    wid.resize(3,-1); // wire position
    double dpos[3];
    for (int p=0; p<3; p++) dpos[p] = currentpos[p];
    for (int p=0; p<3; p++)
      wid[p] = (int)::larutil::Geometry::GetME()->WireCoordinate( dpos, p );
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    tick = (int)(currentpos[0]/cm_per_tick) + 3200;
  }

  void StopMuTracker::getClosestHitsInPlane( const int clusterid, const std::vector<int>& test_pos, 
					     const ::dbscan::dbPoints& src_data, const ::dbscan::dbscanOutput& cluster_info, const larcv::ImageMeta& meta, 
					     std::vector< std::pair<int,double> >& hitlist ) {
    // a wrapper around dbscanOutput::closestHitsInCluster
    // have to convert the 2D pixel col/row information into double test_point
    std::vector<double> dtest_pos(test_pos.size(),0);
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float cm_per_wire = 0.3;
    
    dtest_pos[0] = (double)meta.col( test_pos[0] );
    dtest_pos[1] = (double)meta.row( test_pos[1] );
    cluster_info.closestHitsInCluster( clusterid, dtest_pos, src_data, meta, cm_per_tick, cm_per_wire, hitlist, 5 );
  }

  // ------------------------------------------------------------------------------------------------
  // *** Re-orientation Routines ****
  // ------------------------------------------------------------------------------------------------

  bool StopMuTracker::findNewDirection() {
    // for now, act as a dummy function and kill the stepper by returning false
    return false;
  }

  
  // ------------------------------------------------------------------------------------------------
  // *** UTILITY FUNCTIONS ****
  // ------------------------------------------------------------------------------------------------
  
  float StopMuTracker::_norm( std::vector<float>& vec ) {
    float norm = 0;
    for (size_t i=0; i<vec.size(); i++) norm += vec[i]*vec[i];
    norm = sqrt(norm);
    for (size_t i=0; i<vec.size(); i++) vec[i] /= norm;
    return norm;
  }

  void StopMuTracker::_wire2pixel( const int tick, const std::vector<int>& wid, const larcv::ImageMeta& meta, std::vector<int>& pixel_col, int& pixel_row ) {
    pixel_col.resize(3,0.0);
    if ( pixel_col.size()!=wid.size() ) {
      throw std::runtime_error("StopMuTracker::_wire2pixel: provided wires for more or less than 3 planes");
    }
    for (int p=0; p<3; p++)
      pixel_col[p] = meta.col( wid[p] );
    pixel_row = meta.row( tick );
  }
  
}
