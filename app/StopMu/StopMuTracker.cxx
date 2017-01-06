#include "StopMuTracker.h"

#include <stdio.h>
#include <time.h>
#include <sstream>

#ifndef __CINT__
//#include <opencv2/opencv.hpp>
//#include <opencv2/core/core.hpp>
//#include "CVUtil/CVUtil.h"
#endif

// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "StopMuSkeleton.h"

namespace larlitecv {

  StopMuTracker::StopMuTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, int verbosity ) {
    
    setVerbosity(verbosity);
				
    
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
      dbscan::dbscanOutput cluster = dbalgo.scan( data, 5, 50.0 );
      m_imghits.emplace_back( data );
      m_clusters.emplace_back( cluster );
      //std::cout << "number of clusters on plane " << p << ": " << m_clusters.at(p).clusters.size() << std::endl;
    }//loop over clusters
    
    time_t cluster_finished = time(NULL);
    //double dt_clusters = difftime(cluster_start,cluster_finished);
    double dt_clusters = cluster_finished-cluster_start;
    //std::cout << "clustered in " << dt_clusters << " seconds." << std::endl;
    
    // initialize hit list
    for (int p=0; p<3; p++) {
      current_hit[p] = -1;
    }

    // for (int p=0; p<3; p++) {
    //   cv::Mat imgmat = larcv::as_mat_greyscale2bgr( skel_v.at(p), 0, 1.0);
    //   std::stringstream ss;
    //   ss << "skel_p" << p << ".jpg";
    //   cv::imwrite( ss.str().c_str(), imgmat );
    // }
  }

  void StopMuTracker::trackStopMu( const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
				   const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d, Step3D& trackstart ) {

    // parameters: move these to configuration file later
    float fCloseHitThreshold_cm = 0.5; // 2 wire disagreement
    float fStepSize_cm = 0.1;

    // find matching cluster
    std::vector<int> clusterid;

    const larcv::ImageMeta& meta = skel_v.at(0).meta();
    larcv::Image2D img_cluster( meta );
    img_cluster.paint(0.0);

    std::vector<Hit2DList> hitlists(3);

    for (int p=0; p<3; p++) {
      std::vector<double> testpoint(2);
      testpoint[0] = start2d.at(p)[0]; // X
      testpoint[1] = start2d.at(p)[1]; // Y

      std::cout << "plane " << p << " testpoint: (col,row)=(" << testpoint[0] << "," << testpoint[1] << ") " 
		<< " (wire,tick)=(" << meta.pos_x( testpoint[0] ) << "," << meta.pos_y( testpoint[1] ) << ")" << std::endl;
      int match = m_clusters.at(p).findMatchingCluster( testpoint, m_imghits.at(p), 10.0 );
      clusterid.push_back(match);
      std::cout << "plane=" << p 
		<< " matching cluster index=" << match 
		<< " size=" << m_clusters.at(p).clusters.at(match).size() << std::endl;

      // we make a sorted list of pixels by distance
      // we also need an initial direction

      for (int ihit=0; ihit<m_clusters.at(p).clusters.at(match).size(); ihit++) {
	int hitidx = m_clusters.at(p).clusters.at(match).at(ihit);
	int x = m_imghits.at(p).at(hitidx)[0];
	int y = m_imghits.at(p).at(hitidx)[1];

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
	if ( norm==0 ) {
	  hit.distance = norm;
	  cosine = 1.0;
	}
	if ( cosine>0 ) {
	  hitlists[p].emplace(std::move(hit));
	  hitlists[p].sort();
	}
	
	img_cluster.set_pixel((int)y,(int)x,250.0);
      }

      
      std::cout << "plane " << p << " number of hits=" << hitlists[p].size() << std::endl;
      for (int ihit=0; ihit<(int)hitlists[p].size(); ihit++) {
	std::cout << " [#" << ihit << "] (" << hitlists[p].at(ihit)[0] << "," << hitlists[p].at(ihit)[1] << ")"
		  << " " << hitlists[p].at(ihit).distance << std::endl; 
      }
      
    }//end of loop over planes for sorting hits
    
    // cv::Mat imgmat = larcv::as_mat( img_cluster );
    // std::stringstream ss;
    // ss << "baka.jpg";
    // cv::imwrite( ss.str().c_str(), imgmat );

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
    // get the first step's plane position
    // (to do, not important at the moment)
    
    Step3D* current_step = &trackstart;

    while ( !isfinished && istep<350 ) {
      
      std::cout << "===============================================================" << std::endl;
      std::cout << " Step " << istep << std::endl;

      // update 3D step
      makeProposedPos( current_pos, current_dir, proposed_pos, fStepSize_cm );

      std::cout << " from current pos=(" << current_pos[0] << "," << current_pos[1] << "," << current_pos[2] << ")"
		<< " and current dir=(" << current_dir[0] << "," << current_dir[1] << "," << current_dir[2] << ")" << std::endl;
      std::cout << " stepsize=" << fStepSize_cm << " cm" << std::endl;
      std::cout << " made proposal: (" << proposed_pos[0] << "," << proposed_pos[1] << "," << proposed_pos[2] << ")" << std::endl;

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
	//getClosestHitsInPlane( clusterid[p], test_pos, imghits.at(p), clusters.at(p), meta, closest_pixels.at(p) );
	getClosestHitsInList( test_pos, hitlists.at(p), meta, closest_pixels.at(p) );
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
	  std::cout << " [ok: p=" << p << ", (" << hitlists.at(p)[closest_pixels.at(p).at(0).first][0] << "," << hitlists.at(p)[closest_pixels.at(p).at(0).first][1] << ")" 
		    << "," << closest_pixels.at(p).at(0).second << " cm] ";
	}
	else
	  std::cout << " [not ok: p=" << p << ", " << closest_pixels.at(p).at(0).second << " cm] ";
      }
      std::cout << std::endl;

      bool continue_search = true;
      if ( ngoodplanes>=2 || istep<10 ) {
	//  everything is fine continue
	current_pos = proposed_pos;
	std::cout << "number of good planes=" << ngoodplanes << ". continue tracker" << std::endl;

	// update tracker
	// prepare wire info
	std::vector<Point2D_t> closest_hits_list(3);
	std::vector<Point2D_t> plane_pos_list(3);
	for (int p=0; p<3; p++) {
	  std::vector<int> theclosesthit(2);
	  if ( closest_pixels.at(p).size()>0 ) {
	    int hitidx = closest_pixels.at(p).at(0).first; 
	    for (size_t i=0; i<2; i++)
	      theclosesthit[i] = m_imghits.at(p).at(hitidx)[i];
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
      else {
	//  current path is no good. need to re-orient or stop
	std::cout << "need to re-direct track" << std::endl;
	Step3D* proposed_step = new Step3D;
	continue_search = findNewDirection( *current_step, meta, closest_pixels, hitlists, *proposed_step );
	std::cout << "re-directed step returned." << std::endl;
	if ( continue_search )
	  current_step = proposed_step;
	current_pos = current_step->pos;
	current_dir = current_step->dir;
	std::cin.get();
      }

      if ( !continue_search ) {
	std::cout << "Stopping condition met." << std::endl;
	break;
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

  void StopMuTracker::getClosestHitsInList( const std::vector<int>& test_pos, const Hit2DList& src_list, const larcv::ImageMeta& meta, 
					    std::vector< std::pair<int,double> >& hitlist ) {
    // a wrapper around dbscanOutput::closestHitsInCluster
    // have to convert the 2D pixel col/row information into double test_point
    std::vector<double> dtest_pos(test_pos.size(),0);
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float cm_per_wire = 0.3;
    
    dtest_pos[0] = (double)meta.col( test_pos[0] );
    dtest_pos[1] = (double)meta.row( test_pos[1] );
    src_list.closestHits( dtest_pos, meta, cm_per_tick, cm_per_wire, hitlist, 5 );
  }

  // ------------------------------------------------------------------------------------------------
  // *** Re-orientation Routines ****
  // ------------------------------------------------------------------------------------------------

  bool StopMuTracker::findNewDirection( Step3D& start_step, const larcv::ImageMeta& meta,
					const std::vector< std::vector< std::pair<int,double> > >& closest_pixels, const std::vector<Hit2DList>& sorted_hits, 
					Step3D& proposed_step ) {
    // we are going to use the sorted_hits list to find a new direction forward
    int nhits_back = 3;
    int nhits_forward = 10;
    int nhits_vary = 3;
    
    // we keep backing up until we can get a new direction... (how not to get stuck?)
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << " StopMuTracker::findNewDirection" << std::endl;
    int ntries = 0;
    Step3D& current_step = start_step;
    while ( !current_step.isStart() && ntries<1 ) {
      std::cout << " try #" << ntries << std::endl;
      // we get the 2d hit positions downstream of the clusters
      std::vector<int> next_planewire(3,-1); // wire of hit
      std::vector<int> hit_index(3,-1); // position in Hit2D list
      std::vector<float> intersectZY;
      double tri_area = 0.0;
      int crosses = 0;
      int ngood_planes = 0;
      std::cout << "get plane wires" << std::endl;
      float ave_row = 0;
      float time_weights = 0;
      for (int p=0; p<3; p++) {
	std::cout << "closest hits plane=" << p << ": " << closest_pixels.at(p).size() << std::endl;
	if ( closest_pixels.at(p).size()>0 && closest_pixels.at(p).at(0).first>=0 ) {
	  hit_index[p] = closest_pixels.at(p).at(0).first+nhits_forward;
	  if ( hit_index[p]>=sorted_hits.at(p).size() )
	    hit_index[p] = sorted_hits.at(p).size()-1;
	  std::cout << "hit_index[" << p << "]=" << hit_index[p] << " (" << sorted_hits.at(p).at( hit_index[p] )[0] << "," << sorted_hits.at(p).at( hit_index[p] )[1] << ")" << std::endl;
	  next_planewire[p] = meta.pos_x( sorted_hits.at(p).at( hit_index[p] )[0] );
	  ngood_planes++;
	  
	  float w = 1.0/(closest_pixels.at(p).at(0).second+0.1);
	  ave_row += sorted_hits.at(p).at( hit_index[p] )[1]*w;
	  time_weights += w;
	}
      }
      ave_row /= time_weights;
      std::cout << "ave. row=" << ave_row << std::endl;

      if ( ngood_planes==2 ) {
	std::cout << "number of good planes=" << ngood_planes << std::endl;
	// only have two good planes. find the intersection point, and determine the 3rd plane's intersection point
	int goodplane_wire[2] = {0};
	int goodplane[2] = {0};
	int iplane = 0;
	int badplane = -1;
	for (int p=0;p<3; p++) {
	  if ( hit_index[p]>=0 ) {
	    goodplane[iplane] = p;
	    goodplane_wire[iplane] = next_planewire[p];// convert col to wireID
	    iplane++;
	  }
	  else {
	    badplane = p;
	  }
	}
	larcv::UBWireTool::wireIntersection( goodplane[0], goodplane_wire[0], goodplane[1], goodplane_wire[1], intersectZY, crosses );
	double dpos[3];
	dpos[0] = 100.0;
	dpos[1] = intersectZY[2];
	dpos[2] = intersectZY[1];
	next_planewire[badplane] = (int)larutil::Geometry::GetME()->WireCoordinate( dpos, badplane );
	tri_area = 0.0; // by definition, this is a "good" intersection point (i.e. 3D consistent)
      }
      else if (ngood_planes==3) {
	// we test the 3D consistency of the pixels
	std::cout << "three good hits: get intersection." << std::endl;
	larcv::UBWireTool::wireIntersection( next_planewire, intersectZY, tri_area, crosses );
      }
      else {
	throw std::runtime_error("only 1 good plane. screwed.");
	break;
      }

      std::cout << " new point wire positions: U=" << next_planewire[0] << " V=" << next_planewire[1] << " Y=" << next_planewire[2] << std::endl;
      std::cout << " intersection triangle area: " << tri_area << std::endl;
      if ( tri_area<0.0 ) {
	std::cout << "consistent enough 3D position." << std::endl;
      }
      else {
	// scan pixels and find a good 3D point
	std::cout << "looking for a more consistent 3D space-point that tri-area=" << tri_area << std::endl;
	int dhits[3] = { -nhits_vary, -nhits_vary, -nhits_vary }; // O(nhits_vary^3) -- can't be very big
	std::vector<int> best_wires(3);
	double best_area = 1.0e6;
	double best_ave_row = 0.;
	int crosses = -1;
	std::vector<float> best_intersection(2,0.0);
	int num_valid_wire_combos = 0;
	float min_invalid_diff = -1.0;
	for (int ivary=0; ivary<nhits_vary*nhits_vary*nhits_vary; ivary++) {

	  std::vector<int> test_wires(3);
	  std::vector<float> test_rows(3);
	  float test_ave_tick = 0.;
	  float test_tick_weight = 0.;
	  for (int p=0; p<3; p++) {
	    int hitidx = hit_index[p] + dhits[p];
	    if ( hitidx<0 )                      hitidx = 0;
	    if ( hitidx>=sorted_hits.at(p).size() ) hitidx = sorted_hits.at(p).size()-1;
	    test_wires[p] = sorted_hits.at(p).at(hitidx)[0];
	    test_rows[p]  = sorted_hits.at(p).at(hitidx)[1];
	    if ( test_wires[p]<0 ) test_wires[p] = 0;
	    if ( test_wires[p]>=(int)larutil::Geometry::GetME()->Nwires(p) ) test_wires[p] = (int)larutil::Geometry::GetME()->Nwires(p)-1;
	  }
	  
	  float max_time_diff = 0.;
	  for (int p1=0; p1<3; p1++) {
	    for (int p2=p1+1; p2<3; p2++) {
	      float tdiff = fabs( test_rows[p1]-test_rows[p2] );
	      if ( tdiff>max_time_diff )
		max_time_diff = tdiff;
	    }
	  }
	  if ( max_time_diff>2 ) {
	    if ( min_invalid_diff<0 || min_invalid_diff>max_time_diff ) 
	      min_invalid_diff = max_time_diff;
	    continue; // not really valid
	  }
	  
	  num_valid_wire_combos++;

	  std::vector<float> test_intersection;
	  double test_area;
	  int test_crosses;
	  larcv::UBWireTool::wireIntersection( test_wires, test_intersection, test_area, test_crosses );
	  if ( test_crosses==1 && test_area<best_area ) {
	    best_area = test_area;
	    best_intersection[0] = test_intersection[0];
	    best_intersection[1] = test_intersection[1];
	    best_wires = test_wires;
	    best_ave_row = 0.0;
	    for (int p=0; p<3; p++)
	      best_ave_row += test_rows[p];
	    best_ave_row /= 3.0;
	  }
	  int move_index = ivary%3;
	  dhits[move_index]++;
	}//end of loop over variations

	std::cout << "best tri-area after scanning: " << best_area << std::endl;
	if ( best_area<tri_area ) {
	  tri_area = best_area;
	  next_planewire = best_wires;
	  intersectZY = best_intersection;
	  ave_row = best_ave_row;
	}
	std::cout << "new best wires: U=" << next_planewire[0] << " V=" << next_planewire[1] << " Y=" << next_planewire[2] << std::endl;
	std::cout << "new best tick=" << ave_row << std::endl;
	std::cout << "new intersection (y,z)=(" << intersectZY[1] << "," << intersectZY[0] << ")" << std::endl;
	std::cout << "number of valid wire combinations tested: " << num_valid_wire_combos << std::endl;
	std::cout << "min. invalid time-diff: " << min_invalid_diff << std::endl;
      }

      // turn average row into x position
      float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
      float new_x = (meta.pos_y(ave_row) - 3200)*cm_per_tick;
      std::vector<float> new_pos(3);
      new_pos[0] = new_x;
      new_pos[1] = intersectZY[1];
      new_pos[2] = intersectZY[0];

      // use the new point to get a new direction to step
      std::vector<float> newdir(3);
      for (int p=0; p<3; p++) {
	newdir[p] = new_pos[p] - current_step.pos[p];
      }
      _norm( newdir );
      
      std::cout << "make plane positions for new step" << std::endl;
      std::vector<Point2D_t> planepositions(3);
      for (int p=0; p<3; p++) {
	std::vector<int> step_hit_pos(2);
	step_hit_pos[0] = meta.col( next_planewire[p] );
	step_hit_pos[1] = (int)ave_row;
	planepositions[p] = step_hit_pos;
      }

      std::cout << "fill out proposed step, connect it to current step" << std::endl;
      proposed_step.pos.resize(3,0);
      for (int p=0; p<3; p++) 
	proposed_step.pos[p] = new_pos[p];//current_step.pos[p] + newdir[p];
      proposed_step.dir = newdir;
      proposed_step.closesthits = planepositions;
      proposed_step.planepositions = planepositions;
      proposed_step.prev = &current_step;
      proposed_step.next = NULL;
      current_step.next = &proposed_step;

      // this gives us a new direction to step
      break;
    }//loop over tries
    
    
    return true;
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

  std::vector<larcv::Image2D> StopMuTracker::fillSortedHit2Dlist( const larcv::ImageMeta& meta, 
								  const std::vector< std::vector<int> >& start2d,  const std::vector< std::vector<float> >& start_dir2d,
								  std::vector<Hit2DList>& hitlists, std::vector<int>& clusterid ) {
    // given a starting pixel, uses stored clusters (m_clusters) and list of dbPoints (m_imghits) 
    // to provide a list of Hit2D objects sorted by distance from start point
    // also returns an image containing which contains the hits that form the cluster in question

    std::cout << "StopMuTracker::fillSortedHit2Dlist" << std::endl;

    std::vector<larcv::Image2D> img_clusters;

    clusterid.resize(3,-1);
    hitlists.clear();
    hitlists.resize(3);
    for (int p=0; p<3; p++) {

      larcv::Image2D img_cluster( meta );
      img_cluster.paint(0.0);
    
      std::vector<double> testpoint(2);
      testpoint[0] = start2d.at(p)[0]; // X
      testpoint[1] = start2d.at(p)[1]; // Y
      std::cout << "plane " << p << " testpoint: (col,row)=(" << testpoint[0] << "," << testpoint[1] << ")" 
		<< " (wire,tick)=(" << meta.pos_x( testpoint[0] ) << "," << meta.pos_y( testpoint[1] ) << ")" << std::endl;
      int match = m_clusters.at(p).findMatchingCluster( testpoint, m_imghits.at(p), 5.0 );
      clusterid.push_back(match);
      std::cout << "plane=" << p << " start_dir2d=(" << start_dir2d.at(p)[0] << "," << start_dir2d.at(p)[1] << ") "
		<< " matching cluster index=" << match 
		<< " size=" << m_clusters.at(p).clusters.at(match).size()
		<< std::endl;

      // we make a sorted list of pixels by distance
      // we also need an initial direction

      for (int ihit=0; ihit<m_clusters.at(p).clusters.at(match).size(); ihit++) {
	int hitidx = m_clusters.at(p).clusters.at(match).at(ihit);
	int x = m_imghits.at(p).at(hitidx)[0];
	int y = m_imghits.at(p).at(hitidx)[1];

	// dir from start to point. dir2d derives from cm-scales
	std::vector<float> dir(2);
	dir[0] =  (x-start2d.at(p)[0])*meta.pixel_width()*0.3; // cm
	dir[1] = -(y-start2d.at(p)[1])*meta.pixel_height()*0.5*::larutil::LArProperties::GetME()->DriftVelocity(); // cm
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
	if ( norm==0 ) {
	  hit.distance = 0;
	  cosine = 1;
	}

	if ( cosine>0 ) {
	  hitlists[p].emplace(std::move(hit));
	  hitlists[p].sort();
	}
	
	//std::cout << " hitidx=" << hitidx << ": (" << hit[0] << "," << hit[1] << ") cosine=" << cosine << std::endl;
	img_cluster.set_pixel((int)y,(int)x,250.0);
      }

      
      std::cout << "plane " << p << " number of hits=" << hitlists[p].size() << std::endl;
//       for (int ihit=0; ihit<(int)hitlists[p].size(); ihit++) {
// 	std::cout << " [#" << ihit << "] (" << hitlists[p].at(ihit)[0] << "," << hitlists[p].at(ihit)[1] << ")"
// 		  << " " << hitlists[p].at(ihit).distance << std::endl; 
//       }

      img_clusters.emplace_back( std::move(img_cluster) );
      
    }//end of loop over planes for sorting hits
    
    return img_clusters;
  }// StopMuTracker::fillSortedHit2Dlist

  // ------------------------------------------------------------------------------------------------
  // *** Hit2DList Functions ****
  // ------------------------------------------------------------------------------------------------

  void Hit2DList::closestHits( std::vector<double>& test_pos, const larcv::ImageMeta& meta, const float cm_per_tick, const float cm_per_wire,
			       std::vector< std::pair<int,double> >& hitlist, const int max_nhits, const int ignore_marked ) const {
    struct mycompare_t {
      bool operator() ( std::pair<int,double>& lhs, std::pair<int,double>& rhs ) { 
	if ( lhs.second<rhs.second ) 
	  return true;
	return false;
      }
    } mycompare;

    hitlist.clear();
    for (int ihit=0; ihit<(int)size(); ihit++) {
      double hit_dist  = 0.0;
      const Hit2D& hitpos = at(ihit);
      if ( ignore_marked==1 && hitpos.marked==1 ) continue;
      double dt = ((double)hitpos[1]-test_pos[1])*meta.pixel_height()*cm_per_tick;
      double dw = ((double)hitpos[0]-test_pos[0])*meta.pixel_width()*cm_per_wire;
      hit_dist = sqrt( dt*dt + dw*dw ); // in cm
      std::pair<int,double> test_hit( ihit, hit_dist );
      if ( max_nhits>0 ) {
	// we have to care about the number of hits in the list
	if ( hitlist.size()<max_nhits ) {
	  // but not now
	  hitlist.emplace_back( test_hit );
	}
	else {
	  bool isbetter = false;
	  for (size_t ilisthit=0; ilisthit<hitlist.size(); ilisthit++) {
	    if ( hitlist.at(ilisthit).second > hit_dist ) {
	      isbetter = true;
	      break;
	    }
	  }
	  if ( isbetter )
	    hitlist.emplace_back( test_hit );
	}
      }
      else {
	// dont care about the number of hits. add it to hit list.
	hitlist.emplace_back( test_hit );
      }
      
      // sort the current hitlist vector
      std::sort( hitlist.begin(), hitlist.end(), mycompare );
      // truncate the end
      if ( max_nhits>0 && max_nhits<(int)hitlist.size() )
	hitlist.resize(max_nhits);
            
    }//end of loop over hits

    if ( _verbose_ ) {
      std::cout << "sorted distances from test point (total hits=" << size() << "):" << std::endl;
      for (size_t i=0; i<hitlist.size(); i++) {
	  std::cout << " #" << i << ": " << hitlist.at(i).first << " " << hitlist.at(i).second << std::endl;
      }
    }

  }

  // ------------------------------------------------------------------------------------------------
  // *** String Minimizer ****
  // ------------------------------------------------------------------------------------------------

  FitDataHack* FitDataHack::_global_instance = NULL;

  double stop_mu_position_score(const double *xx )
  {
    // fit score comes from these components:
    // 1) how close it is to a pixel with charge above threshold
    // 2) how straight the step is with respect to the last step
    // 3) does the charge expected along the step consistent with that seen?
    
    // needs to find the closest hit
    Double_t cosz = xx[0];
    Double_t phi  = xx[1];
    bool _verbose_ = false;
    
    const larcv::ImageMeta& meta = *(FitDataHack::getMe()->data_meta);
    const std::vector<larcv::Image2D>& img_v = *(FitDataHack::getMe()->data_images);
    const std::vector<Hit2DList>& hitlists = *(FitDataHack::getMe()->data_hitlist);
    const std::vector<float>& anchor   = *(FitDataHack::getMe()->data_anchor);
    const std::vector<float>& prev_dir = *(FitDataHack::getMe()->data_prev_dir);
    const float dist_weight   = FitDataHack::getMe()->dist_weight;
    const float bend_weight   = FitDataHack::getMe()->bend_weight;
    const float charge_weight = FitDataHack::getMe()->charge_weight;
    const float step_size     = FitDataHack::getMe()->step_size_cm;

    // first make point and direction
    std::vector<float> dir(3,0.0);
    std::vector<float> pos(3,0.0);
    float rphi = sqrt(1-cosz*cosz);
    dir[2] = cosz;
    dir[1] = rphi*sin(phi);
    dir[0] = rphi*cos(phi);
    for (int i=0; i<3; i++) pos[i] = anchor[i] + step_size*dir[i];

    // get plane positions
    int tick;
    std::vector<int> wids;
    larlitecv::StopMuTracker::imagePositions( pos, tick, wids );
    if ( tick<meta.min_y() )  tick = meta.min_y();
    if ( tick>=meta.max_y() ) tick = meta.max_y()-1;
    
    int tick_anchor;
    std::vector<int> wids_anchor;
    larlitecv::StopMuTracker::imagePositions( anchor, tick_anchor, wids_anchor );
    if ( tick_anchor<meta.min_y() )  tick_anchor = meta.min_y();
    if ( tick_anchor>=meta.max_y() ) tick_anchor = meta.max_y()-1;

    // find the closest point on each plane
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float cm_per_wire = 0.3;
    std::vector<float> closest_dists(3,1.0e3);
    float total_dist = 0.;
    if ( _verbose_ )
      std::cout << "test (" <<  cosz << "," << phi << ") distances from (" << pos[0] << "," << pos[1] << "," << pos[2] << "): ";
    for (int p=0; p<3; p++) {
      const Hit2DList& planehits = hitlists.at(p);
      std::vector< std::pair<int,double> > closest;
      std::vector<double> testpoint(2,0.0);
      testpoint[0] = meta.col( wids[p] );
      testpoint[1] = meta.row( tick );
      if ( testpoint[0]< 0 ) testpoint[0];
      if ( testpoint[0]>=meta.cols() ) testpoint[0] = meta.cols()-1;
      if ( testpoint[1]< 0 ) testpoint[1] = 0;
      if ( testpoint[1]>=meta.rows() ) testpoint[1] = meta.rows()-1;
      planehits.closestHits( testpoint, meta, cm_per_tick, cm_per_wire, closest, 3, 1 );
      if ( closest.size()>0 && closest.front().second<closest_dists[p]) {
	closest_dists[p] = closest.front().second;
	if ( _verbose_ )
	  std::cout << " p" << p << "=" << closest.front().second;
      }
      total_dist+= closest_dists[p]*closest_dists[p];
    }

    // cosine between previous and current direction
    float cosine = 0.;
    for (int i=0; i<3; i++) cosine += prev_dir[i]*dir[i];
    float cos_score = 1.0-cosine;
    if ( _verbose_ )
      std::cout << " bend cos=" << cosine;
    
    // charge along the path. 
    // we use Y-only because the calorimetry there is more reliable
    // what is the charge per y-pixel
    float mip_adc_per_voxel = 30.0;
    
    // first count the charge along the step
    int n_ypixels = 0;
    float tot_adc_on_yplane = 0.;
    std::vector<float> ydir(2);
    ydir[0] = float(meta.col(wids[2]))-float(meta.col(wids_anchor[2]));
    ydir[1] = float(meta.row(tick)) - float(meta.row(tick_anchor));
    float ynorm = sqrt( ydir[0]*ydir[0] + ydir[1]*ydir[1] );
    if ( fabs(ydir[0])<1.0e-3 && fabs(ydir[1])<1.0e-3 ) {
      int col = meta.col(wids[2]);
      int row = meta.row(tick);
      if ( col<0 ) col = 0;
      if ( col>=meta.cols() ) col = meta.cols()-1;
      if ( row<0 ) row = 0;
      if ( row>=meta.rows() ) row = meta.rows()-1;
      tot_adc_on_yplane = img_v[2].pixel( row, col );
      n_ypixels++;
    }
    else {
      // need to step through pixels
      //std::cout << " ynorm=" << ynorm << " ydir=(" << ydir[0] << "," << ydir[1] << ") dcol=" << meta.col(wids[2]) << "-" << meta.col(wids_anchor[2]) << " ";
      ydir[0] /= ynorm;
      ydir[1] /= ynorm;
      float ncols = fabs( float(meta.col(wids[2]))-float(meta.col(wids_anchor[2])) );
      float nrows = fabs( float(meta.row(tick)) - float(meta.row(tick_anchor)) );
      if ( _verbose_ )
	std::cout << " nrows=" << nrows << " ncols=" << ncols << std::endl;
      if ( fabs(ydir[0])>fabs(ydir[1]) ) {
	// cols faster moving than rows
	for (int icol=0; icol<=(int)ncols; icol++) {
	  int row = float(meta.row(tick_anchor)) + (ydir[1]/ydir[0])*icol;
	  int col = float(meta.col(wids_anchor[2])) + icol;
	  if ( col<0 ) col = 0;
	  if ( col>=meta.cols() ) col = meta.cols()-1;
	  if ( row<0 ) row = 0;
	  if ( row>=meta.rows() ) row = meta.rows()-1;
	  tot_adc_on_yplane += img_v[2].pixel( row, col );
	  n_ypixels++;
	}
      }
      else {
	// rows faster than cols
	for (int irow=0; irow<=(int)nrows; irow++) {
	  int col = float(meta.col(wids_anchor[2])) + (ydir[0]/ydir[1])*irow;
	  int row = float(meta.row(tick_anchor)) + irow;
	  if ( col<0 ) col = 0;
	  if ( col>=meta.cols() ) col = meta.cols()-1;
	  if ( row<0 ) row = 0;
	  if ( row>=meta.rows() ) row = meta.rows()-1;
	  tot_adc_on_yplane += img_v[2].pixel( row, col );
	  n_ypixels++;
	}	
      }
    }

    // now predict the charge for the step
    // basically, to do this, we need to create a voxel whose width is the wire pitch, the height is the dirft distance for one tick
    // and length is the wire. The latter is in principle, but in practice the track could end (since its stopping) any time 
    
    // ok, so we need to know which boundary the step hits first: (yz, zy, xz)
    float sx = (1.0*cm_per_tick)/dir[0];  // crossing over to the next tick
    float sy= (117.0-anchor[1])/dir[1]; // the actual cell
    float sz = (1.0*cm_per_wire)/dir[2];  // crossing over to the next wire
    if ( dir[1]<0 )
      sy = (-117.0-anchor[1])/dir[1]; 
    float s[3] = { sx, sy, sz };
    
    int shortest_dir = -1;
    for (int i=0; i<3; i++) {
      if ( dir[i]!=0 && (shortest_dir==-1 || shortest_dir>fabs(s[i]) ) ) {
	shortest_dir = i;
      }
    }
    
    // where does this point end up?
    std::vector<float> voxel_cross(3,0.0);
    float voxel_dist = 0.;
    for (int i=0; i<3; i++) {
      voxel_cross[i] = fabs(s[shortest_dir])*dir[i];
      voxel_dist += voxel_cross[i]*voxel_cross[i];
    }
    voxel_dist = sqrt(voxel_dist);
    
    float total_voxel_dist = voxel_dist*n_ypixels;
    
    float path_adc = mip_adc_per_voxel*(total_voxel_dist/cm_per_wire); // the ADC scale is assuming the wire crosses perpendicularly under the wire
    float adc_diff = path_adc-tot_adc_on_yplane;
    
    // final score
    if ( _verbose_ ) {
      std::cout << " voxel_dist=" << voxel_dist << " n_ypixels=" << n_ypixels;
      std::cout << " path_adc=" << path_adc << " sum_adc=" << tot_adc_on_yplane << " adc_diff=" << adc_diff << std::endl;
    }

    if ( _verbose_ )
      std::cout << std::endl;
    
    double score = 0.5*dist_weight*(total_dist) + 0.5*bend_weight*cos_score*cos_score + 0.5*charge_weight*adc_diff*adc_diff;
    
    return score;
  }

  void StopMuTracker::stopMuString( const std::vector<larcv::Image2D>& img_v, 
				    const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
				    const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d, Step3D& trackstart ) {

    // parameters: move these to configuration file later
    float fCloseHitThreshold_cm = 0.5; // 2 wire disagreement
    float fStepSize_cm = 0.3; // 3 wires or ~10 microseconds or 20 ticks, which are a little more than 3 rows
    //float fStepSize_cm = 1.0; // 3 wires or ~10 microseconds or 20 ticks, which are a little more than 3 rows
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec]*[usec/tick]
    float cm_per_wire = 0.3;

    // find matching cluster
    std::vector<int> clusterid;
    const larcv::ImageMeta& meta = skel_v.at(0).meta();
    std::vector<Hit2DList> hitlists(3);
    std::vector<larcv::Image2D> img_clusters = fillSortedHit2Dlist( meta, start2d, start_dir2d, hitlists, clusterid );

    // for debug output
    // if ( m_verbosity>2 ) {
    //   for (int p=0; p<3; p++) {
    // 	const larcv::Image2D& img_cluster = img_clusters.at(p);
    // 	cv::Mat imgmat = larcv::as_mat( img_cluster );
    // 	std::stringstream ss;
    // 	ss << "baka_p" << p << ".jpg";
    // 	cv::imwrite( ss.str().c_str(), imgmat );
    //   }
    // }

    // we need to check the quality of the cluster
    bool clusterok = true;
    for (size_t p=0; p<3; p++) {
      if ( hitlists.at(p).size()<1 ) {
	clusterok = false;
	break;
      }
    }

    if ( clusterok==false ) {
      return;
    }

    int istep = 0;
    bool isfinished = false;

    // we build 3D steps by fitting the position.  The step size is fixed, but the 3D angle is chosen. 
    // the potential is the closest distance to a hit + a "stress" term that wants to keep the line straight

    // current step, direction
    std::vector< float > current_pos =  start_pos3d;
    std::vector< float > current_dir =  start_dir3d;
    std::vector< float > proposed_pos(3,0.0);
    std::vector< float > proposed_dir(3,0.0);
    _norm( current_dir );

    // make the first Step
    trackstart.pos = current_pos;
    trackstart.dir = current_dir;
    trackstart.closesthits = start2d;
    trackstart.planepositions = start2d;

    // get the first step's plane position
    // (to do, not important at the moment)
    
    Step3D* current_step = &trackstart;

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
    minimizer->SetMaxFunctionCalls(10000);
    minimizer->SetMaxIterations(10000); 
    minimizer->SetTolerance(0.001);
    minimizer->SetPrintLevel(0);
      
    ROOT::Math::Functor func(&stop_mu_position_score,2);
    double step[2] = {0.1,0.1};

    minimizer->SetFunction(func);

    while ( !isfinished && istep<1000 ) {
     
      if ( m_verbosity > 0 ) {
	std::cout << "===============================================================" << std::endl;
	std::cout << " Step " << istep << " [ptr=" << &(*current_step) << "]" << std::endl;
      }

      current_pos = current_step->pos;
      current_dir = current_step->dir;

      // update 3D step
      makeProposedPos( current_pos, current_dir, proposed_pos, fStepSize_cm );

      if ( m_verbosity > 0 ) {
	std::cout << " from current pos=(" << current_pos[0] << "," << current_pos[1] << "," << current_pos[2] << ")"
		  << " and current dir=(" << current_dir[0] << "," << current_dir[1] << "," << current_dir[2] << ")" << std::endl;
	std::cout << " stepsize=" << fStepSize_cm << " cm" << std::endl;
	std::cout << " made proposal: " << proposed_pos[0] << "," << proposed_pos[1] << "," << proposed_pos[2] << ")" << std::endl;
      }

      double min_value = 0.;
      double min_cosine = 2.0;
      if ( istep>=1 ) {
	// we need to minimize it
	
	// pass in the data we need
	FitDataHack::getMe()->data_images = &img_v;
	FitDataHack::getMe()->data_meta = &meta;
	FitDataHack::getMe()->data_hitlist = &hitlists;
	FitDataHack::getMe()->data_anchor = &current_pos;
	FitDataHack::getMe()->data_prev_dir = &( current_step->GetPrev().dir );
	FitDataHack::getMe()->dist_weight = 1.0;
	FitDataHack::getMe()->bend_weight = 0.1;
	FitDataHack::getMe()->charge_weight = 0.00;
	FitDataHack::getMe()->step_size_cm = 2.0;

	// set the initial values
	double variable[2] = { 0.98, 0.0 };
	// get the direction into cosz,phi
	variable[0] = current_dir[2];
	variable[1] = atan2(current_dir[1],current_dir[0]);

	minimizer->SetLimitedVariable(0, "cosz", variable[0], step[0], -1.0, 1.0 );
	minimizer->SetLimitedVariable(1,  "phi", variable[1], step[1], -3.14159, 3.14159 );

	// do the minimization
	if ( m_verbosity > 0 ) {
	  std::cout << "initialize minimizer with cosz=" << variable[0] << " phi=" << variable[1] << std::endl;
	  std::cout << "Run Minimizer." << std::endl;
	}
	minimizer->Minimize();

	const double *xs = minimizer->X();
	if ( m_verbosity > 0 )
	  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << minimizer->MinValue()  << std::endl;
	proposed_dir[2] = xs[0];
	proposed_dir[0] = sqrt(1.0-xs[0]*xs[0])*cos(xs[1]);
	proposed_dir[1] = sqrt(1.0-xs[0]*xs[0])*sin(xs[1]);
	_norm(proposed_dir);
	for (int i=0; i<3; i++) proposed_pos[i] = current_pos[i] + proposed_dir[i]*fStepSize_cm;
	if ( m_verbosity > 0 ) {
	  std::cout << "post-min pos: (" << proposed_pos[0] << "," << proposed_pos[1] << "," << proposed_pos[2] << ")" << std::endl;
	  std::cout << "post-min dir: (" << proposed_dir[0] << "," << proposed_dir[1] << "," << proposed_dir[2] << ")" << std::endl;
	}
	min_value = minimizer->MinValue();
	min_cosine = 0.;
	for (int i=0; i<3; i++) min_cosine += current_dir[i]*proposed_dir[i];
      }
      else {
	// we skip the minimization for the first step in order to get a previous direction for input
	proposed_dir = current_dir;
      }

      // stopping condition test
      if ( m_verbosity > 0 )
	std::cout << "post-min cosine: " << min_cosine << std::endl;
      if ( min_cosine<-0.2 ) {
	if ( m_verbosity > 0 )
	  std::cout << "Stopping condition met." << std::endl;
	break;
      }

      // calculate 2D positions
      int tick;
      std::vector<int> wid;
      imagePositions( proposed_pos, tick, wid );
      if ( tick<meta.min_y() )  tick = meta.min_y();
      if ( tick>=meta.max_y() ) tick = meta.max_y()-1;

      if ( m_verbosity > 0 )
	std::cout << "proposal location on the 2D planes: tick=" << tick << " planes U=" << wid[0] << " V=" << wid[1] << " Y=" << wid[2] << std::endl;

      std::vector<int> pixel_cols;
      int pixel_row;
      _wire2pixel( tick, wid, meta, pixel_cols, pixel_row );
      if ( m_verbosity > 0 ) 
	std::cout << "proposal location on the image: row=" << pixel_row << " plane col U=" << pixel_cols[0] << " V=" << pixel_cols[1] << " Y=" << pixel_cols[2] << std::endl;


      // update tracker
      // prepare wire info
      if ( m_verbosity > 0 )
	std::cout << "update tracker/prepare wire info" << std::endl;
      std::vector<Point2D_t> closest_hits_list(3);
      std::vector<Point2D_t> plane_pos_list(3);
      bool close_enough = true;
      for (int p=0; p<3; p++) {
	std::vector<int> theplanehit(2);
	theplanehit[1] = pixel_row;
	theplanehit[0] = pixel_cols[p];
	plane_pos_list[p] = theplanehit;

	std::vector<double> dhit(2,0);
	dhit[0] = theplanehit[0];
	dhit[1] = theplanehit[1];
	std::vector<int> theclosesthit(2,0);
	std::vector< std::pair<int,double> > close_list;
	if ( m_verbosity>1 )
	  hitlists[p]._verbose_ = true;
	else
	  hitlists[p]._verbose_ = false;
	hitlists[p].closestHits( dhit, meta, cm_per_tick, cm_per_wire, close_list, 3, 1 );
	hitlists[p]._verbose_ = false;
	if ( m_verbosity > 0 )
	  std::cout << "mark hit up to: ";
	if ( close_list.size()>0 ) {
	  int hitidx = close_list.front().first;
	  hitlists[p].markUpTo( hitidx );
	  theclosesthit[0] = hitlists[p].at( hitidx )[0];
	  theclosesthit[1] = hitlists[p].at( hitidx )[1];
	  closest_hits_list[p] = theclosesthit;
	  if ( m_verbosity > 0 ) {
	    std::cout << " p" << p << "=" << hitidx 
		      << " closest hit=(" << theclosesthit[0] << "," << theclosesthit[1] << ")"
		      << std::endl;
	  }
	  if ( close_list.front().second>0.3*20.0 ) {
	    close_enough = false;
	  }
	}
	else {
	  close_enough = false;
	}
      }
      
      Step3D* next_step = new Step3D( proposed_pos, proposed_dir, closest_hits_list, plane_pos_list, *current_step );
      current_step = next_step;
      if ( m_verbosity > 0 )
	std::cout << "updated Step3D list." << std::endl;

      bool hits_left = true;
      for (int p=0; p<3; p++) {
	if ( hitlists[p].back().marked==1 ) {
	  hits_left = false;
	  break;
	}
      }
	
      if ( !hits_left ) {
	if ( m_verbosity > 0 )
	  std::cout << "no more hits in the cluster." << std::endl;
	break;
      }
      if ( !close_enough ) {
	if ( m_verbosity > 0 )
	  std::cout << "we have lost our way." << std::endl;
	break;
      }

	
      istep++;
    }
    
    if ( m_verbosity>0 ) {
      std::cout << "[End of stopmu tracker. numebr of steps=" << istep << ".]" << std::endl;
      std::cout << "===============================================================" << std::endl;    
    }
  }
  
}
