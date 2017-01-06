#include <iostream>
#include <sstream>
#include <exception>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

#ifndef __CINT__
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/EventROI.h"
#include "CVUtil/CVUtil.h"

// larlitecv
#include "ThruMu/EmptyChannelAlgo.h"

// larlite/stopmu
#include "StopMuAlgoTypes.h"
#include "StopMuFilterSpacePoints.h"
#include "StopMuStart.h"
#include "StopMuTracker.h"
#include "StopMuSkeleton.h"


int main( int nargs, char** argv ) {

  std::cout << "Test the stop mu tracker." << std::endl;

  // config file
  std::string cfg_file = "stopmu.cfg";

  // data coordinator: load an example data file
  larlitecv::DataCoordinator dataco;
  //dataco.add_inputfile( "output_larcv_testextbnb_001.root", "larcv" );
  //dataco.add_inputfile( "~/data/larbys/cosmic_tagger/mcc7_bnbcosmic/thrumu_testmcbnbcosmic_signalnumu.root","larcv");
  dataco.add_inputfile( "output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "StopMu" );
  dataco.initialize();

  // --------------------------------------------
  // Instatiate Algos
  
  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet stopmu_cfg = cfg.get<larcv::PSet>("StopMu");
  
  // filter out end points we are unlikely to be interested in: duplicates and those already used as thru-mu end points
  larlitecv::StopMuFilterSpacePointsConfig stopmu_filter_cfg;
  stopmu_filter_cfg.filter_thrumu_by_tag = true;
  stopmu_filter_cfg.duplicate_radius_cm = 20.0;
  stopmu_filter_cfg.pixel_threshold = 10.0;
  stopmu_filter_cfg.track_width_pixels = 5;
  stopmu_filter_cfg.row_tag_neighborhood = 10;
  stopmu_filter_cfg.col_tag_neighborhood = 10;
  larlitecv::StopMuFilterSpacePoints stopmu_filterpts(stopmu_filter_cfg);
  
  // start point direction
  larlitecv::StopMuStart start_finder_algo;
  start_finder_algo.setVerbose(1);
  float fThreshold = 10.0;
  int rneighbor = 10;
  int cneighbor = 10;


  int nentries = dataco.get_nentries("larcv");
  int start_entry = 0;
  //nentries = 1;


  for (int ientry=start_entry; ientry<start_entry+nentries; ientry++) {

    std::cout << "=====================================================" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << " Entry " << ientry << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // go to some entry
    dataco.goto_entry(ientry, "larcv");
    
    // for the test, we target a top-passing, stop muon
    // at time around tick 3880
    larcv::EventImage2D* imgs            = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
    larcv::EventImage2D* marked3d        = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
    larcv::EventROI* rois                = (larcv::EventROI*)dataco.get_larcv_data( larcv::kProductROI, "tpc" );
    
    // make the bad channel image
    larlitecv::EmptyChannelAlgo emptyalgo;
    larcv::EventChStatus* ev_status      = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
    std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status ); //< fill this out!
    std::cout << "number of bad ch imgs: " << badch_v.size() << std::endl;
    
    // get the imgs and the thru-mu tagged images
    const std::vector<larcv::Image2D>& img_v    = imgs->Image2DArray();
    const std::vector<larcv::Image2D>& thrumu_v = marked3d->Image2DArray();
    const larcv::ImageMeta& meta = img_v.at(0).meta();
    
    // make a list of the EventPixel2D containers
    std::vector< larcv::EventPixel2D* > ev_pixs;
    std::vector<std::string> endpt_list = stopmu_cfg.get< std::vector<std::string> >( "EndPointProducers" );
    int tot_endpts = 0;
    std::cout << "End points: " << std::endl;
    for ( auto &producer : endpt_list ) {
      larcv::EventPixel2D* evpix = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, producer );
      std::cout << "  " << producer << "=" << evpix->Pixel2DArray(0).size() << std::endl;
      tot_endpts += evpix->Pixel2DArray(0).size();
      ev_pixs.push_back( evpix );
    }
    std::cout << "total end points=" << tot_endpts << std::endl;
    
    // --------------------------------------------
    // Output Data objects
  
    // output stopmu-tagged pixels
    std::vector<larcv::Image2D> stopmu_v;
    for (size_t p=0; p<img_v.size(); p++) {
      larcv::Image2D stopmu_img( img_v.at(p).meta() );
      stopmu_img.paint(0);
      stopmu_v.emplace_back( std::move(stopmu_img) );
    }
    
    // --------------------------------------------
    // Algo Prep
    std::vector< std::vector< const larcv::Pixel2D* > > stopmu_candidate_endpts = stopmu_filterpts.filterSpacePoints( ev_pixs, thrumu_v, badch_v );

    // --------------------------------------------
    // stop mu tracker
    
    larlitecv::StopMuTracker sttracker( img_v, thrumu_v, stopmu_candidate_endpts );
    sttracker.setVerbosity(0);

    std::cout << " Number of candidate stop-mu start points: " << stopmu_candidate_endpts.size() << std::endl;

    // make a opencv object for drawing's sake
    std::vector< cv::Mat > cvimgs_v;
    for (int p=0; p<3; p++) {
      cv::Mat imgmat = larcv::as_mat_greyscale2bgr( img_v.at(p), 5.0, 50.0 );
      
      // and label thru-mu pixels
      for (int r=0; r<img_v.at(p).meta().rows(); r++) {
	for (int c=0; c<img_v.at(p).meta().cols(); c++) {
	  if ( thrumu_v.at(p).pixel(r,c)>0 ) {
	    imgmat.at< cv::Vec3b >( cv::Point(c,r) )[0] = (unsigned char)200;
	    imgmat.at< cv::Vec3b >( cv::Point(c,r) )[1] = (unsigned char)0;
	    imgmat.at< cv::Vec3b >( cv::Point(c,r) )[2] = (unsigned char)0;
	  }
	}
      }

      // draw interaction BBox
      larcv::draw_bb( imgmat, img_v.at(p).meta(), rois->ROIArray().at(0).BB().at(p), 0, 200, 0, 1 );
      
      cvimgs_v.emplace_back( std::move(imgmat) );
    }

    bool run_tracker = true;
    if ( stopmu_candidate_endpts.size()>stopmu_cfg.get<int>("StopMuEndPointLimits") ) {
      // this is a bad situation for this code.  we skip it for now, because we need to bring bigger guns
      std::cout << "too many candidates. this is probably an indication of an event with too much cosmic ray activity for current code." << std::endl;
      run_tracker = false;
    }
    
    // OK now we have a list of end points to loop through
    int ipix = 0;
    for ( auto &pix_v : stopmu_candidate_endpts ) {
      std::cout << "==============================================================================" << std::endl;
      std::cout << "Candidate Pixel #" << ipix << std::endl;
      
      // make starting point object
      std::vector<larcv::Pixel2D> start;
      for (size_t p=0; p<pix_v.size(); p++) {
	larcv::Pixel2D pix( *pix_v.at(p) );
	start.emplace_back( pix );
      }

      std::cout << "Pixel position in plane views: tick=" << img_v.at(0).meta().pos_y( start.at(0).Y() ) 
		<< " wids=(" << img_v.at(0).meta().pos_x( start.at(0).X() ) << ","
		<< img_v.at(1).meta().pos_x(start.at(1).X() ) << ","
		<< img_v.at(2).meta().pos_x(start.at(2).X() ) << ")" << std::endl;

      
      // starting dir and position
      std::cout << "--- getting starting point's direction and 3D position ---" << std::endl;
      std::vector<float> start_spacepoint;
      std::vector< std::vector<float> > start_dir2d;
      std::vector<float> start_dir3d;
      try {
	start_finder_algo.getStartDirectionV( img_v, badch_v, start, rneighbor, cneighbor, fThreshold, start_spacepoint, start_dir2d, start_dir3d );
      }
      catch (const std::exception& e) {
	std::cout << "candidate ipixel=" << ipix << " doesn't return a good start direction. error: " << e.what() << std::endl;
	for (int p=0; p<3; p++) {
	  cv::Mat imgmat = cvimgs_v.at(p);
	  // draw the start position
	  cv::circle(imgmat,cv::Point(start.at(p).X(),start.at(p).Y()), 5, cv::Scalar(0,255,0),-1);//
	}
	ipix++;
	continue;
      }

      std::cout << "3D position=(" << start_spacepoint[0] << "," << start_spacepoint[1] << "," << start_spacepoint[2] << ") " << std::endl;
      std::cout << "3D direction=(" << start_dir3d[0] << "," << start_dir3d[1] << ","<< start_dir3d[2] << ")" << std::endl;
      std::cout << "Plane directions: "
		<< " p0=(" << start_dir2d.at(0)[0] << "," << start_dir2d.at(0)[1] << ") "
		<< " p1=(" << start_dir2d.at(1)[0] << "," << start_dir2d.at(1)[1] << ") "
		<< " p2=(" << start_dir2d.at(2)[0] << "," << start_dir2d.at(2)[1] << ") " 
		<< std::endl;
      
      // translate start point into simpler object
      std::vector< std::vector<int> > start2d_pos;
      for (int p=0; p<3; p++) {
	std::vector<int> pos(2);
	pos[0] = start.at(p).X();
	pos[1] = start.at(p).Y();
	start2d_pos.emplace_back( pos );
      }

      std::cout << "--- use start point to track stop-mu  ---" << std::endl;
      bool tracked = false;
      larlitecv::Step3D start_track; // a linked list
      enum { ok=0, unknown, clustererror };
      int tracker_err = ok;
      if ( run_tracker ) {

	try {
	  sttracker.stopMuString( img_v, start2d_pos, start_dir2d, start_spacepoint, start_dir3d, start_track );
	  tracked = true;
	  tracker_err = ok;
	}
	catch (const std::exception& e) {
	  std::cout << "went pear shaped. error: " << e.what() << std::endl;
	  tracked = false;
	  tracker_err = unknown;
	  if ( e.what()==std::string("no cluster error") )
	    tracker_err = clustererror;
	}
      }
      
      // add dot information to cv image
      for (int p=0; p<3; p++) {
	cv::Mat imgmat = cvimgs_v.at(p);
	// draw the start position
	cv::Scalar dotcolor(255,0,255);
	if ( tracker_err==unknown )
	  dotcolor = cv::Scalar(0,200,0);
	else if (tracker_err==clustererror)
	  dotcolor = cv::Scalar(200,200,0);
	cv::circle(imgmat,cv::Point(start.at(p).X(),start.at(p).Y()), 5, dotcolor,-1);//
      }
      
      // tag stopmu pixels
      //for (int p=0; p<3; p++) {
      //const larcv::Image2D& img = img_v.at(p);
      //larcv::Image2D& stopmu = stopmu_v.at(p);
      //}
      
      int nsteps = 0;
      int non_3plane_steps = 0;
      larlitecv::Step3D* current_step = &start_track;
      do { 
	//std::cout << " nstep=" << nsteps << " addr=" << current_step << std::endl;
	if ( current_step->planepositions.size()!=3 ) {
	  // why?
	  non_3plane_steps++;
	  nsteps++;
	  if ( !current_step->isEnd() )
	    current_step = &(current_step->GetNext());
	  continue;
	}

	// label track
	for (int p=0; p<3; p++) {
	  cv::Mat& cvimg = cvimgs_v.at(p);
	  cv::circle( cvimg, cv::Point(current_step->planepositions.at(p)[0], current_step->planepositions.at(p)[1]), 1, cv::Scalar(0,0,200), -1 );
	}
	
	if ( !current_step->isEnd() )
	  current_step = &(current_step->GetNext());
	nsteps++;
      } while ( !current_step->isEnd() );
      
      std::string status = "good";
      if ( !tracked ) status = "bad";
      std::cout << "Pixel #" << ipix << ": produced a " << status << " track with " << nsteps << " steps. (" << non_3plane_steps << " non-3 plane steps)" << std::endl;
      
      // make/store 3D trajectory
      ipix++;
      
      // destroy linked list
      std::cout << "clean up track" << std::endl;
      while ( !current_step->isStart() ) {
	//std::cout << " current=" << current_step << " prev=" << &(current_step->GetPrev()) << std::endl;
	current_step = &(current_step->GetPrev());
	current_step->removeNext();
      }

    }//end of candidate stopmu end points

    for (int p=0; p<3; p++) {
      cv::Mat& imgmat = cvimgs_v.at(p);
      std::stringstream ss;
      ss << "stopmu_n" << ientry << "_r" << dataco.run() << "_s" << dataco.subrun() << "_e" << dataco.event() << "_p" << p << ".png";
      cv::imwrite( ss.str(), imgmat );
    }
    
    // store
    
    larcv::EventImage2D* stopmu_eventimgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stopmu" );
    stopmu_eventimgs->Emplace( std::move(stopmu_v) );
    
    dataco.save_entry();

    //if ( ientry>=10 )
    //break;
    
  }// loop over entries
  
  dataco.finalize();
  
  return 0;
}
