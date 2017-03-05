#include <iostream>
#include <sstream>
#include <exception>

// ROOT
#include "TVector3.h"

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
#include "DataFormat/track.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/Pixel2DCluster.h"
#include "CVUtil/CVUtil.h"

// larlitecv
#include "GapChs/EmptyChannelAlgo.h"

// larlite/stopmu
#include "StopMuAlgoTypes.h"
#include "StopMuFilterSpacePoints.h"
#include "StopMuStart.h"
#include "StopMuClusterConfig.h"
#include "StopMuCluster.h"


int main( int nargs, char** argv ) {
  std::cout << "Test the stop mu tracker." << std::endl;

  if ( nargs!=4 && nargs!=2) {
    std::cout << "Usage:" << std::endl;
    std::cout << "./stopmu [config file] [larcv input filelist] [larlite input filelist]" << std::endl;
    std::cout << " -- or -- " << std::endl;
    std::cout << "./stopmu [configfile with input filelists specified]" << std::endl;
    return 0;
  }
  

  // config file
  std::string cfg_file         = argv[1];
 
  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet stopmu_cfg = cfg.get<larcv::PSet>("StopMu");

  // --------------------------------------------
  // load data
  std::string flist_larcv;  
  std::string flist_larlite;
  if ( nargs==2 ) {
    flist_larcv   = stopmu_cfg.get<std::string>("LArCVInputList");
    flist_larlite = stopmu_cfg.get<std::string>("LArLiteInputList");
  }
  else if ( nargs==4 ) {
    flist_larcv = argv[2];
    flist_larlite = argv[3];
  }
  else {
    std::cout << "wrong number of arguments." << std::endl;
    return 0;
  }
  
  // data coordinator setup
  larlitecv::DataCoordinator dataco;
  dataco.set_filelist( flist_larcv,   "larcv" );
  dataco.set_filelist( flist_larlite, "larlite" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "StopMu" );
  dataco.initialize();

  // --------------------------------------------
  // Instatiate Algos
    
  // filter out end points we are unlikely to be interested in: duplicates and those already used as thru-mu end points
  larcv::PSet stopmu_filter_pset = stopmu_cfg.get<larcv::PSet>( "StopMuSpacePointsFilter" );
  larlitecv::StopMuFilterSpacePointsConfig stopmu_filter_cfg = larlitecv::MakeStopMuFilterSpacePointsConfigFromPSet( stopmu_filter_pset );
  larlitecv::StopMuFilterSpacePoints stopmu_filterpts(stopmu_filter_cfg);

  larlitecv::StopMuClusterConfig stopmu_cluster_config = larlitecv::makeStopMuClusterConfigFromPSet( stopmu_cfg.get<larcv::PSet>("StopMuCluster") );
  larlitecv::StopMuCluster smcluster( stopmu_cluster_config);

  // start point direction
  larlitecv::StopMuStart start_finder_algo;
  start_finder_algo.setVerbose(1);
  float fThreshold = 10.0;
  int rneighbor = 10;
  int cneighbor = 10;

  // tagging parameters
  int tagged_stopmu_pixelradius = 2;

  int nentries = dataco.get_nentries("larcv");
  int user_nentries =   stopmu_cfg.get<int>("NumEntries",-1);
  int user_startentry = stopmu_cfg.get<int>("StartEntry",-1);
  int startentry = 0;
  if ( user_startentry>=0 ) {
    startentry = user_startentry;
  }
  int endentry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+startentry<nentries )
      endentry = user_nentries+startentry;
  }
  
  for (int ientry=startentry; ientry<endentry; ientry++) {

    std::cout << "=====================================================" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << " Entry " << ientry << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // go to some entry
    dataco.goto_entry(ientry, "larcv");
    int run,subrun,event;
    dataco.get_id(run,subrun,event);
    
    larcv::EventImage2D* imgs            = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
    larcv::EventImage2D* marked3d        = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
    larcv::EventROI* rois                = (larcv::EventROI*)dataco.get_larcv_data( larcv::kProductROI, "tpc" );
    
    // make the bad channel image
    larlitecv::EmptyChannelAlgo emptyalgo;
    larcv::EventChStatus* ev_status      = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
    std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status ); //< fill this out!
    std::cout << "number of bad ch imgs: " << badch_v.size() << std::endl;

    int maxgap = 200;
    std::vector< larcv::Image2D> gapchimgs_v = emptyalgo.findMissingBadChs( imgs->Image2DArray(), badch_v, 5, maxgap );
    // combine with badchs
    for ( size_t p=0; p<badch_v.size(); p++ ) {
      larcv::Image2D& gapchimg = gapchimgs_v.at(p);
      gapchimg += badch_v.at(p);
    }
    
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
  
    // output: stopmu-tagged pixels
    std::vector<larcv::Image2D> stopmu_v;
    for (size_t p=0; p<img_v.size(); p++) {
      larcv::Image2D stopmu_img( img_v.at(p).meta() );
      stopmu_img.paint(0);
      stopmu_v.emplace_back( std::move(stopmu_img) );
    }
    // output: pixel clusters for eah stopmu track
    larcv::EventPixel2D* ev_stopmu_pixels = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "stopmupixels" );
    // output: 3D trajectory points from stopmu tracker
    larlite::event_track* ev_stopmu_tracks = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "stopmutracks" );

    // --------------------------------------------
    // Algo Prep
    std::vector< std::vector< const larcv::Pixel2D* > > stopmu_candidate_endpts = stopmu_filterpts.filterSpacePoints( ev_pixs, thrumu_v, badch_v );
    std::cout << " Number of candidate stop-mu start points: " << stopmu_candidate_endpts.size() << std::endl;

    std::stringstream ss;
    ss << "smcluster_" << ientry << "_r" << run << "_s" << subrun << "_e" << event;

    smcluster.setOpenCVImageStemName( ss.str() );
    smcluster.findStopMuTracks( img_v, gapchimgs_v, thrumu_v, stopmu_candidate_endpts );

    larcv::EventImage2D* stopmu_eventimgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stopmu" );
    stopmu_eventimgs->Emplace( std::move(stopmu_v) );
    
    dataco.save_entry();

    //if ( ientry>=10 )
    //break;
    
  }// loop over entries
  
  dataco.finalize();
  
  return 0;
}
