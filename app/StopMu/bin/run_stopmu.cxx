#include <iostream>
#include <sstream>

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
  dataco.add_inputfile( "output_larcv_testextbnb_001.root", "larcv" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "StopMu" );
  dataco.initialize();

  // go to some entry
  dataco.goto_entry(0, "larcv");

  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet stopmu_cfg = cfg.get<larcv::PSet>("StopMu");

  // STOP MU START PARAMETERS
  float fThreshold = 10.0;
  int rneighbor = 10;
  int cneighbor = 10;

  // for the test, we target a top-passing, stop muon
  // at time around tick 3880
  larcv::EventImage2D* imgs            = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
  larcv::EventImage2D* marked3d        = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
  
  // make the bad channel image
  larlitecv::EmptyChannelAlgo emptyalgo;
  larcv::EventChStatus* ev_status      = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
  std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status ); //< fill this out!
  std::cout << "number of bad ch imgs: " << badch_v.size() << std::endl;
  
  const std::vector<larcv::Image2D>& img_v    = imgs->Image2DArray();
  const std::vector<larcv::Image2D>& thrumu_v = marked3d->Image2DArray();
  const larcv::ImageMeta& meta = img_v.at(0).meta();

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

  larlitecv::StopMuFilterSpacePointsConfig stopmu_filter_cfg;
  stopmu_filter_cfg.duplicate_radius_cm = 10.0;
  stopmu_filter_cfg.pixel_threshold = 10.0;
  stopmu_filter_cfg.track_width_pixels = 5;
  stopmu_filter_cfg.row_tag_neighborhood = 5;
  stopmu_filter_cfg.col_tag_neighborhood = 5;

  larlitecv::StopMuFilterSpacePoints stopmu_filterpts(stopmu_filter_cfg);
  stopmu_filterpts.filterSpacePoints( ev_pixs, thrumu_v, badch_v );
  
  return 0;
}
