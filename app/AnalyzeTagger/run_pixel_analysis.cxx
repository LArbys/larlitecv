/*
  TAGGER PIXEL ANALYSIS

  Inputs:

  larcv source file
  -----------------
  original images. (set in cfg with InputLArCVImages)
  segmentation imgs. (only if IsMC: true )
  chstatus object.
  badch images. (can specify to get from file or make from chstatus images)
  
  larlite source file
  -------------------
  mctruth from "generator"
  mctrack from "mcreco"
  mcshower from "mcreco"
  trigger from "triggersim" or "trigger" (set in cfg with TriggerProducerName)

  larcv tagger/CROI file
  ----------------------
  pixel2dclusters for thrumu tracks "thrumupixels"
  pixel2dclusters for stopmu tracks "stopmupixels"
  pixel2dclusters for untagged clusters "untaggedpixels"
  pixel2dclusters for contained ROI "croipixels"
  ROIs for contained ROIs "croi"
  boundary space points (only relevant if MC info available, IsMC:true)1

  larlite tagger/CROI file
  -------------------------
  not used

 */



#include <iostream>
#include <string>
#include <sstream>

// ROOT
#include "TFile.h"
#include "TTree.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventChStatus.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "UBWireTool/UBWireTool.h"

// larlite
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctrajectory.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/simch.h"
#include "DataFormat/trigger.h"
#include "DataFormat/track.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "TaggerTypes/dwall.h"
#include "MCTruthTools/crossingPointsAnaMethods.h"
#include "MCTruthTools/extractTruthMethods.h"



// OpenCV
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif


int main( int nargs, char** argv ) {

  TFile* output_file         = new TFile("output_file_with_croi_info.root", "RECREATE");
  TTree* tree_of_croi_values = new TTree("tree_of_croi_values", "Tree For Reading In Info About Collection Plane CROIs");

  int    croi_run;
  int    croi_subrun;
  int    croi_event;
  int    plane_num;
  double min_y;
  double max_y;
  double min_x;
  double max_x;

  tree_of_croi_values->Branch("croi_run", &croi_run, "croi_run/I");
  tree_of_croi_values->Branch("croi_subrun", &croi_subrun, "croi_subrun/I");
  tree_of_croi_values->Branch("croi_event", &croi_event, "croi_event/I");
  tree_of_croi_values->Branch("plane_num", &plane_num, "plane_num/I");
  tree_of_croi_values->Branch("min_y", &min_y, "min_y/D");
  tree_of_croi_values->Branch("max_y", &max_y, "max_y/D");
  tree_of_croi_values->Branch("min_x", &min_x, "min_x/D");
  tree_of_croi_values->Branch("max_x", &max_x, "max_x/D");

  // PIXEL ANALYSIS

  if ( nargs!=2 ) {
    std::cout << "usage: ./run_pixel_analysis [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("AnalyzeTagger");

  enum SourceTypes_t { kSource=0, kCROIfile, kNumSourceTypes };
  std::string source_param[2] = { "InputSourceFilelist", "InputCROIFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "LOADING " << source_param[isrc] << " FILES" << std::endl;
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "AnalyzeTagger" );
    dataco[isrc].initialize();
  }

  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }
  
  int nentries = dataco[kCROIfile].get_nentries("larcv");

  int user_nentries =   pset.get<int>("NumEntries",-1);
  int user_startentry = pset.get<int>("StartEntry",-1);
  int startentry = 0;
  if ( user_startentry>=0 ) {
    startentry = user_startentry;
  }
  int endentry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+startentry<nentries )
      endentry = user_nentries+startentry;
  }
  if ( startentry>=nentries ) {
    std::cout << "Starting beyond end of file. Nothing to do." << std::endl;
    return 0;
  }
  std::cout << "Start Entry: " << startentry << std::endl;
  std::cout << "End Entry: " << endentry-1 << std::endl;
  std::cout << "Buckle up!" << std::endl;
  
  for (int ientry=startentry; ientry<endentry; ientry++) {

    // load the entry
    dataco[kCROIfile].goto_entry(ientry,"larcv");
    dataco[kCROIfile].get_id(croi_run,croi_subrun,croi_event);

    std::cout << "run = " << croi_run << "." << std::endl;
    std::cout << "subrun = " << croi_subrun << "." << std::endl;
    std::cout << "event = " << croi_event << "." << std::endl;

    // sync up the other files
    for ( int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
      if ( isrc!=kCROIfile ) {
        dataco[isrc].goto_event(croi_run,croi_subrun,croi_event, "larcv");
      }
    }
    
    // get the result of the contained ROI analysis
    larcv::EventROI* ev_contained_roi = NULL;
    try {
      ev_contained_roi = (larcv::EventROI*)dataco[kCROIfile].get_larcv_data(larcv::kProductROI,"croi");
    }
    catch (...) {
      ev_contained_roi = NULL;
    }
    std::vector<larcv::ROI> containedrois_v;
    if ( ev_contained_roi ) {
      containedrois_v = ev_contained_roi->ROIArray();
      std::cout << "====ROIs===========================" << std::endl;
      for ( auto& roi : containedrois_v ) {
	std::cout << " roi: " << roi.dump();
      }
      std::cout << "===================================" << std::endl;
    }
      
    for ( auto& roi : containedrois_v ) {
      for (size_t p=0; p<3; p++ ) {
	if ( p != 2 ) continue;
	const larcv::ImageMeta& bb = roi.BB( (larcv::PlaneID_t)p );

	// Set the values that you want to fill in the tree.
	plane_num = p;
	min_y     = bb.min_y();
	max_y     = bb.max_y();
	min_x     = bb.min_x();
	max_x     = bb.max_x();

	std::cout << "Plane Number = " << p << "." << std::endl;
	std::cout << "min_y = " << min_y << "." << std::endl;
	std::cout << "max_y = " << max_y << "." << std::endl;
	std::cout << "min_x = " << min_x << "." << std::endl;
	std::cout << "max_x = " << max_x << "." << std::endl;
	

	tree_of_croi_values->Fill();

      } // End of the loop over the planes.
    } // End of the loop over the contained ROIs.
    
  }//end of entry loop

  output_file->Write();
  output_file->Close();

  return 0;

}//end of main
