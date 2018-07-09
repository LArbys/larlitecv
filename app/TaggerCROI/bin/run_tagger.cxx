#include <iostream>
#include <string>
#include <vector>

// larlite
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/user_info.h"
#include "larlite/UserDev/BasicTool/FhiclLite/PSet.h"
#include "larlite/UserDev/SelectionTool/LEEPreCuts/LEEPreCut.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "ANN/ANNAlgo.h"

// larlitecv
#include "GapChs/EmptyChannelAlgo.h"
#include "UnipolarHack/UnipolarHackAlgo.h"
#include "TaggerCROI/TaggerCROIAlgoConfig.h"
#include "TaggerCROI/TaggerCROIAlgo.h"
#include "TaggerCROI/TaggerCROITypes.h"
#include "TaggerCROI/PayloadWriteMethods.h"

int main(int nargs, char** argv ) {

  // pre-amble/parse arguments
  std::cout << "Cosmic Muon Tagger and Contained ROI Selection" << std::endl;

  if ( nargs!=2 ) {
  	std::cout << "usage: ./run_tagger [config file]" << std::endl;
  	return 0;
  }

  std::string cfg_file = argv[1];

  // tagger routine configuration
  larcv::PSet cfg  = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg.get<larcv::PSet>("TaggerCROI");
  
  // Print Statement
  std::cout << "Making the PSet into a Configuration file." << std::endl;
  
  // Configure application
  std::string larcv_image_producer = pset.get<std::string>("LArCVImageProducer");
  std::string larcv_chstatus_producer = pset.get<std::string>("LArCVChStatusProducer");  
  bool DeJebWires = pset.get<bool>("DeJebWires");
  float jebwiresfactor = pset.get<float>("JebWiresFactor");
  std::vector<float> emptych_thresh = pset.get< std::vector<float> >("EmptyChannelThreshold");
  std::string chstatus_datatype = pset.get<std::string>("ChStatusDataType");
  std::vector<std::string> opflash_producers = pset.get< std::vector<std::string> >( "OpflashProducers" );
  bool RunPreCuts   = pset.get<bool>("RunPrecuts");   // Run the algo
  bool ApplyPrecuts = pset.get<bool>("ApplyPrecuts"); // Apply the cuts. If precuts fail, we do not produce CROIs.
  bool RunThruMu  = pset.get<bool>("RunThruMu");
  bool RunStopMu  = pset.get<bool>("RunStopMu");
  bool RunCROI    = pset.get<bool>("RunCROI");
  bool save_thrumu_space = pset.get<bool>("SaveThruMuSpace", false);
  bool save_stopmu_space = pset.get<bool>("SaveStopMuSpace", false);
  bool save_croi_space   = pset.get<bool>("SaveCROISpace",   false);
  bool save_mc           = pset.get<bool>("SaveMC",false);
  bool skip_empty_events = pset.get<bool>("SkipEmptyEvents",false);
  bool apply_unipolar_hack = pset.get<bool>("ApplyUnipolarHack",false);

  // Configure DL PMT PreCuts
  fcllite::PSet tmp( "tmp",pset.get<larcv::PSet>("LEEPreCut").data_string() );  // convert to fcllite version of pset
  fcllite::PSet precutcfg( tmp.get<fcllite::PSet>("LEEPreCut") );
  larlite::LEEPreCut  precutalgo;
  precutalgo.configure( precutcfg ); // setup the algo
  precutalgo.initialize(); // actually does nothing  
  
  // Configure Tagger Algos
  larlitecv::TaggerCROIAlgoConfig tagger_cfg = larlitecv::TaggerCROIAlgoConfig::makeConfigFromFile( cfg_file );
  
  // Setup Input Data Coordindator
  larlitecv::DataCoordinator dataco;
  dataco.set_filelist( pset.get<std::string>("LArCVInputFilelist"),   "larcv"   );
  dataco.set_filelist( pset.get<std::string>("LArLiteInputFilelist"), "larlite" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "TaggerCROI" );
  dataco.initialize();

  // Output Data Coordinator
  larlitecv::DataCoordinator dataco_out;
  dataco_out.configure( cfg_file, "StorageManagerOut", "IOManagerOut", "TaggerCROI" );
  dataco_out.initialize();
  

  // Get run entries
  int nentries = dataco.get_nentries("larcv");
  std::cout << "NUMBER OF ENTRIES: LARCV=" << nentries << " LARLITE=" << dataco.get_nentries("larlite") << std::endl;

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

  // SETUP THE ALGOS
  larlitecv::EmptyChannelAlgo emptyalgo;
  larlitecv::UnipolarHackAlgo unihackalgo;  
  larlitecv::TaggerCROIAlgo tagger_algo( tagger_cfg );

  for (int ientry=startentry; ientry<endentry; ientry++) {

    std::cout << "Now looping over event #" << ientry << "." << std::endl;

    dataco.goto_entry(ientry,"larcv");
    int run,subrun,event;
    dataco.get_id(run,subrun,event);
    dataco_out.set_id(run,subrun,event);

    std::cout << "=====================================================" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << " Entry " << ientry << " : " << run << " " << subrun << " " << event << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // ------------------------------------------------------------------------------------------//
    // PREPARE THE INPUT
    larlitecv::InputPayload input_data = tagger_algo.loadInput( dataco );

    // -------------------------------------------------------------------------------------------//
    // RUN ALGOS
    std::cout << "[ RUN ALGOS ]" << std::endl;

    if ( RunPreCuts ) {
      precutalgo.analyze( &(dataco.get_larlite_io()) );
      // save data in larlite user data structure
      larlite::event_user* ev_precutresults = (larlite::event_user*)dataco.get_larlite_data( larlite::data::kUserInfo, "precutresults" );
      larlite::user_info precut_results;
      precut_results.store( "pass",    precutalgo.passes() );
      precut_results.store( "vetoPE",  precutalgo.vetoPE() );
      precut_results.store( "beamPE",  precutalgo.beamPE() );
      precut_results.store( "maxFrac", precutalgo.maxFrac() );
      ev_precutresults->emplace_back( std::move(precut_results) );      
    }
    else {
      // create dummy values
      larlite::event_user* ev_precutresults = (larlite::event_user*)dataco.get_larlite_data( larlite::data::kUserInfo, "precutresults" );
      larlite::user_info precut_results;      
      precut_results.store( "pass",     1   );
      precut_results.store( "vetoPE",  -1.0 );
      precut_results.store( "beamPE",  -1.0 );
      precut_results.store( "maxFrac", -1.0 );
      ev_precutresults->emplace_back( std::move(precut_results) );
    }


    if ( RunThruMu ) {
      larlitecv::ThruMuPayload thrumu_data = tagger_algo.runThruMu( input_data );
      if ( save_thrumu_space )
        thrumu_data.saveSpace();

      if ( RunStopMu ) {
        larlitecv::StopMuPayload stopmu_data = tagger_algo.runStopMu( input_data, thrumu_data );
        // we skip stopmu. we mimic instead.
        // larlitecv::StopMuPayload stopmu_data;
        // // make empty tagged pixel images
        // for ( auto const& img : input_data.img_v ) {
        //   larcv::Image2D stopmu_fake( img.meta() );
        //   stopmu_fake.paint(0);
        //   stopmu_data.stopmu_v.emplace_back( std::move(stopmu_fake) );
        // }

        if ( save_stopmu_space )
          stopmu_data.saveSpace();

        if ( RunCROI ) {
          larlitecv::CROIPayload croi_data = tagger_algo.runCROISelection( input_data, thrumu_data, stopmu_data );
          if ( save_croi_space )
            croi_data.saveSpace();
          WriteCROIPayload( croi_data, input_data, tagger_cfg, dataco_out );
        }
        WriteStopMuPayload( stopmu_data, input_data, tagger_cfg, dataco_out );
      }//end of if stopmu

      WriteThruMuPayload( thrumu_data, input_data, tagger_cfg, dataco_out );
    }

    // -------------------------------------------------------------------------------------------//
    // SAVE DATA

    WriteInputPayload( input_data, tagger_cfg, dataco_out );

    // SAVE MC DATA, IF SPECIFIED, TRANSFER FROM INPUT TO OUTPUT
    if ( save_mc && pset.get<larcv::PSet>("MCWriteConfig").get<bool>("WriteSegment") ) {
      std::cout << "WRITE SEGMENTATION IMAGE INFORMATION" << std::endl;
      larcv::PSet mcwritecfg = pset.get<larcv::PSet>("MCWriteConfig");
      larcv::EventImage2D* ev_segment  = (larcv::EventImage2D*)dataco.get_larcv_data(     larcv::kProductImage2D, mcwritecfg.get<std::string>("SegmentProducer") );
      larcv::EventImage2D* out_segment = (larcv::EventImage2D*)dataco_out.get_larcv_data( larcv::kProductImage2D, mcwritecfg.get<std::string>("SegmentProducer") );
      for ( auto const& img : ev_segment->Image2DArray() ) out_segment->Append( img );
    }

    if ( save_mc && pset.get<larcv::PSet>("MCWriteConfig").get<bool>("WriteTrackShower") ) {
      std::cout << "WRITE MC TRACK SHOWER INFORMATION" << std::endl;
      larcv::PSet mcwritecfg = pset.get<larcv::PSet>("MCWriteConfig");
      larlite::event_mctrack* event_mctrack = (larlite::event_mctrack*)dataco.get_larlite_data(     larlite::data::kMCTrack, mcwritecfg.get<std::string>("MCTrackShowerProducer") );
      larlite::event_mctrack* evout_mctrack = (larlite::event_mctrack*)dataco_out.get_larlite_data( larlite::data::kMCTrack, mcwritecfg.get<std::string>("MCTrackShowerProducer") );
      for ( auto const& track : *event_mctrack ) evout_mctrack->push_back( track  );

      larlite::event_mcshower* event_mcshower = (larlite::event_mcshower*)dataco.get_larlite_data(     larlite::data::kMCShower, mcwritecfg.get<std::string>("MCTrackShowerProducer") );
      larlite::event_mcshower* evout_mcshower = (larlite::event_mcshower*)dataco_out.get_larlite_data( larlite::data::kMCShower, mcwritecfg.get<std::string>("MCTrackShowerProducer") );
      for ( auto const& shower : *event_mcshower ) evout_mcshower->push_back( shower  );
    }

    if ( save_mc && pset.get<larcv::PSet>("MCWriteConfig").get<bool>("WriteMCROI") ) {
      std::cout << "WRITE MC TRACK SHOWER INFORMATION" << std::endl;
      larcv::PSet mcwritecfg = pset.get<larcv::PSet>("MCWriteConfig");
      larcv::EventROI* event_roi = (larcv::EventROI*)dataco.get_larcv_data(     larcv::kProductROI, mcwritecfg.get<std::string>("MCROIProducer") );
      larcv::EventROI* evout_roi = (larcv::EventROI*)dataco_out.get_larcv_data( larcv::kProductROI, mcwritecfg.get<std::string>("MCROIProducer") );
      evout_roi->Set( event_roi->ROIArray() );
    }
    
    dataco_out.save_entry();

    ann::ANNAlgo::cleanup();
  }
  
  tagger_algo.printTimeTracker( endentry-startentry );
  
  dataco.finalize();
  dataco_out.finalize();

  return 0;

}
