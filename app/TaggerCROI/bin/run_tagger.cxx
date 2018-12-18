#include <iostream>
#include <string>
#include <vector>

// larlite
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/user_info.h"
#include "DataFormat/ophit.h"
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
  bool RunPreCuts   = pset.get<bool>("RunPreCuts");   // Run the algo
  bool ApplyPreCuts = pset.get<bool>("ApplyPreCuts"); // Apply the cuts. If precuts fail, we do not produce CROIs.
  bool ApplyEndPointLimit = pset.get<bool>("ApplyEndPointLimit");
  bool RunThruMu  = pset.get<bool>("RunThruMu");
  bool RunStopMu  = pset.get<bool>("RunStopMu");
  bool RunCROI    = pset.get<bool>("RunCROI");
  bool save_thrumu_space = pset.get<bool>("SaveThruMuSpace", false);
  bool save_stopmu_space = pset.get<bool>("SaveStopMuSpace", false);
  bool save_croi_space   = pset.get<bool>("SaveCROISpace",   false);
  bool save_mc           = pset.get<bool>("SaveMC",false);
  bool skip_empty_events = pset.get<bool>("SkipEmptyEvents",false);
  bool apply_unipolar_hack = pset.get<bool>("ApplyUnipolarHack",false);
  int fEndPointLimit     = pset.get<int>("EndPointLimit");
  bool isMCC9            = pset.get<bool>("IsMCC9",false);
  std::vector<float> mcc9scale_factors(3,1.0/200.0);
  if ( isMCC9 ) {
    mcc9scale_factors = pset.get< std::vector<float> >("MCC9ScaleFactors");
  }

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
    // MCC9 Beta Hack(s)

    bool image_empty=false;
    if ( isMCC9 ) {
      // we need to scale the image down;
      for ( auto& img : input_data.img_v ) {
	float mcc9scale = mcc9scale_factors[ img.meta().plane() ];
	for (int r=0; r<(int)img.meta().rows(); r++) {
	  for (int c=0; c<(int)img.meta().cols(); c++) {
	    float pixval = img.pixel(r,c)*mcc9scale;
	    img.set_pixel(r,c,pixval);
	  }
	}
      }

      // check if image is empty
      int abovethreshpix = 0;
      for ( auto& img : input_data.img_v ) {
	for ( size_t r=0; r<img.meta().rows(); r++ ) {
	  for ( size_t c=0; c<img.meta().cols(); c++ ) {
	    if ( img.pixel(r,c)>10.0 ) abovethreshpix++;
	  }
	}
      }
      
      if ( abovethreshpix<10 )
	image_empty = true;
    }
    
    // -------------------------------------------------------------------------------------------//
    // RUN ALGOS
    std::cout << "[ RUN ALGOS ]" << std::endl;

    larlite::event_ophit* ophit_v = (larlite::event_ophit*)dataco.get_larlite_data( larlite::data::kOpHit, "ophitBeam" );
    std::cout << "ophit: " << ophit_v << std::endl;
    if ( ophit_v )
      std::cout << "ophit entries: " << ophit_v->size() << std::endl;

    bool continue_tagger = true;
    
    if ( RunPreCuts ) {
      precutalgo.analyze( &(dataco.get_larlite_io()) );
      // save data in larlite user data structure
      larlite::event_user* ev_precutresults = (larlite::event_user*)dataco_out.get_larlite_data( larlite::data::kUserInfo, "precutresults" );
      larlite::user_info precut_results;
      std::cout << "==== PRECUT RESULTS ===============" << std::endl;
      std::cout << " PASS: "   << precutalgo.passes() << std::endl;
      std::cout << " vetope:   " << precutalgo.vetoPE() << std::endl;
      std::cout << " beampe:   " << precutalgo.beamPE() << std::endl;
      std::cout << " maxfrac:  " << precutalgo.maxFrac() << std::endl;
      std::cout << " beamTick: " << precutalgo.beamFirstTick() << std::endl;
      std::cout << " vetoTick: " << precutalgo.vetoFirstTick() << std::endl;
      std::cout << "===================================" << std::endl;
	
      precut_results.store( "pass",    precutalgo.passes() );
      precut_results.store( "vetoPE",  precutalgo.vetoPE() );
      precut_results.store( "beamPE",  precutalgo.beamPE() );
      precut_results.store( "maxFrac", precutalgo.maxFrac() );
      precut_results.store( "beamFirstTick", precutalgo.beamFirstTick() );
      precut_results.store( "vetoFirstTick", precutalgo.vetoFirstTick() );      
      ev_precutresults->emplace_back( std::move(precut_results) );
    }
    else {
      // create dummy values
      larlite::event_user* ev_precutresults = (larlite::event_user*)dataco_out.get_larlite_data( larlite::data::kUserInfo, "precutresults" );
      larlite::user_info precut_results;      
      precut_results.store( "pass",     1   );
      precut_results.store( "vetoPE",  -1.0 );
      precut_results.store( "beamPE",  -1.0 );
      precut_results.store( "maxFrac", -1.0 );
      precut_results.store( "beamFirstTick", -1 );
      precut_results.store( "vetoFirstTick", -1 );
      ev_precutresults->emplace_back( std::move(precut_results) );
    }

    if ( (ApplyPreCuts && precutalgo.passes()==0) || image_empty ) {

      std::cout << "[PRECUTS FAILED (or) IMAGE IS EMPTY]" << std::endl;
      
      // precuts failed.
      // we want to create an empty entry still
      larlitecv::ThruMuPayload thrumu_data;
      larlitecv::StopMuPayload stopmu_data;
      larlitecv::CROIPayload   croi_data;
      bool fillempty = true;
      WriteCROIPayload(   croi_data,   input_data, tagger_cfg, dataco_out, fillempty );
      WriteStopMuPayload( stopmu_data, input_data, tagger_cfg, dataco_out, fillempty );
      WriteThruMuPayload( thrumu_data, input_data, tagger_cfg, dataco_out, fillempty );
      continue_tagger = false;
    }

    // if we are not applying the precuts (only running the algo)
    // or if we passed the precuts (yay!)
    // then run the rest of the tagger code

    // Boundary End Points
    int num_boundary_points = 0;
    larlitecv::ThruMuPayload thrumu_data;
    if ( continue_tagger ) {

      // first start off by findnig end points
      larlitecv::ThruMuPayload endpoint_thrumu_data = tagger_algo.runBoundaryPointFinder( input_data );
      int num_boundary_points = 0;
      num_boundary_points += endpoint_thrumu_data.side_filtered_v.size();
      num_boundary_points += endpoint_thrumu_data.anode_filtered_v.size();
      num_boundary_points += endpoint_thrumu_data.cathode_filtered_v.size();
      num_boundary_points += endpoint_thrumu_data.imgends_filtered_v.size();
            
      if ( ApplyEndPointLimit && num_boundary_points>=fEndPointLimit ) {
	// we will stop the tagger here, so write empty data
	larlitecv::StopMuPayload stopmu_data;
	larlitecv::CROIPayload   croi_data;
	bool fillempty = true;
	WriteCROIPayload(   croi_data,   input_data, tagger_cfg, dataco_out, fillempty );
	WriteStopMuPayload( stopmu_data, input_data, tagger_cfg, dataco_out, fillempty );
	WriteThruMuPayload( thrumu_data, input_data, tagger_cfg, dataco_out, fillempty );
	continue_tagger = false;
      }

      std::swap( thrumu_data, endpoint_thrumu_data ); // use swap to avoid a copy
    }
    
    larlite::event_user* ev_endpointresults = (larlite::event_user*)dataco_out.get_larlite_data( larlite::data::kUserInfo, "endpointresults" );
    larlite::user_info endpoint_results;
    endpoint_results.store("numendpoints", num_boundary_points);
    ev_endpointresults->emplace_back( std::move(endpoint_results) );
    
    // ThruMu Tracker/StopMu/Contained + FlashAnalysis
    if ( continue_tagger ) {
      
      if ( RunThruMu ) {
	tagger_algo.runThruMu( input_data, thrumu_data );
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
    }// if dont apply precuts or precuts pass

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
