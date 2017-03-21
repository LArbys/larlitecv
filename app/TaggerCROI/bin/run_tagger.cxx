#include <iostream>
#include <string>
#include <vector>

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"

// larlitecv
#include "GapChs/EmptyChannelAlgo.h"
#include "TaggerCROI/TaggerCROIAlgoConfig.h"
#include "TaggerCROI/TaggerCROIAlgo.h"
#include "TaggerCROI/TaggerCROITypes.h"
#include "TaggerCROI/PayloadWriteMethods.h"

int main(int nargs, char** argv ) {

  std::cout << "Cosmic Muon Tagger and Contained ROI Selection" << std::endl;

  if ( nargs!=2 ) {
  	std::cout << "usage: ./run_tagger [config file]" << std::endl;
  	return 0;
  }

  std::string cfg_file = argv[1];

  // configuration
  larcv::PSet cfg  = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg.get<larcv::PSet>("TaggerCROI");
  larlitecv::TaggerCROIAlgoConfig tagger_cfg = larlitecv::TaggerCROIAlgoConfig::makeConfigFromFile( cfg_file );
  std::string larcv_image_producer = pset.get<std::string>("LArCVImageProducer");
  bool DeJebWires = pset.get<bool>("DeJebWires");
  float jebwiresfactor = pset.get<float>("JebWiresFactor");
  std::vector<float> emptych_thresh = pset.get< std::vector<float> >("EmptyChannelThreshold");
  std::string chstatus_datatype = pset.get<std::string>("ChStatusDataType");
  std::vector<std::string> opflash_producers = pset.get< std::vector<std::string> >( "OpflashProducers" );
  bool RunThruMu = pset.get<bool>("RunThruMu");
  bool RunStopMu = pset.get<bool>("RunStopMu");
  bool RunCROI   = pset.get<bool>("RunCROI");    
  bool save_thrumu_space = pset.get<bool>("SaveThruMuSpace", true);
  bool save_stopmu_space = pset.get<bool>("SaveStopMuSpace", true);
  bool save_croi_space   = pset.get<bool>("SaveCROISpace", true);

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

  // SETUP THE ALGOS
  larlitecv::EmptyChannelAlgo emptyalgo;  
  larlitecv::TaggerCROIAlgo tagger_algo( tagger_cfg );

  for (int ientry=startentry; ientry<endentry; ientry++) {

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

    larlitecv::InputPayload input_data;

    // get images (from larcv)
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    if ( event_imgs->Image2DArray().size()==0) {
    	std::cout << "  Number of images=0. Skipping Entry." << std::endl;
      continue;
    }
    else if ( event_imgs->Image2DArray().size()!=3 ) {
      throw std::runtime_error("  Number of Images!=3. Weird.");
    }

    // ------------------------------------------------------------------------------------------//
    // FILL IMAGES
    
    for ( auto const &img : event_imgs->Image2DArray() ) {
      larcv::Image2D dejebbed( img );
      if ( DeJebWires ) {
      	const larcv::ImageMeta& meta = img.meta();
        for ( int col=0; col<(int)meta.cols(); col++) {
          for (int row=0; row<(int)meta.rows(); row++) {
            float val = img.pixel( row, col );
            if (meta.plane()==0 ) { 
              if ( (col>=2016 && col<=2111) || (col>=2176 && 2212) ) {
                val *= jebwiresfactor;
              }
            }
            dejebbed.set_pixel(row,col,val);
          }
        }
      }
      input_data.img_v.emplace_back( dejebbed );
    }
    if ( input_data.img_v.size()!=3 ) 
    	throw std::runtime_error("Number of Images incorrect.");
    event_imgs->clear();
    
    // ------------------------------------------------------------------------------------------//
    // LABEL EMPTY CHANNELS
    
    std::vector< larcv::Image2D > emptyimgs;
    for ( auto const &img : event_imgs->Image2DArray() ) {
      int p = img.meta().plane();
      larcv::Image2D emptyimg = emptyalgo.labelEmptyChannels( emptych_thresh.at(p), img );
      emptyimgs.emplace_back( emptyimg );
    }

    // ------------------------------------------------------------------------------------------//
    // LABEL BAD CHANNELS

    input_data.badch_v.clear();
    if ( chstatus_datatype=="LARLITE" ) {
      larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco.get_larlite_data( larlite::data::kChStatus, "chstatus" );
      input_data.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
      ev_status->clear();
    }
    else if ( chstatus_datatype=="LARCV" ) {
      larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
      input_data.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
      ev_status->clear();
    }
    else if ( chstatus_datatype=="NONE" ) {
      for ( auto const& img : input_data.img_v ) {
	larcv::Image2D badch( img.meta() );
	badch.paint(0.0);
	input_data.badch_v.emplace_back( std::move(badch) );
      }
    }
    else {
      throw std::runtime_error("ERROR: ChStatusDataType must be either LARCV or LARLITE");
    }

    if ( input_data.badch_v.size()!=3 ) {
      throw std::runtime_error("Number of Bad Channels not correct.");
    }

    // ------------------------------------------------------------------------------------------//
    // MAKE GAP CHANNEL IMAGE

    int maxgap = 200;
    std::vector< larcv::Image2D> gapchimgs_v = emptyalgo.findMissingBadChs( input_data.img_v, input_data.badch_v, 5, maxgap );
    // combine with badchs
    for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
      larcv::Image2D& gapchimg = gapchimgs_v.at(p);
      gapchimg += input_data.badch_v.at(p);
    }

    // Set BadCh in input data
    for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
      input_data.gapch_v.emplace_back( std::move(gapchimgs_v.at(p)) );
    }
    if ( input_data.gapch_v.size()!=3 ) 
    	throw std::runtime_error("Number of Gap Channels not correct.");

    // -------------------------------------------------------------------------------------------//
    // COLLECT FLASH INFORMATION

    for ( auto &flashproducer : opflash_producers ) {
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
      input_data.opflashes_v.push_back( opdata );
    }   

    // -------------------------------------------------------------------------------------------//
    // RUN ALGOS
    std::cout << "[ RUN ALGOS ]" << std::endl;
    
    if ( RunThruMu ) {
      larlitecv::ThruMuPayload thrumu_data = tagger_algo.runThruMu( input_data );
      if ( save_thrumu_space )
	thrumu_data.saveSpace();

      if ( RunStopMu ) {
        larlitecv::StopMuPayload stopmu_data = tagger_algo.runStopMu( input_data, thrumu_data );
	if ( save_stopmu_space )
	  stopmu_data.saveSpace();

	if ( RunCROI ) {
	  larlitecv::CROIPayload croi_data = tagger_algo.runCROISelection( input_data, thrumu_data, stopmu_data );
	  if ( save_croi_space )
	    croi_data.saveSpace();
	  WriteCROIPayload( croi_data, tagger_cfg, dataco_out );
	}
	
	WriteStopMuPayload( stopmu_data, tagger_cfg, dataco_out );
      }
      
      WriteThruMuPayload( thrumu_data, tagger_cfg, dataco_out );
    }


    // -------------------------------------------------------------------------------------------//
    // SAVE DATA

    WriteInputPayload( input_data, tagger_cfg, dataco_out );

    dataco_out.save_entry();

  }

  dataco.finalize();
  dataco_out.finalize();

  return 0;

}
