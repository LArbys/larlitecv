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
  std::vector<float> emptych_thresh = pset.get< std::vector<float> >("EmptryChannelThreshold");
  std::string chstatus_datatype = pset.get<std::string>("ChStatusDataType");
  std::vector<std::string> opflash_producers = pset.get< std::vector<std::string> >( "OpflashProducers" );


  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  dataco.set_filelist( pset.get<std::string>("LArCVInputFilelist"),   "larcv"   );
  dataco.set_filelist( pset.get<std::string>("LArLiteInputFilelist"), "lalrite" );

  // Configure
  dataco.configure( cfg_file, "StorageManager", "IOManager", "TaggerCROI" );
  
  // Initialize
  dataco.initialize();

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

  // Setup the TaggerCROIAlgo
  larlitecv::TaggerCROIAlgo tagger_algo( tagger_cfg );

  for (int ientry=startentry; ientry<endentry; ientry++) {

    dataco.goto_entry(ientry,"larcv");

    // ------------------------------------------------------------------------------------------//
    // PREPARE THE INPUT

    larlitecv::InputPayload input_data;

    // get images (from larcv)
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    std::cout << "get data: number of images=" << event_imgs->Image2DArray().size() << std::endl;
    if ( event_imgs->Image2DArray().size()==0 )
      continue;

    // ------------------------------------------------------------------------------------------//
    // Fill Images
    
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
    
    // ------------------------------------------------------------------------------------------//
    // LABEL EMPTY CHANNELS
    
    larlitecv::EmptyChannelAlgo emptyalgo;
    std::vector< larcv::Image2D > emptyimgs;
    for ( auto const &img : event_imgs->Image2DArray() ) {
      int p = img.meta().plane();
      larcv::Image2D emptyimg = emptyalgo.labelEmptyChannels( emptych_thresh.at(p), img );
      emptyimgs.emplace_back( emptyimg );
    }

    // ------------------------------------------------------------------------------------------//
    // LABEL BAD CHANNELS
    std::vector< larcv::Image2D > badchimgs;
    if ( chstatus_datatype=="LARLITE" ) {
      larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco.get_larlite_data( larlite::data::kChStatus, "chstatus" );
      badchimgs = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    }
    else if ( chstatus_datatype=="LARCV" ) {
      larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "chstatus" );
      badchimgs = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    }
    else {
      std::cout << "ERROR: ChStatusDataType must be either LARCV or LARLITE" << std::endl;
      return 0;
    }
    std::cout << "number of bad ch imgs: " << badchimgs.size() << std::endl;

    // ------------------------------------------------------------------------------------------//
    // LABEL GAP CHANNELS

    int maxgap = 200;
    std::vector< larcv::Image2D> gapchimgs_v = emptyalgo.findMissingBadChs( event_imgs->Image2DArray(), badchimgs, 5, maxgap );
    // combine with badchs
    for ( size_t p=0; p<badchimgs.size(); p++ ) {
      larcv::Image2D& gapchimg = gapchimgs_v.at(p);
      gapchimg += badchimgs.at(p);
    }

    // Set BadCh in input data
    for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
      input_data.badch_v.emplace_back( std::move(gapchimgs_v.at(p)) );
    }

    // -------------------------------------------------------------------------------------------//
    // COLLECT FLASH INFORMATION

    for ( auto &flashproducer : opflash_producers ) {
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
      input_data.opflashes_v.push_back( opdata );
    }   

    // -------------------------------------------------------------------------------------------//
    // RUN ALGOS

    larlitecv::ThruMuPayload thrumu_data = tagger_algo.runThruMu( input_data );

  }



  return 0;

}