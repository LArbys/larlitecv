#include <iostream>
#include <cmath>
#include <utility>
#include <assert.h>
#include <string>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"

// larelitecv
#include "GapChs/EmptyChannelAlgo.h"
#include "GapChs/GapChProcessor.h"



int main( int nargs, char** argv ) {
  
  std::cout << "[BOUNDARY MUON TAGGER]" << std::endl;

  if (nargs<4) {
    std::cout << "usage: thrumu [config file] [larcv input list] [larlite input list]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];
  std::string input_larcv = argv[2];
  std::string input_larlite = argv[3];

  // configuration
  larcv::PSet cfg = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet bmt = cfg.get<larcv::PSet>("EstimateBadChannels");

  std::string larcv_image_producer = bmt.get<std::string>("InputLArCVImages");
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  dataco.set_filelist( input_larlite, "larlite" );
  dataco.set_filelist( input_larcv, "larcv" );

  // configure
  dataco.configure( cfg_file, "StorageManager", "IOManager", "EstimateBadChannels" );
  
  // initialize
  dataco.initialize();

  larlitecv::EmptyChannelAlgo emptyalgo;
  larlitecv::GapChProcessor gapfinder;


  int nentries = dataco.get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++ ) {

    dataco.goto_entry( ientry, "larcv" );

    // get imgs
    larcv::EventImage2D* ev_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    const std::vector<larcv::Image2D>& img_v = ev_imgs->Image2DArray();

    // make badch img from chstatus
    larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco.get_larlite_data( larlite::data::kChStatus, "chstatus" );
    std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );

    // pass them into the gap channel processor
    gapfinder.addEventImages( img_v, badch_v );
    if ( ientry>=20 )
      break;
  }
  

  std::vector<larcv::Image2D> gapch_v = gapfinder.makeGapChImage();

  return 0;
}
