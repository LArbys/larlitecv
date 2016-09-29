#include <iostream>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"

// larcv data
#include "DataFormat/EventImage2D.h"

// larelite
#include "ThruMu/BoundaryMuonTaggerAlgo.h"

// algos



int main( int nargs, char** argv ) {
  
  std::cout << "[BOUNDARY MUON TAGGER]" << std::endl;
  
  // configuration
  larcv::PSet cfg = larcv::CreatePSetFromFile( "bmt.cfg" );
  larcv::PSet bmt = cfg.get<larcv::PSet>("BoundaryMuonTagger");
  larcv::PSet sidetagger_pset  = bmt.get<larcv::PSet>("BMTSideTagger");
  larcv::PSet flashtagger_pset = bmt.get<larcv::PSet>("BMTFlashTagger");

  std::string larcv_image_producer = bmt.get<std::string>("InputLArCVImages");
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  // larlite
  dataco.add_inputfile( "data/data_samples/v05/data_extbnb_to_larcv_v00_p00/larlite_opdigit_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/data_extbnb_to_larcv_v00_p00/larlite_wire_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/data_extbnb_to_larcv_v00_p00/larlite_opreco_0000.root", "larlite" );

  // larcv
  dataco.add_inputfile( "data/data_samples/v05/data_extbnb_to_larcv_v00_p00/supera_data_0000.root", "larcv" );

  // configure
  dataco.configure( "bmt.cfg", "StorageManager", "IOManager", "BoundaryMuonTagger" );
  
  // initialize
  dataco.initialize();


  // Configure Algorithms
  // side-tagger
  larlitecv::BoundaryMuonTaggerAlgo sidetagger;
  larlitecv::ConfigBoundaryMuonTaggerAlgo sidetagger_cfg;
  sidetagger_cfg.neighborhoods = sidetagger_pset.get< std::vector<int> >("Neighborhoods");
  sidetagger_cfg.thresholds    = sidetagger_pset.get< std::vector<float> >( "Thresholds" );
  sidetagger.configure(sidetagger_cfg);

  // flash-tagger
  std::vector<std::string> opflash_producers = flashtagger_pset.get< std::vector<std::string> >( "OpflashProducers" );


  // Start Event Loop
  //int nentries = dataco.get_nentries("larcv");
  int nentries = 5;
  for (int ientry=0; ientry<nentries; ientry++) {
    
    dataco.goto_entry(ientry,"larcv");

    // get images (from larcv)
    std::cout << "get data" << std::endl;
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    larcv::EventImage2D* boundary_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, sidetagger_pset.get<std::string>("OutputMatchedPixelImage") );

    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > outhits;
    sidetagger.searchforboundarypixels( event_imgs->Image2DArray(), outhits );

    // next we run flash-based end-tagger
    //dataco.get_larlite_data(larlite::data::kOpFlash, flashtagger_pset.get<std::string>( );

    // save items

    // from side tagger
    boundary_imgs->Emplace( std::move(outhits) );    
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
