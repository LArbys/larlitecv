#include <iostream>
#include <cmath>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"

// larelite
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"

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
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_mcinfo_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_wire_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opdigit_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opreco_0000.root", "larlite" );

  // larcv
  dataco.add_inputfile( "data/data_samples/v05/spoon/larcv/spoon_larcv_out_0000.root", "larcv" );

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
  larlitecv::FlashMuonTaggerAlgo anode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kAnode );
  larlitecv::FlashMuonTaggerAlgo cathode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kCathode );
  larlitecv::FlashMuonTaggerConfig flashtagger_cfg;
  flashtagger_cfg.pixel_value_threshold        = flashtagger_pset.get< std::vector<float> >( "ChargeThreshold" ); // one for each plane
  flashtagger_cfg.clustering_time_neighborhood = flashtagger_pset.get< std::vector<int> >( "ClusteringTimeWindow" );
  flashtagger_cfg.clustering_wire_neighborhood = flashtagger_pset.get< std::vector<int> >( "ClusteringWireWindow" );
  flashtagger_cfg.clustering_minpoints         = flashtagger_pset.get< std::vector<int> >( "ClusteringMinPoints" );
  flashtagger_cfg.clustering_radius            = flashtagger_pset.get< std::vector<double> >( "ClusteringRadius" );
  flashtagger_cfg.endpoint_time_neighborhood   = flashtagger_pset.get< std::vector<int> >( "EndpointTimeNeighborhood" );
  flashtagger_cfg.verbosity                    = flashtagger_pset.get< int >( "Verbosity", 2 );
  flashtagger_cfg.trigger_tick                 = flashtagger_pset.get< float >( "TriggerTick", 3200.0 );
  flashtagger_cfg.usec_per_tick                = flashtagger_pset.get< float >( "MicrosecondsPerTick", 0.5 );
  flashtagger_cfg.drift_distance               = flashtagger_pset.get< float >( "DriftDistance", 250.0 );
  flashtagger_cfg.drift_velocity               = flashtagger_pset.get< float >( "DriftVelocity", 0.114 );
  anode_flash_tagger.configure( flashtagger_cfg );
  cathode_flash_tagger.configure(flashtagger_cfg);

  // Start Event Loop
  //int nentries = dataco.get_nentries("larcv");
  //int nentries = 5;
  int nentries = 1;
  
  for (int ientry=0; ientry<nentries; ientry++) {
    
    dataco.goto_entry(ientry,"larcv");

    // get images (from larcv)
    std::cout << "get data" << std::endl;
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );

    // ------------------------------------------------------------------------------------------//
    // SIDE TAGGER //
    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > outhits;
    sidetagger.searchforboundarypixels( event_imgs->Image2DArray(), outhits );

    // ------------------------------------------------------------------------------------------//
    // FLASH TAGGER //

    // create storage for new images
    std::vector< larcv::Image2D > flashtagger_hits;

    // new image for flash hits
    std::vector<larcv::Image2D> stage1_annode_hits;  // all in-time hits (non-clustered, non-edged)
    std::vector<larcv::Image2D> stage1_cathode_hits; // all in-time hits (non-clustered, non-edged)
    for ( auto &tpc_img : event_imgs->Image2DArray() ) {
      larcv::Image2D annode_img( tpc_img.meta() );
      larcv::Image2D cathode_img( tpc_img.meta() );
      annode_img.paint(0.0);
      cathode_img.paint(0.0);
      stage1_annode_hits.emplace_back( annode_img );
      stage1_cathode_hits.emplace_back( cathode_img );
    }

    std::vector<larcv::Image2D> annode_hits;
    std::vector<larcv::Image2D> cathode_hits;
    
    // loop through flash producers, get event_opflash ptrs
    std::vector< larlite::event_opflash* > opflash_containers;
    for ( auto &flashproducer : flashtagger_pset.get< std::vector<std::string> >( "OpflashProducers" ) ) {
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
      opflash_containers.push_back( opdata );
    }

    for ( auto &tpc_img : event_imgs->Image2DArray() ) {
      larcv::Image2D& annode_img = stage1_annode_hits.at( (int)tpc_img.meta().plane() );
      larcv::Image2D& cathode_img = stage1_cathode_hits.at( (int)tpc_img.meta().plane() );

      std::vector< std::vector<int> > trackendpts_anode; // we don't use them for now
      anode_flash_tagger.findTrackEnds( opflash_containers, tpc_img, trackendpts_anode, annode_img );

      std::vector< std::vector<int> > trackendpts_cathode; // we don't use them for now
      cathode_flash_tagger.findTrackEnds( opflash_containers, tpc_img, trackendpts_cathode, cathode_img );
      
    }
    
    // ------------------------------------------------------------------------------------------//
    // SAVE IMAGES //
    
    // from side tagger
    larcv::EventImage2D* boundary_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, sidetagger_pset.get<std::string>("OutputMatchedPixelImage") );
    boundary_imgs->Emplace( std::move(outhits) );    
    
    // flash tagger
    larcv::EventImage2D* stage1_annode_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_annode" );
    larcv::EventImage2D* stage1_cathode_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_cathode" );
    stage1_annode_imgs->Emplace( std::move(stage1_annode_hits) );
    stage1_cathode_imgs->Emplace( std::move(stage1_cathode_hits) );
    
    // go to tree
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
