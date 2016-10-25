#include <iostream>
#include <cmath>
#include <utility>
#include <assert.h>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"
#include "DataFormat/track.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"

// larelitecv
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/BoundaryEndPt.h"
#include "ThruMu/BMTrackCluster2D.h"
#include "ThruMu/BMTrackCluster3D.h"

// algos

int main( int nargs, char** argv ) {
  
  std::cout << "[Process 2D tracks to 3D tracks]" << std::endl;
  
  // configuration
  larcv::PSet cfg = larcv::CreatePSetFromFile( "tracker3d.cfg" );
  larcv::PSet bmt = cfg.get<larcv::PSet>("BMTracker3D");
  larcv::PSet sidetagger_pset  = bmt.get<larcv::PSet>("BMTSideTagger");
  
  std::string larcv_image_producer = bmt.get<std::string>("InputLArCVImages");
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  // ----------------------------------------------------------------------------------------------------
  // larcv
  dataco.add_inputfile( "output_larcv_testextbnb_track2d_10evts.root", "larcv" );
  // ----------------------------------------------------------------------------------------------------


  // configure
  dataco.configure( "tracker3d.cfg", "StorageManager", "IOManager", "BMTracker3D" );
  
  // initialize
  dataco.initialize();


  // Configure Algorithms
  // side-tagger
  larlitecv::BoundaryMuonTaggerAlgo sidetagger;
  larlitecv::ConfigBoundaryMuonTaggerAlgo sidetagger_cfg;
  sidetagger_cfg.neighborhoods  = sidetagger_pset.get< std::vector<int> >("Neighborhoods");
  sidetagger_cfg.thresholds     = sidetagger_pset.get< std::vector<float> >( "Thresholds" );
  sidetagger_cfg.emptych_thresh = sidetagger_pset.get< std::vector<float> >( "EmptyChannelThrehsold" );
  sidetagger_cfg.edge_win_wires = sidetagger_pset.get< std::vector<int> >( "EdgeWinWires" );
  sidetagger_cfg.edge_win_times = sidetagger_pset.get< std::vector<int> >( "EdgeWinTimes" );
  sidetagger_cfg.edge_win_hitthresh = sidetagger_pset.get< std::vector<float> >( "EdgeWinHitThreshold" );
  sidetagger_cfg.boundary_cluster_minpixels = sidetagger_pset.get< std::vector<int> >( "BoundaryClusterMinPixels" );
  sidetagger_cfg.boundary_cluster_radius    = sidetagger_pset.get< std::vector<float> >( "BoundaryClusterRadius" );
  sidetagger_cfg.astar_thresholds = sidetagger_pset.get< std::vector<float> >( "AStarThresholds" );
  sidetagger_cfg.astar_neighborhood = sidetagger_pset.get< std::vector<int> >( "AStarNeighborhood" );
  sidetagger.configure(sidetagger_cfg);

  // Start Event Loop
  //int nentries = dataco.get_nentries("larcv");
  //int nentries = 20;
  int nentries = 10;
  
  for (int ientry=0; ientry<nentries; ientry++) {

    std::cout << "[ENTRY " << ientry << "]" << std::endl;
    dataco.goto_entry(ientry,"larcv");

    // DATA

    // from mod channel imgs
    larcv::EventImage2D* ev_mod_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
    
    // from empty channel imgs
    larcv::EventImage2D* ev_empty_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "emptychs" );
    
    // from empty channel imgs
    larcv::EventImage2D* ev_badch_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "badchs" );
    

    const std::vector< larcv::Image2D >& imgs   = ev_mod_imgs->Image2DArray();
    const std::vector< larcv::Image2D >& badchs = ev_badch_imgs->Image2DArray();
    
    // ------------------------------------------------------------------------------------------//
    // Get 2D Track Objects

    larcv::EventPixel2D* event_tracks = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "thrumu" );
    
    std::cout << "reconstitute BMTrack2DCluster from saved data" << std::endl;
    std::vector< std::vector< larlitecv::BMTrackCluster2D > > track2d_v;
    for ( int itrack=0; itrack<(int)event_tracks->Pixel2DClusterArray(0).size(); itrack++ ) {
      std::vector< larlitecv::BMTrackCluster2D > temp_track2d_container;
      for (int p=0; p<3; p++) {
	const larcv::Pixel2DCluster& cluster = event_tracks->Pixel2DClusterArray(p).at(itrack);
	larlitecv::BMTrackCluster2D track2d;
	track2d.plane = p;
	track2d.pixelpath = cluster; // copy?
	temp_track2d_container.emplace_back( std::move(track2d) );
      }
      track2d_v.emplace_back( std::move(temp_track2d_container) );
    }
    
    // process and make 3d tracks
    std::cout << "process 2d tracks" << std::endl;
    std::vector< larlitecv::BMTrackCluster3D > tracks3d;
    sidetagger.process2Dtracks( track2d_v, ev_mod_imgs->Image2DArray(), ev_badch_imgs->Image2DArray(), tracks3d );
    
    // ------------------------------------------------------------------------------------------//
    // SAVE OUTPUT //

    // save 3D track object
    larlite::event_track* ev_tracks = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "track3d" );
    
    // convert BMTrackCluster3D
    int id = 0;
    for ( auto const& track3d : tracks3d ) {
      larlite::track lltrack;
      lltrack.set_track_id( id );
      for ( auto const& point3d : track3d.path3d ) {
	TVector3 vec( point3d[0], point3d[1], point3d[2] );
	lltrack.add_vertex( vec );
      }
      ev_tracks->emplace_back( std::move(lltrack) );
    }
    
    // go to tree
    std::cout << "save entry" << std::endl;
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
