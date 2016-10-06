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
#include "DataFormat/EventPixel2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"

// larelitecv
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/BoundaryEndPt.h"

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
  sidetagger_cfg.neighborhoods  = sidetagger_pset.get< std::vector<int> >("Neighborhoods");
  sidetagger_cfg.thresholds     = sidetagger_pset.get< std::vector<float> >( "Thresholds" );
  sidetagger_cfg.edge_win_wires = sidetagger_pset.get< std::vector<int> >( "EdgeWinWires" );
  sidetagger_cfg.edge_win_times = sidetagger_pset.get< std::vector<int> >( "EdgeWinTimes" );
  sidetagger_cfg.edge_win_hitthresh = sidetagger_pset.get< std::vector<float> >( "EdgeWinHitThreshold" );
  sidetagger_cfg.astar_thresholds = sidetagger_pset.get< std::vector<float> >( "AStarThresholds" );
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
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    std::cout << "get data: number of images=" << event_imgs->Image2DArray().size() << std::endl;
    if ( event_imgs->Image2DArray().size()==0 )
      continue;

    // ------------------------------------------------------------------------------------------//
    // SIDE TAGGER //
    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > outhits;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > end_points;
    sidetagger.searchforboundarypixels( event_imgs->Image2DArray(), outhits );
    sidetagger.clusterBoundaryPixels( event_imgs->Image2DArray(), outhits, end_points );
    if ( end_points.size()>0 ) {
      std::cout << "[[ Side Tagger End Points ]]" << std::endl;
      for (int p=0; p<3; p++) {
	std::cout << "  [Plane " << p << "]" << std::endl;
	for (int ch=0; ch<4; ch++) {
	  std::cout << "    channel " << ch << ": " << end_points.at( p*4 + ch ).size() << std::endl; 
	}
      }
    }
    else {
      end_points.resize(12);
      std::cout << "No end points?" << std::endl;
    }

    // here we take those images and do some clustering.

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

    std::vector< std::vector< larlitecv::BoundaryEndPt > > trackendpts_anode;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > trackendpts_cathode;
    for ( auto &tpc_img : event_imgs->Image2DArray() ) {

      larcv::Image2D& annode_img = stage1_annode_hits.at( (int)tpc_img.meta().plane() );
      larcv::Image2D& cathode_img = stage1_cathode_hits.at( (int)tpc_img.meta().plane() );

      std::vector< larlitecv::BoundaryEndPt > anode_ends; 
      anode_flash_tagger.findTrackEnds( opflash_containers, tpc_img, anode_ends, annode_img );
      trackendpts_anode.emplace_back( std::move(anode_ends) );

      std::vector< larlitecv::BoundaryEndPt > cathode_ends;
      cathode_flash_tagger.findTrackEnds( opflash_containers, tpc_img, cathode_ends, cathode_img );      
      trackendpts_cathode.emplace_back( std::move(cathode_ends) );
    }

    // ------------------------------------------------------------------------------------------//
    // Make tracks using end points

    
    
    
    // ------------------------------------------------------------------------------------------//
    // SAVE OUTPUT //
    
    // from side tagger
    larcv::EventImage2D* boundary_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, sidetagger_pset.get<std::string>("OutputMatchedPixelImage") );
    boundary_imgs->Emplace( std::move(outhits) );
    enum { toppt=0, botpt, uppt, dnpt, nchs };
    larcv::EventPixel2D* endpoints[4];
    endpoints[toppt] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("TopEndpoints") );
    endpoints[botpt] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("BottomEndpoints") );
    endpoints[uppt]  = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("UpstreamEndpoints") );
    endpoints[dnpt]  = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("DownstreamEndpoints") );
    std::vector< larlitecv::BoundaryEndPt > be_endpoints[nchs][3];
    for (int p=0; p<3; p++) {
      for (int ich=0; ich<nchs; ich++) {
	const std::vector< larlitecv::BoundaryEndPt >& pts = end_points.at( p*nchs + ich );
	int npts = (int)pts.size();
	for (int endpt=0; endpt<npts; endpt++) {
	  larcv::Pixel2D pixel( pts.at(endpt).w, pts.at(endpt).t );
	  pixel.Intensity( ich );
	  endpoints[ich]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
	  be_endpoints[ich][p].push_back( pts.at(endpt) );
	}
      }
    }
    
    // flash tagger
    larcv::EventImage2D* stage1_annode_imgs  = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_annode" );
    larcv::EventImage2D* stage1_cathode_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_cathode" );
    enum { flanode=0, flcathode, nflashends };
    larcv::EventPixel2D* flashends[2];
    flashends[flanode]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("AnodeEndpointProducer") );
    flashends[flcathode] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("CathodeEndpointProducer") );
    stage1_annode_imgs->Emplace( std::move(stage1_annode_hits) );
    stage1_cathode_imgs->Emplace( std::move(stage1_cathode_hits) );
    std::vector< larlitecv::BoundaryEndPt > fbe_endpoints[nflashends][3];
    for (int p=0; p<3; p++) {
      std::vector< larlitecv::BoundaryEndPt >& anode_pts   = trackendpts_anode.at(p);
      int anode_npts = (int)anode_pts.size();
      for (int ich=0; ich<anode_npts; ich++) {
	larcv::Pixel2D pixel( anode_pts.at(ich).w, anode_pts.at(ich).t );
	flashends[flanode]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
	fbe_endpoints[flanode][p].push_back( anode_pts.at(ich) );
      }
      std::vector< larlitecv::BoundaryEndPt >& cathode_pts = trackendpts_cathode.at(p);
      int cathode_npts = (int)cathode_pts.size();
      for (int ich=0; ich<cathode_npts; ich++) {
	larcv::Pixel2D pixel( cathode_pts.at(ich).w, cathode_pts.at(ich).t );
	flashends[flcathode]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
	fbe_endpoints[flcathode][p].push_back( cathode_pts.at(ich) );
      }
    }

    
    // now track clustering per plane
    larcv::EventPixel2D* event_tracks = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "thrumu" );
    for (int p=0; p<3; p++) {
      std::cout << "PLANE " << p << " TRACKING" << std::endl;
      const larcv::Image2D& img = event_imgs->Image2DArray().at(p);
      const larcv::ImageMeta& meta = img.meta();
      larcv::Image2D badchimg(meta);
      badchimg.paint(0.0);
      std::vector< larcv::Pixel2DCluster > trackcluster;
      sidetagger.makePlaneTrackCluster( img, badchimg, 
					be_endpoints[toppt][p],
					be_endpoints[botpt][p],
					be_endpoints[uppt][p],
					be_endpoints[dnpt][p],
					fbe_endpoints[flanode][p],
					fbe_endpoints[flcathode][p],
					trackcluster );
      for ( auto &track : trackcluster ) {
	std::cout << " plane=" << p << " track. length=" << track.size() << std::endl;
	event_tracks->Emplace( (larcv::PlaneID_t)p, std::move(track) );
      }
    }
    
    // save tracking
    

    // go to tree
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
