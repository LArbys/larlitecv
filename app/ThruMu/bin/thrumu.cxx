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
#include "../../larlite/core/DataFormat/chstatus.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"

// larelitecv
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/BoundaryEndPt.h"
#include "ThruMu/EmptyChannelAlgo.h"

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

  /*
  // ----------------------------------------------------------------------------------------------------
  // SPOON (MC)
  // larlite
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_mcinfo_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_wire_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opdigit_0000.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opreco_0000.root", "larlite" );

  // larcv
  dataco.add_inputfile( "data/data_samples/v05/spoon/larcv/spoon_larcv_out_0000.root", "larcv" );
  // ----------------------------------------------------------------------------------------------------
  */

  // ----------------------------------------------------------------------------------------------------
  // Example EXTBNB
  // larlite
  dataco.add_inputfile( "data/data_samples/v05/extbnb/larlite_opdigit_20161020_192212_034523.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/extbnb/larlite_opreco_20161020_192212_037150.root", "larlite" );
  dataco.add_inputfile( "data/data_samples/v05/extbnb/larlite_wire_20161020_192212_031437.root", "larlite" );
  // larcv
  dataco.add_inputfile( "data/data_samples/v05/extbnb/supera_extbnb.root", "larcv" );
  // ----------------------------------------------------------------------------------------------------

  // larcv
  

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
  sidetagger_cfg.emptych_thresh = sidetagger_pset.get< std::vector<float> >( "EmptyChannelThrehsold" );
  sidetagger_cfg.edge_win_wires = sidetagger_pset.get< std::vector<int> >( "EdgeWinWires" );
  sidetagger_cfg.edge_win_times = sidetagger_pset.get< std::vector<int> >( "EdgeWinTimes" );
  sidetagger_cfg.edge_win_hitthresh = sidetagger_pset.get< std::vector<float> >( "EdgeWinHitThreshold" );
  sidetagger_cfg.boundary_cluster_minpixels = sidetagger_pset.get< std::vector<int> >( "BoundaryClusterMinPixels" );
  sidetagger_cfg.boundary_cluster_radius    = sidetagger_pset.get< std::vector<float> >( "BoundaryClusterRadius" );
  sidetagger_cfg.astar_thresholds = sidetagger_pset.get< std::vector<float> >( "AStarThresholds" );
  sidetagger_cfg.astar_neighborhood = sidetagger_pset.get< std::vector<int> >( "AStarNeighborhood" );
  sidetagger.configure(sidetagger_cfg);

  // flash-tagger
  larlitecv::FlashMuonTaggerAlgo anode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kAnode );
  larlitecv::FlashMuonTaggerAlgo cathode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kCathode );
  larlitecv::FlashMuonTaggerAlgo imgends_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kOutOfImage );
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
  imgends_flash_tagger.configure(flashtagger_cfg);

  // Start Event Loop
  //int nentries = dataco.get_nentries("larcv");
  //int nentries = 20;
  int nentries = 1;
  
  for (int ientry=0; ientry<nentries; ientry++) {
    
    dataco.goto_entry(ientry,"larcv");

    // get images (from larcv)
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    std::cout << "get data: number of images=" << event_imgs->Image2DArray().size() << std::endl;
    if ( event_imgs->Image2DArray().size()==0 )
      continue;

    // ------------------------------------------------------------------------------------------//
    // CORRECT JEB WIRES
    std::vector< larcv::Image2D > imgs;
    for ( auto const &img : event_imgs->Image2DArray() ) {
      const larcv::ImageMeta& meta = img.meta();
      larcv::Image2D dejebbed(meta);
      for ( int col=0; col<meta.cols(); col++) {
	for (int row=0; row<meta.rows(); row++) {
	  float val = img.pixel( row, col );
	  if (meta.plane()==0 ) { 
	    if ( (col>=2016 && col<=2111) || (col>=2176 && 2212) ) {
	      val *= 5.0;
	    }
	  }
	  dejebbed.set_pixel(row,col,val);
	}
      }
      imgs.emplace_back( dejebbed );
    }
    
    // ------------------------------------------------------------------------------------------//
    // LABEL EMPTY CHANNELS
    
    larlitecv::EmptyChannelAlgo emptyalgo;
    std::vector< larcv::Image2D > emptyimgs;
    for ( auto const &img : event_imgs->Image2DArray() ) {
      int p = img.meta().plane();
      larcv::Image2D emptyimg = emptyalgo.labelEmptyChannels( sidetagger_cfg.emptych_thresh.at(p), img );
      emptyimgs.emplace_back( emptyimg );
    }

    // ------------------------------------------------------------------------------------------//
    // LABEL BAD CHANNELS
    
    larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco.get_larlite_data( larlite::data::kChStatus, "chstatus" );
    std::cout << "ch status planes: " << ev_status->size() << std::endl;
    std::vector< larcv::Image2D > badchimgs = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    std::cout << "number of bad ch imgs: " << badchimgs.size() << std::endl;

    // ------------------------------------------------------------------------------------------//
    // SIDE TAGGER //
    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > outhits;
    std::vector< larcv::Image2D > boundarypixels;
    std::vector< larcv::Image2D > realspacehits;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > space_points;
    sidetagger.searchforboundarypixels3D( imgs, badchimgs, space_points, boundarypixels, realspacehits );
    if ( space_points.size()>0 ) {
      std::cout << "[[ Side Tagger End Points 3D=" << space_points.size() << " ]]" << std::endl;
      int ncrossings[4] = {0};
      for (int i=0; i<space_points.size(); i++) {
	ncrossings[ space_points.at(i).at(0).type ]++;
      }
      std::cout << "  top=" << ncrossings[0] << std::endl;
      std::cout << "  bottom=" << ncrossings[1] << std::endl;
      std::cout << "  upstream=" << ncrossings[2] << std::endl;
      std::cout << "  downstream=" << ncrossings[3] << std::endl;
    }

    // here we take those images and do some clustering.

    // ------------------------------------------------------------------------------------------//
    // FLASH TAGGER //
    std::cout << "[Run Flash Tagger]" << std::endl;

    // create storage for new images
    std::vector< larcv::Image2D > flashtagger_hits;

    // new image for flash hits
    std::vector<larcv::Image2D> stage1_annode_hits;  // all in-time hits (non-clustered, non-edged)
    std::vector<larcv::Image2D> stage1_cathode_hits; // all in-time hits (non-clustered, non-edged)
    std::vector<larcv::Image2D> stage1_imgends_hits; // all in-time hits (non-clustered, non-edged)
    for ( auto &tpc_img : imgs ) {
      larcv::Image2D annode_img( tpc_img.meta() );
      larcv::Image2D cathode_img( tpc_img.meta() );
      larcv::Image2D imgends_img( tpc_img.meta() );
      annode_img.paint(0.0);
      cathode_img.paint(0.0);
      imgends_img.paint(0.0);
      stage1_annode_hits.emplace_back( annode_img );
      stage1_cathode_hits.emplace_back( cathode_img );
      stage1_imgends_hits.emplace_back( imgends_img );
    }

    std::vector<larcv::Image2D> annode_hits;
    std::vector<larcv::Image2D> cathode_hits;
    std::vector<larcv::Image2D> imgends_hits;
    
    // loop through flash producers, get event_opflash ptrs
    std::vector< larlite::event_opflash* > opflash_containers;
    for ( auto &flashproducer : flashtagger_pset.get< std::vector<std::string> >( "OpflashProducers" ) ) {
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
      opflash_containers.push_back( opdata );
    }

    std::vector< std::vector< larlitecv::BoundaryEndPt > > trackendpts_anode;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > trackendpts_cathode;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > trackendpts_imgends;
    anode_flash_tagger.flashMatchTrackEnds( opflash_containers, imgs, trackendpts_anode, stage1_annode_hits );
    cathode_flash_tagger.flashMatchTrackEnds( opflash_containers, imgs, trackendpts_cathode, stage1_cathode_hits );
    imgends_flash_tagger.findImageTrackEnds( imgs, trackendpts_imgends, stage1_imgends_hits );
    
    std::cout << "[[ Flash Tagger End Points ]]" << std::endl;
    std::cout << "    anode: " << trackendpts_anode.size() << std::endl;
    std::cout << "    cathode: " << trackendpts_cathode.size() << std::endl;
    std::cout << "    imgends: " << trackendpts_imgends.size() << std::endl;

//     if ( use_pos_flash_matching ) {
//       std::cout << "skipping ahead" << std::endl;
//       continue;
//     }

    // ------------------------------------------------------------------------------------------//
    // MAKE TRACKS USING COLLECTED END POINTS

    std::vector< const std::vector< larlitecv::BoundaryEndPt >* > all_endpoints;
    
    // gather endpoints from space points
    for (int isp=0; isp<space_points.size(); isp++) {
      const std::vector< larlitecv::BoundaryEndPt >* pts = &(space_points.at( isp ));
      all_endpoints.push_back( pts );
    }

    // gather boundary points: anode/cathode/imageends
    for (int isp=0; isp<trackendpts_anode.size(); isp++) {
      const std::vector< larlitecv::BoundaryEndPt >* pts = &(trackendpts_anode.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<trackendpts_cathode.size(); isp++) {
      const std::vector< larlitecv::BoundaryEndPt >* pts = &(trackendpts_cathode.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<trackendpts_imgends.size(); isp++) {
      const std::vector< larlitecv::BoundaryEndPt >* pts = &(trackendpts_imgends.at(isp));
      all_endpoints.push_back( pts );
    }
    

    // do track building
    std::vector< std::vector< larlitecv::BMTrackCluster2D > > trackclusters;
    std::vector< larlitecv::BMTrackCluster2D > trackcluster;
    sidetagger.makeTrackClusters3D( imgs, badchimgs, all_endpoints, trackclusters );


    // mark image with pixels around the track
    std::vector< larcv::Image2D > track_marked_imgs;
    //sidetagger.markImageWithTrackClusters( imgs, plane_trackclusters, track_marked_imgs );
    
    // ------------------------------------------------------------------------------------------//
    // STAGE ONE MATCHING
    //   end point match by type, being fairly conservative, looking across all three planes
//     std::vector< larlitecv::BMTrackCluster3D > tracks3d;
//     std::vector< std::vector<larlitecv::BMTrackCluster2D>* > input_tracks2d;
//     input_tracks2d.push_back( &(plane_trackclusters.at(0)) );
//     input_tracks2d.push_back( &(plane_trackclusters.at(1)) );
//     input_tracks2d.push_back( &(plane_trackclusters.at(2)) );
//     //sidetagger.matchTracksStage1( imgs, input_tracks2d, tracks3d );
//     std::cout << "Number of matched tracks returned: " << tracks3d.size() << std::endl;
//     for (int i=0; i<tracks3d.size(); i++) {
//       const larlitecv::BMTrackCluster3D& track3d = tracks3d.at(i);
//       std::cout << " [track 3d] (" << track3d.trackidx[0] << "," << track3d.trackidx[1] << "," << track3d.trackidx[2] << ") "
// 		<< " start=" << track3d.tick_start << " type=" << track3d.start_type << " --> "
// 		<< " end=" << track3d.tick_end << " type=" << track3d.end_type << std::endl;
//     }

    // Stage 2 matching
    //  end point match by type again, but look over two planes, and going into the third to match

    // Stage 3 matching
    //  maybe it doesn't come to this?  Go into the existing ones and break up tracks at intersections and kinks

    // Stage 4 matching
    //  despearation

    // Stage 5
    //  acceptance

    // ------------------------------------------------------------------------------------------//
    // SAVE OUTPUT //

    // from mod channel imgs
    larcv::EventImage2D* ev_mod_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
    ev_mod_imgs->Emplace( std::move( imgs ) );

    // from empty channel imgs
    larcv::EventImage2D* ev_empty_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "emptychs" );
    ev_empty_imgs->Emplace( std::move( emptyimgs ) );

    // from empty channel imgs
    larcv::EventImage2D* ev_badch_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "badchs" );
    ev_badch_imgs->Emplace( std::move( badchimgs ) );
    
    // side tagger -- real space hits
    larcv::EventImage2D* realspace_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "realspacehits" );
    realspace_imgs->Emplace( std::move(realspacehits) );
    larcv::EventImage2D* boundarypixels_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "boundarypixels" );
    boundarypixels_imgs->Emplace( std::move(boundarypixels) );
    enum { toppt=0, botpt, uppt, dnpt, nendpts };
    larcv::EventPixel2D* realspace_endpts[nendpts];
    realspace_endpts[toppt] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("TopSpacePts") );
    realspace_endpts[botpt] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("BottomSpacePts") );
    realspace_endpts[uppt]  = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("UpstreamSpacePts") );
    realspace_endpts[dnpt]  = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("DownstreamSpacePts") );
    for (int i=0; i<space_points.size(); i++) {
      const std::vector< larlitecv::BoundaryEndPt >& sp_v = space_points.at(i);
      int sptype = (int)sp_v.at(0).type;
      for (int p=0; p<3; p++) {
	const larlitecv::BoundaryEndPt& sp = sp_v.at(p);
	larcv::Pixel2D pixel( sp.w, sp.t );
	pixel.Intensity( sptype );
	realspace_endpts[sptype]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }

    // flash tagger
    larcv::EventImage2D* stage1_annode_imgs  = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_annode" );
    larcv::EventImage2D* stage1_cathode_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_cathode" );
    larcv::EventImage2D* stage1_imgends_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stage1_imgends" );
    stage1_annode_imgs->Emplace( std::move(stage1_annode_hits) );
    stage1_cathode_imgs->Emplace( std::move(stage1_cathode_hits) );
    stage1_imgends_imgs->Emplace( std::move(stage1_imgends_hits) );

    enum { flanode=0, flcathode, flimgends, nflashends };
    larcv::EventPixel2D* flashends[nflashends];
    flashends[flanode]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("AnodeEndpointProducer") );
    flashends[flcathode] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("CathodeEndpointProducer") );
    flashends[flimgends] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("ImageEndpointProducer") );        
    for (int i=0; i<(int)trackendpts_anode.size(); i++) {
      const std::vector< larlitecv::BoundaryEndPt >& anode_pts   = trackendpts_anode.at(i);
      for (int p=0; p<3; p++) {
	larcv::Pixel2D pixel( anode_pts.at(p).w, anode_pts.at(p).t );
	flashends[flanode]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }
    for (int i=0; i<(int)trackendpts_cathode.size(); i++) {
      const std::vector< larlitecv::BoundaryEndPt >& cathode_pts   = trackendpts_cathode.at(i);
      for (int p=0; p<3; p++) {
	larcv::Pixel2D pixel( cathode_pts.at(p).w, cathode_pts.at(p).t );
	flashends[flcathode]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }
    for (int i=0; i<(int)trackendpts_imgends.size(); i++) {
      const std::vector< larlitecv::BoundaryEndPt >& imgends_pts   = trackendpts_imgends.at(i);
      for (int p=0; p<3; p++) {
	larcv::Pixel2D pixel( imgends_pts.at(p).w, imgends_pts.at(p).t );
	flashends[flimgends]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }
    
    // clustering
    std::cout << "SAVE 2D TRACKS as PIXEL2D OBJECTS" << std::endl;
    larcv::EventPixel2D* event_tracks = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "thrumu" );
    for (int i=0; i<(int)trackclusters.size(); i++) {
      std::vector< larlitecv::BMTrackCluster2D >& trackcluster = trackclusters.at(i);
      std::cout << "save track cluster" << std::endl;
      for (int p=0; p<3; p++) {
	larlitecv::BMTrackCluster2D& track = trackcluster.at(p);
	larcv::Pixel2DCluster cluster;
	std::swap( cluster, track.pixelpath );
	std::cout << " plane=" << p << " track. length=" << cluster.size() << std::endl;
	event_tracks->Emplace( (larcv::PlaneID_t)p, std::move(cluster) );
      }
    }
    
    // track cluster img
    larcv::EventImage2D* ev_trackcluster_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "trackcluster" );
    ev_trackcluster_imgs->Emplace( std::move( track_marked_imgs ) );

    // go to tree
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
