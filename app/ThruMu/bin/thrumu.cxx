#include <iostream>
#include <cmath>
#include <utility>

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
    // SIDE TAGGER //
    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > outhits;
    std::vector< larcv::Image2D > realspacehits;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > end_points;
    std::vector< std::vector< larlitecv::BoundaryEndPt > > space_points;
    sidetagger.searchforboundarypixels( imgs, outhits );
    sidetagger.clusterBoundaryPixels( imgs, outhits, end_points );
    sidetagger.searchforboundarypixels3D( imgs, realspacehits );
    sidetagger.clusterBoundaryPixels3D( realspacehits, space_points );
    if ( end_points.size()>0 ) {
      std::cout << "[[ Side Tagger End Points 2D ]]" << std::endl;
      for (int p=0; p<3; p++) {
	std::cout << "  [Plane " << p << "]" << std::endl;
	for (int ch=0; ch<4; ch++) {
	  std::cout << "    channel " << ch << ": " << end_points.at( p*4 + ch ).size() << std::endl; 
	}
      }
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
    else {
      end_points.resize(12);
      std::cout << "No end points?" << std::endl;
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
    bool use_pos_flash_matching = true;
    if ( !use_pos_flash_matching ) {
//       for ( auto &tpc_img : imgs ) {
// 	larcv::Image2D& annode_img  = stage1_annode_hits.at( (int)tpc_img.meta().plane() );
// 	larcv::Image2D& cathode_img = stage1_cathode_hits.at( (int)tpc_img.meta().plane() );
// 	larcv::Image2D& imgends_img = stage1_imgends_hits.at( (int)tpc_img.meta().plane() );
	
// 	std::vector< larlitecv::BoundaryEndPt > anode_ends; 
// 	anode_flash_tagger.findTrackEnds( opflash_containers, tpc_img, anode_ends, annode_img );
// 	trackendpts_anode.emplace_back( std::move(anode_ends) );
	
// 	std::vector< larlitecv::BoundaryEndPt > cathode_ends;
// 	cathode_flash_tagger.findTrackEnds( opflash_containers, tpc_img, cathode_ends, cathode_img );      
// 	trackendpts_cathode.emplace_back( std::move(cathode_ends) );
	
// 	std::vector< larlitecv::BoundaryEndPt > imgends_ends;
// 	imgends_flash_tagger.findImageBoundaryEnds( tpc_img, imgends_ends, imgends_img );      
// 	trackendpts_imgends.emplace_back( std::move(imgends_ends) );
//       }
    }
    else {
      anode_flash_tagger.flashMatchTrackEnds( opflash_containers, imgs, trackendpts_anode, stage1_annode_hits );
      cathode_flash_tagger.flashMatchTrackEnds( opflash_containers, imgs, trackendpts_cathode, stage1_cathode_hits );
      imgends_flash_tagger.findImageTrackEnds( imgs, trackendpts_imgends, stage1_imgends_hits );
//       for (auto &tpc_img : imgs ) {
// 	larcv::Image2D& imgends_img = stage1_imgends_hits.at( (int)tpc_img.meta().plane() );
// 	std::vector< larlitecv::BoundaryEndPt > imgends_ends;
//         imgends_flash_tagger.findImageBoundaryEnds( tpc_img, imgends_ends, imgends_img );
//         trackendpts_imgends.emplace_back( std::move(imgends_ends) );
//       }
    }
    
    std::cout << "[[ Flash Tagger End Points ]]" << std::endl;
    for (int p=0; p<3; p++) {
      std::cout << "  [Plane " << p << "]" << std::endl;
      std::cout << "    anode: " << trackendpts_anode.at(p).size() << std::endl;
      std::cout << "    cathode: " << trackendpts_cathode.at(p).size() << std::endl;
      std::cout << "    imgends: " << trackendpts_imgends.at(p).size() << std::endl;
    }

//     if ( use_pos_flash_matching ) {
//       std::cout << "skipping ahead" << std::endl;
//       continue;
//     }

    // ------------------------------------------------------------------------------------------//
    // MAKE TRACKS USING COLLECTED END POINTS

    enum { toppt=0, botpt, uppt, dnpt, nchs };

    // gather boundary points: top/down/upstream/downstream
    std::vector< larlitecv::BoundaryEndPt > be_endpoints[nchs][3];
    for (int p=0; p<3; p++) {
      for (int ich=0; ich<nchs; ich++) {
	const std::vector< larlitecv::BoundaryEndPt >& pts = end_points.at( p*nchs + ich );
	int npts = (int)pts.size();
	for (int endpt=0; endpt<npts; endpt++) {
	  be_endpoints[ich][p].push_back( pts.at(endpt) );
	}
      }
    }

    // gather boundary points: anode/cathode/imageends
    enum { flanode=0, flcathode, flimgends, nflashends };
    std::vector< larlitecv::BoundaryEndPt > fbe_endpoints[nflashends][3];
    for (int p=0; p<3; p++) {
      std::vector< larlitecv::BoundaryEndPt >& anode_pts   = trackendpts_anode.at(p);
      int anode_npts = (int)anode_pts.size();
      for (int ich=0; ich<anode_npts; ich++) {
	larcv::Pixel2D pixel( anode_pts.at(ich).w, anode_pts.at(ich).t );
	fbe_endpoints[flanode][p].push_back( anode_pts.at(ich) );
      }
      std::vector< larlitecv::BoundaryEndPt >& cathode_pts = trackendpts_cathode.at(p);
      int cathode_npts = (int)cathode_pts.size();
      std::cout << "plane " << p << " cathod end points: " << cathode_npts << std::endl;
      for (int ich=0; ich<cathode_npts; ich++) {
	larcv::Pixel2D pixel( cathode_pts.at(ich).w, cathode_pts.at(ich).t );
	fbe_endpoints[flcathode][p].push_back( cathode_pts.at(ich) );
      }
      std::vector< larlitecv::BoundaryEndPt >& imgends_pts = trackendpts_imgends.at(p);
      int imgends_npts = (int)imgends_pts.size();
      for (int ich=0; ich<imgends_npts; ich++) {
	fbe_endpoints[flimgends][p].push_back( imgends_pts.at(ich) );
      }
    }
    
    // do track building
    std::vector< std::vector< larlitecv::BMTrackCluster2D > > plane_trackclusters;
    for (int p=0; p<3; p++) {
      const larcv::Image2D& img = imgs.at(p);
      const larcv::Image2D& badchimg = emptyimgs.at(p);
      std::vector< larlitecv::BMTrackCluster2D > trackcluster;
      sidetagger.makePlaneTrackCluster( img, badchimg, 
					be_endpoints[toppt][p],
					be_endpoints[botpt][p],
					be_endpoints[uppt][p],
					be_endpoints[dnpt][p],
					fbe_endpoints[flanode][p],
					fbe_endpoints[flcathode][p],
					fbe_endpoints[flimgends][p],
					trackcluster );
      plane_trackclusters.emplace_back( std::move(trackcluster) );
    }

    // mark image with pixels around the track
    std::vector< larcv::Image2D > track_marked_imgs;
    sidetagger.markImageWithTrackClusters( imgs, plane_trackclusters, track_marked_imgs );
    
    // ------------------------------------------------------------------------------------------//
    // STAGE ONE MATCHING
    //   end point match by type, being fairly conservative, looking across all three planes
    std::vector< larlitecv::BMTrackCluster3D > tracks3d;
    std::vector< std::vector<larlitecv::BMTrackCluster2D>* > input_tracks2d;
    input_tracks2d.push_back( &(plane_trackclusters.at(0)) );
    input_tracks2d.push_back( &(plane_trackclusters.at(1)) );
    input_tracks2d.push_back( &(plane_trackclusters.at(2)) );
    sidetagger.matchTracksStage1( imgs, input_tracks2d, tracks3d );
    std::cout << "Number of matched tracks returned: " << tracks3d.size() << std::endl;
    for (int i=0; i<tracks3d.size(); i++) {
      const larlitecv::BMTrackCluster3D& track3d = tracks3d.at(i);
      std::cout << " [track 3d] (" << track3d.trackidx[0] << "," << track3d.trackidx[1] << "," << track3d.trackidx[2] << ") "
		<< " start=" << track3d.tick_start << " type=" << track3d.start_type << " --> "
		<< " end=" << track3d.tick_end << " type=" << track3d.end_type << std::endl;
    }

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

    // from empty channel imgs
    larcv::EventImage2D* ev_empty_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "emptychs" );
    ev_empty_imgs->Emplace( std::move( emptyimgs ) );
    
    // from side tagger
    larcv::EventImage2D* boundary_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, sidetagger_pset.get<std::string>("OutputMatchedPixelImage") );
    boundary_imgs->Emplace( std::move(outhits) );
    larcv::EventPixel2D* endpoints[4];
    endpoints[toppt] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("TopEndpoints") );
    endpoints[botpt] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("BottomEndpoints") );
    endpoints[uppt]  = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("UpstreamEndpoints") );
    endpoints[dnpt]  = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("DownstreamEndpoints") );
    for (int p=0; p<3; p++) {
      for (int ich=0; ich<nchs; ich++) {
	const std::vector< larlitecv::BoundaryEndPt >& pts = end_points.at( p*nchs + ich );
	int npts = (int)pts.size();
	for (int endpt=0; endpt<npts; endpt++) {
	  larcv::Pixel2D pixel( pts.at(endpt).w, pts.at(endpt).t );
	  pixel.Intensity( ich );
	  endpoints[ich]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
	}
      }
    }

    // side tagger -- real space hits
    larcv::EventImage2D* realspace_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "realspacehits" );
    realspace_imgs->Emplace( std::move(realspacehits) );
    larcv::EventPixel2D* realspace_endpts[4];
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
    larcv::EventPixel2D* flashends[nflashends];
    flashends[flanode]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("AnodeEndpointProducer") );
    flashends[flcathode] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("CathodeEndpointProducer") );
    flashends[flimgends] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("ImageEndpointProducer") );
    stage1_annode_imgs->Emplace( std::move(stage1_annode_hits) );
    stage1_cathode_imgs->Emplace( std::move(stage1_cathode_hits) );
    stage1_imgends_imgs->Emplace( std::move(stage1_imgends_hits) );
    for (int p=0; p<3; p++) {
      std::vector< larlitecv::BoundaryEndPt >& anode_pts   = trackendpts_anode.at(p);
      int anode_npts = (int)anode_pts.size();
      for (int ich=0; ich<anode_npts; ich++) {
	larcv::Pixel2D pixel( anode_pts.at(ich).w, anode_pts.at(ich).t );
	flashends[flanode]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
      std::vector< larlitecv::BoundaryEndPt >& cathode_pts = trackendpts_cathode.at(p);
      int cathode_npts = (int)cathode_pts.size();
      for (int ich=0; ich<cathode_npts; ich++) {
	larcv::Pixel2D pixel( cathode_pts.at(ich).w, cathode_pts.at(ich).t );
	flashends[flcathode]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
      std::vector< larlitecv::BoundaryEndPt >& imgends_pts = trackendpts_imgends.at(p);
      int imgends_npts = (int)imgends_pts.size();
      for (int ich=0; ich<imgends_npts; ich++) {
	larcv::Pixel2D pixel( imgends_pts.at(ich).w, imgends_pts.at(ich).t );
	flashends[flimgends]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }

    
    // clustering
    std::cout << "SAVE 2D TRACKS as PIXEL2D OBJECTS" << std::endl;
    larcv::EventPixel2D* event_tracks = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "thrumu" );
    for (int p=0; p<3; p++) {
      std::cout << "Save plane " << p << " 2D tracks." << std::endl;
      std::vector< larlitecv::BMTrackCluster2D >& trackcluster = plane_trackclusters.at(p);
      for ( auto &track : trackcluster ) {
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
