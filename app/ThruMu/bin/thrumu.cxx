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
#include "DataFormat/opflash.h"
#include "DataFormat/track.h"
#include "../../larlite/core/DataFormat/chstatus.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"

// larelitecv
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/BoundarySpacePoint.h"
#include "ThruMu/BoundaryEndPt.h"
#include "ThruMu/EndPointFilter.h"
#include "GapChs/EmptyChannelAlgo.h"



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
  larcv::PSet bmt = cfg.get<larcv::PSet>("BoundaryMuonTagger");
  larcv::PSet sidetagger_pset  = bmt.get<larcv::PSet>("BMTSideTagger");
  larcv::PSet flashtagger_pset = bmt.get<larcv::PSet>("BMTFlashTagger");
  larcv::PSet thrumutrack_pset = bmt.get<larcv::PSet>("ThruMuTracker");

  std::string larcv_image_producer = bmt.get<std::string>("InputLArCVImages");
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  dataco.set_filelist( input_larlite, "larlite" );
  dataco.set_filelist( input_larcv, "larcv" );

  // configure
  dataco.configure( cfg_file, "StorageManager", "IOManager", "BoundaryMuonTagger" );
  
  // initialize
  dataco.initialize();

  // Configure Algorithms
  // side-tagger
  larlitecv::BoundaryMuonTaggerAlgo sidetagger;
  larlitecv::ConfigBoundaryMuonTaggerAlgo sidetagger_cfg = larlitecv::MakeConfigBoundaryMuonTaggerAlgoFromPSet( sidetagger_pset );
  sidetagger.configure(sidetagger_cfg);

  sidetagger_cfg.print();

  // flash-tagger
  larlitecv::FlashMuonTaggerAlgo anode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kAnode );
  larlitecv::FlashMuonTaggerAlgo cathode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kCathode );
  larlitecv::FlashMuonTaggerAlgo imgends_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kOutOfImage );
  larlitecv::FlashMuonTaggerConfig flashtagger_cfg = larlitecv::MakeFlashMuonTaggerConfigFromPSet( flashtagger_pset );
  anode_flash_tagger.configure( flashtagger_cfg );
  cathode_flash_tagger.configure(flashtagger_cfg);
  imgends_flash_tagger.configure(flashtagger_cfg);

  // end point filter algo
  larlitecv::EndPointFilter endptfilter;

  // other parameters 
  int max_endpts_to_process = bmt.get<int>("MaxEndPointsToProcess"); // safety valve
  int jebwiresfactor = bmt.get<float>("JebWiresFactor",5.0);

  // Start Event Loop
  int nentries = dataco.get_nentries("larcv");
  int user_nentries =   bmt.get<int>("NumEntries",-1);
  int user_startentry = bmt.get<int>("StartEntry",-1);
  int startentry = 0;
  if ( user_startentry>=0 ) {
    startentry = user_startentry;
  }
  int endentry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+startentry<nentries )
      endentry = user_nentries+startentry;
  }
  
  for (int ientry=startentry; ientry<endentry; ientry++) {
    
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
    std::vector< larcv::Image2D > badchimgs;
    if ( bmt.get<std::string>("ChStatusDataType")=="LARLITE" ) {
      larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco.get_larlite_data( larlite::data::kChStatus, "chstatus" );
      badchimgs = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    }
    else if (bmt.get<std::string>("ChStatusDataType")=="LARCV") {
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
    //int maxgap = 16;
    int maxgap = 200;
    std::vector< larcv::Image2D> gapchimgs_v = emptyalgo.findMissingBadChs( event_imgs->Image2DArray(), badchimgs, 5, maxgap );
    // combine with badchs
    for ( size_t p=0; p<badchimgs.size(); p++ ) {
      larcv::Image2D& gapchimg = gapchimgs_v.at(p);
      gapchimg += badchimgs.at(p);
    }

    // ------------------------------------------------------------------------------------------//
    // SIDE TAGGER //
    // first we run the boundary muon tagger algo that matches top/side/upstream/downstream
    std::cout << "search for boundary pixels" << std::endl;
    std::vector< larcv::Image2D > boundarypixels; // hits in real space cooridinaes
    std::vector< larcv::Image2D > realspacehits;  // 
    std::vector< larlitecv::BoundarySpacePoint > space_points;
    std::vector< std::vector<int> > boundary_endpt_usedtag(7); // [top,bot,up,down,anode,cathode,imgend]
    sidetagger.searchforboundarypixels3D( imgs, badchimgs, space_points, boundarypixels, realspacehits );
    if ( space_points.size()>0 ) {
      int nsidepts[4] = {0};
      std::cout << "[[ Side Tagger End Points 3D=" << space_points.size() << " ]]" << std::endl;
      for (int i=0; i<(int)space_points.size(); i++) {
        nsidepts[ space_points.at(i).at(0).type ]++;
      }
      std::cout << "  top=" << nsidepts[0] << std::endl;
      std::cout << "  bottom=" << nsidepts[1] << std::endl;
      std::cout << "  upstream=" << nsidepts[2] << std::endl;
      std::cout << "  downstream=" << nsidepts[3] << std::endl;
      for (int i=0; i<4; i++)
        boundary_endpt_usedtag.at(i).resize( nsidepts[i], 0 );
    }

    // here we take those images and do some clustering.
    ann::ANNAlgo::cleanup();

    // ------------------------------------------------------------------------------------------//
    // FLASH TAGGER //
    std::cout << "[Run Flash Tagger]" << std::endl;

    // loop through flash producers, get event_opflash ptrs
    std::vector< larlite::event_opflash* > opflash_containers;
    for ( auto &flashproducer : flashtagger_pset.get< std::vector<std::string> >( "OpflashProducers" ) ) {
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
      opflash_containers.push_back( opdata );
    }

    std::vector< larlitecv::BoundarySpacePoint > trackendpts_anode;
    std::vector< larlitecv::BoundarySpacePoint > trackendpts_cathode;
    std::vector< larlitecv::BoundarySpacePoint > trackendpts_imgends;
    anode_flash_tagger.flashMatchTrackEnds( opflash_containers, imgs, badchimgs, trackendpts_anode  );
    cathode_flash_tagger.flashMatchTrackEnds( opflash_containers, imgs, badchimgs, trackendpts_cathode );
    imgends_flash_tagger.findImageTrackEnds( imgs, badchimgs, trackendpts_imgends );
    
    std::cout << "[[ Flash Tagger End Points ]]" << std::endl;
    std::cout << "    anode: " << trackendpts_anode.size() << std::endl;
    std::cout << "    cathode: " << trackendpts_cathode.size() << std::endl;
    std::cout << "    imgends: " << trackendpts_imgends.size() << std::endl;
    boundary_endpt_usedtag.at(4).resize( trackendpts_anode.size(), 0 );
    boundary_endpt_usedtag.at(5).resize( trackendpts_cathode.size(), 0 );
    boundary_endpt_usedtag.at(6).resize( trackendpts_imgends.size(), 0 );

//     if ( use_pos_flash_matching ) {
//       std::cout << "skipping ahead" << std::endl;
//       continue;
//     }

    // ------------------------------------------------------------------------------------------//
    // MAKE TRACKS USING COLLECTED END POINTS

    /// pixels of different types are unrolled
    std::vector< const larlitecv::BoundarySpacePoint* > all_endpoints;
    
    // gather endpoints from space points
    for (int isp=0; isp<(int)space_points.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(space_points.at( isp ));
      all_endpoints.push_back( pts );
    }

    // gather boundary points: anode/cathode/imageends
    for (int isp=0; isp<(int)trackendpts_anode.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(trackendpts_anode.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)trackendpts_cathode.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(trackendpts_cathode.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)trackendpts_imgends.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(trackendpts_imgends.at(isp));
      all_endpoints.push_back( pts );
    }

    // filter end points
    std::vector<int> endpoint_passes( all_endpoints.size(), 1 );
    endptfilter.removeBoundaryAndFlashDuplicates( all_endpoints, imgs, gapchimgs_v, endpoint_passes );
    //endptfilter.removeSameBoundaryDuplicates( all_endpoints, imgs, gapchimgs_v, endpoint_passes );
    //endptfilter.removeDiffBoundaryDuplicates( all_endpoints, imgs, gapchimgs_v, endpoint_passes );    

    // remove the filtered end points
    std::vector< const larlitecv::BoundarySpacePoint* > filtered_endpoints;
    for ( size_t idx=0; idx<endpoint_passes.size(); idx++ ) {
      if ( endpoint_passes.at(idx)==1 ) {
        filtered_endpoints.push_back( all_endpoints.at(idx) );
      }
    }

    // ------------------------------------------------------------------------------------------//
    // Check that the number of end points is not excessive
    bool runtracker = true;
    if ( (int)filtered_endpoints.size()>max_endpts_to_process ) {
      std::cout << "The number of total end points is excessive: " << filtered_endpoints.size() << ". The event is probably too busy." << std::endl;
      runtracker = false;
      continue;
    }

    // ------------------------------------------------------------------------------------------//
    // Form 3D tracks on each plane from boundary points
    std::vector< larcv::Image2D > tagged_v;
    std::vector< larlitecv::BMTrackCluster3D > tracks3d;
    std::vector< int > used_filtered_endpoints( filtered_endpoints.size(), 0 );
    if (runtracker) {
      sidetagger.makeTrackClusters3D( imgs, gapchimgs_v, filtered_endpoints, tracks3d, tagged_v, used_filtered_endpoints );
    }

    // ------------------------------------------------------------------------------------------//
    // Filter 2D tracks and form 3D tracks
    
    std::vector<int> goodlist( tracks3d.size(), 1);
    //sidetagger.process2Dtracks( trackclusters, imgs, gapchimgs_v, tracks3d, goodlist );
    
    std::cout << "[NUMBER OF POST-PROCESSED 3D TRAJECTORIES: " << tracks3d.size() << "]" << std::endl;

    // ------------------------------------------------------------------------------------------//


    //std::vector< larcv::Image2D > markedimgs;
    //sidetagger.markImageWithTrackClusters( imgs, gapchimgs_v, tracks3d, goodlist, markedimgs );

    // ------------------------------------------------------------------------------------------//
    // SAVE OUTPUT //

    // from mod channel imgs
    larcv::EventImage2D* ev_mod_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
    ev_mod_imgs->Emplace( std::move( imgs ) );

    // from empty channel imgs
    larcv::EventImage2D* ev_badch_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "badchs" );
    ev_badch_imgs->Emplace( std::move( badchimgs ) );
    
    // from gap channel imgs
    larcv::EventImage2D* ev_gapch_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "gapchs" );
    ev_gapch_imgs->Emplace( std::move( gapchimgs_v ) );
    
    // side tagger -- real space hits
    if ( sidetagger_cfg.save_endpt_images ) {
      larcv::EventImage2D* realspace_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "realspacehits" );
      realspace_imgs->Emplace( std::move(realspacehits) );

      // boundary pixels
      larcv::EventImage2D* boundarypixels_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "boundarypixels" );
      boundarypixels_imgs->Emplace( std::move(boundarypixels) );
    }

    // Save All End-Points
    larcv::EventPixel2D* realspace_endpts[larlitecv::kNumEndTypes];
    realspace_endpts[larlitecv::kTop]        = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("TopSpacePts") );
    realspace_endpts[larlitecv::kBottom]     = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("BottomSpacePts") );
    realspace_endpts[larlitecv::kUpstream]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("UpstreamSpacePts") );
    realspace_endpts[larlitecv::kDownstream] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, sidetagger_pset.get<std::string>("DownstreamSpacePts") );
    realspace_endpts[larlitecv::kAnode]      = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("AnodeEndpointProducer") );
    realspace_endpts[larlitecv::kCathode]    = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("CathodeEndpointProducer") );    
    realspace_endpts[larlitecv::kImageEnd]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, flashtagger_pset.get<std::string>("ImageEndpointProducer") );
    for (int i=0; i<(int)filtered_endpoints.size(); i++) {
      const larlitecv::BoundarySpacePoint& sp_v = *(filtered_endpoints.at(i));
      int sptype = (int)sp_v.at(0).type;
      if ( sptype<0 ) continue; // should
      std::cout << "space point #" << i << " type=" << sptype << std::endl;
      for (int p=0; p<3; p++) {
        const larlitecv::BoundaryEndPt& sp = sp_v.at(p);
        larcv::Pixel2D pixel( sp.col, sp.row );
        pixel.Intensity( sptype ); // using intensity to label pixel
        if ( used_filtered_endpoints.at(i)==1 )
          pixel.Width( 1 ); // using width to mark if used
        else
          pixel.Width( 0 ); // using width to mark if unused
        realspace_endpts[sptype]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }

    // Save Unused End-Points: for downstream steps
    larcv::EventPixel2D* unused_endpts[larlitecv::kNumEndTypes];
    unused_endpts[larlitecv::kTop]        = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+sidetagger_pset.get<std::string>("TopSpacePts") );
    unused_endpts[larlitecv::kBottom]     = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+sidetagger_pset.get<std::string>("BottomSpacePts") );
    unused_endpts[larlitecv::kUpstream]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+sidetagger_pset.get<std::string>("UpstreamSpacePts") );
    unused_endpts[larlitecv::kDownstream] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+sidetagger_pset.get<std::string>("DownstreamSpacePts") );
    unused_endpts[larlitecv::kAnode]      = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+flashtagger_pset.get<std::string>("AnodeEndpointProducer") );
    unused_endpts[larlitecv::kCathode]    = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+flashtagger_pset.get<std::string>("CathodeEndpointProducer") );    
    unused_endpts[larlitecv::kImageEnd]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_"+flashtagger_pset.get<std::string>("ImageEndpointProducer") );
    for (int i=0; i<(int)filtered_endpoints.size(); i++) {
      const larlitecv::BoundarySpacePoint& sp_v = *(filtered_endpoints.at(i));
      int sptype = (int)sp_v.at(0).type;
      if ( used_filtered_endpoints.at(i)==1 )
        continue;
      for (int p=0; p<3; p++) {
        const larlitecv::BoundaryEndPt& sp = sp_v.at(p);
        larcv::Pixel2D pixel( sp.col, sp.row );
        pixel.Intensity( sptype ); // using intensity to label pixel
        pixel.Width( 1 ); // using width to mark if used
        unused_endpts[sptype]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
      }
    }

    // save 2D track objects filtered by good 3d tracks
    larcv::EventPixel2D* ev_tracks2d = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "thrumu2d" );
    for (int i3d=0; i3d<(int)tracks3d.size(); i3d++) {
      larlitecv::BMTrackCluster3D& track3d = tracks3d.at(i3d);
      if ( goodlist.at( i3d )==0 ) continue;
      std::vector< larlitecv::BMTrackCluster2D >& trackcluster2d = track3d.plane_paths;
      for (int p=0; p<3; p++) {
        larlitecv::BMTrackCluster2D& track = trackcluster2d.at(p);
        larcv::Pixel2DCluster cluster;
        std::swap( cluster, track.pixelpath ); // ok to do this, because we will not use these pixels again
        std::cout << " plane=" << p << " track. length=" << cluster.size() << std::endl;
        ev_tracks2d->Emplace( (larcv::PlaneID_t)p, std::move(cluster) );
      }
    }

    // save 3D track object
    larlite::event_track* ev_tracks = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "thrumu3d" );
    
    // convert BMTrackCluster3D to larlite::track
    int id = 0;
    for ( int itrack=0; itrack<(int)tracks3d.size(); itrack++ ) {
      const larlitecv::BMTrackCluster3D& track3d = tracks3d.at(itrack);
      if ( goodlist.at( id )==0 ) continue;
      larlite::track lltrack;
      lltrack.set_track_id( id );
      int istep = 0;
      for ( auto const& point3d : track3d.path3d ) {
        TVector3 vec( point3d[0], point3d[1], point3d[2] );
        lltrack.add_vertex( vec );
        if ( istep+1<(int)track3d.path3d.size() ) {
          TVector3 dir( track3d.path3d.at(istep+1)[0]-point3d[0], track3d.path3d.at(istep+1)[1]-point3d[1], track3d.path3d.at(istep+1)[2]-point3d[2] );
          lltrack.add_direction( dir );
        }
        else {
          TVector3 dir(0,0,0);
          lltrack.add_direction( dir );
        }
      }
      std::cout <<  "storing 3D track with " << lltrack.NumberTrajectoryPoints() << " trajectory points" << std::endl;
      ev_tracks->emplace_back( std::move(lltrack) );
    }

    // Marked images
    larcv::EventImage2D* event_markedimgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
    event_markedimgs->Emplace( std::move(tagged_v) );


    // go to tree
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
