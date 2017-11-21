#include "PayloadWriteMethods.h"

#include <vector>

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/user_info.h"

// larcv
#include "DataFormat/DataFormatTypes.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"

// larlitecv
#include "T3DMerge/T3D2LarliteTrack.h"

namespace larlitecv {

  void WriteInputPayload( const InputPayload& data, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

    // input
    if ( config.input_write_cfg.get<bool>("SaveInputTPC") ) {
      larcv::EventImage2D* evout_img = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimg" );
      for ( auto const& img : data.img_v )
	evout_img->Append( img );
    }

    // opflash
    if ( config.input_write_cfg.get<bool>("SaveOpflash") ) {
      larlite::event_opflash* evout_opflash = (larlite::event_opflash*)dataco.get_larlite_data( larlite::data::kOpFlash, "opflash" );
      for ( auto const& opflash_v : data.opflashes_v ) {
	for ( auto const& opflash : *opflash_v )
	  evout_opflash->push_back( opflash );
      }
    }

    // badch
    if ( config.input_write_cfg.get<bool>("SaveBadChImage") ) {
      larcv::EventImage2D* evout_badch = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "badch" );
      for ( auto const& img : data.badch_v )
	evout_badch->Append( img );
    }

    // gapchs
    if ( config.input_write_cfg.get<bool>("SaveGapChImage") ) {
      larcv::EventImage2D* evout_gapchs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "gapchs" );
      for ( auto const& img : data.gapch_v )
	evout_gapchs->Append( img );
    }

  }


  void WriteThruMuPayload( const ThruMuPayload& data, const InputPayload& input, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

    // side tagger -- real space hits (mostly for debugging)
    if ( config.sidetagger_cfg.save_endpt_images ) {
      if ( config.thrumu_write_cfg.get<bool>("WriteRealSpaceHitsImage") ) {
        larcv::EventImage2D* realspace_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "realspacehits" );
        for ( auto const& img : data.realspacehit_image_v ) realspace_imgs->Append( img );
      }

      // boundary pixels (mostly for debugging)
      if ( config.thrumu_write_cfg.get<bool>("WriteBoundaryPixelsImage") ) {
        larcv::EventImage2D* boundarypixels_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "boundarypixels" );
        for ( auto const& img : data.boundarypixel_image_v ) boundarypixels_imgs->Append( img );
      }
    }

    // Save All pre-filter end points (mostly for debugging)
    if ( config.thrumu_write_cfg.get<bool>("WritePrefilterBoundaryPoints") ) {
      larcv::EventPixel2D* ev_prefilter_endpts = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "prefilterpts" );
      std::vector< const std::vector<BoundarySpacePoint>* > prefiltered_v;
      prefiltered_v.push_back( &(data.side_spacepoint_v) );
      prefiltered_v.push_back( &(data.anode_spacepoint_v) );
      prefiltered_v.push_back( &(data.cathode_spacepoint_v) );
      prefiltered_v.push_back( &(data.imgends_spacepoint_v) );
      for ( auto const& plist : prefiltered_v ) {
        for ( auto const& sp_v : *plist ) {
          int sptype = (int)sp_v.type();
          if ( sptype<0 ) continue; // bad endpoint type
          for (size_t p=0; p<sp_v.size(); p++) {
            const larlitecv::BoundaryEndPt& sp = sp_v.at(p);
            larcv::Pixel2D pixel( sp.col, sp.row );
            pixel.Intensity( sptype ); // using intensity to label pixel
            ev_prefilter_endpts->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
          }
        }
      }
    }

    // Save All (post-filter) End-Points
    if ( config.thrumu_write_cfg.get<bool>("WriteAllBoundaryPoints") ) {
      larcv::EventPixel2D* realspace_endpts[larlitecv::kNumEndTypes];
      realspace_endpts[larlitecv::kTop]        = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "topspacepts" );
      realspace_endpts[larlitecv::kBottom]     = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "botspacepts" );
      realspace_endpts[larlitecv::kUpstream]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "upspacepts" );
      realspace_endpts[larlitecv::kDownstream] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "downspacepts" );
      realspace_endpts[larlitecv::kAnode]      = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "anodepts" );
      realspace_endpts[larlitecv::kCathode]    = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "cathodepts" );
      realspace_endpts[larlitecv::kImageEnd]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "imgendpts" );
      std::vector< const std::vector<BoundarySpacePoint>* > spacepoint_lists;
      // spacepoint_lists.push_back( &(data.used_spacepoint_v) );
      spacepoint_lists.push_back( &(data.side_filtered_v) );
      spacepoint_lists.push_back( &(data.anode_filtered_v) );
      spacepoint_lists.push_back( &(data.cathode_filtered_v) );
      spacepoint_lists.push_back( &(data.imgends_filtered_v) );            

      // spacepoint_lists.push_back( &(data.unused_spacepoint_v) );
      for ( auto const& plist : spacepoint_lists ) {
      	for ( auto const& sp_v : *plist ) {
          int sptype = (int)sp_v.type();
          if ( sptype<0 ) continue; // should
          for (size_t p=0; p<sp_v.size(); p++) {
            const larlitecv::BoundaryEndPt& sp = sp_v.at(p);
            larcv::Pixel2D pixel( sp.col, sp.row );
            pixel.Intensity( sptype ); // using intensity to label pixel
            if ( plist==&(data.used_spacepoint_v) )
              pixel.Width( 1 ); // using width to mark if used
            else
              pixel.Width( 0 ); // using width to mark if unused
            realspace_endpts[sptype]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
          }
        }
      }
    }//end of if write all boundary space points

    // Save Unused End-Points: for downstream steps
    if ( config.thrumu_write_cfg.get<bool>("WriteUnusedBoundaryPoints") ) {
      larcv::EventPixel2D* unused_endpts[larlitecv::kNumEndTypes];
      unused_endpts[larlitecv::kTop]        = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_topspacepts" );
      unused_endpts[larlitecv::kBottom]     = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_botspacepts" );
      unused_endpts[larlitecv::kUpstream]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_upspacepts" );
      unused_endpts[larlitecv::kDownstream] = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_downspacepts" );
      unused_endpts[larlitecv::kAnode]      = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_anodepts" );
      unused_endpts[larlitecv::kCathode]    = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_cathodepts" );
      unused_endpts[larlitecv::kImageEnd]   = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "unused_imgendpts" );
      for ( auto const& sp_v : data.unused_spacepoint_v ) {
        int sptype = (int)sp_v.type();
        for (int p=0; p<3; p++) {
          const larlitecv::BoundaryEndPt& sp = sp_v.at(p);
          larcv::Pixel2D pixel( sp.col, sp.row );
          pixel.Intensity( sptype ); // using intensity to label pixel
          pixel.Width( 0 ); // using width to mark if used
          unused_endpts[sptype]->Emplace( (larcv::PlaneID_t)p, std::move(pixel) );
        }
      }
    }

    // save 2D track objects filtered by good 3d tracks
    if ( config.thrumu_write_cfg.get<bool>("WriteThruMuPixels") ) {
      larcv::EventPixel2D* ev_tracks2d = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "thrumupixels" );
      for (int i3d=0; i3d<(int)data.trackcluster3d_v.size(); i3d++) {
        const larlitecv::BMTrackCluster3D& track3d = data.trackcluster3d_v.at(i3d);
        std::vector< larcv::Pixel2DCluster > trackcluster2d = track3d.getTrackPixelsFromImages( input.img_v, input.badch_v,
          config.thrumu_tracker_cfg.pixel_threshold, config.thrumu_tracker_cfg.tag_neighborhood, 0.3 );
        for (int p=0; p<3; p++) {
          ev_tracks2d->Append( (larcv::PlaneID_t)p, trackcluster2d.at(p), input.img_v.at(p).meta() );
        }
      }
    }

    // save 3D track object
    if ( config.thrumu_write_cfg.get<bool>("WriteThruMuTracks") ) {
      larlite::event_track* ev_tracks = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "thrumu3d" );

      // convert BMTrackCluster3D to larlite::track
      int id = 0;
      for ( int itrack=0; itrack<(int)data.trackcluster3d_v.size(); itrack++ ) {
        const larlitecv::BMTrackCluster3D& track3d = data.trackcluster3d_v.at(itrack);
        larlite::track lltrack = track3d.makeTrack();
        lltrack.set_track_id( id );
        ev_tracks->emplace_back( std::move(lltrack) );
      }
    }

    // Marked images
    if ( config.thrumu_write_cfg.get<bool>("WriteThruMuTaggedImage") ) {
      larcv::EventImage2D* event_markedimgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
      for ( auto const& tagged : data.tagged_v ) event_markedimgs->Append( tagged );
    }

  }

  void WriteStopMuPayload( const StopMuPayload& data, const InputPayload& input, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

    // save 3D track object
    if ( config.stopmu_write_cfg.get<bool>("WriteStopMuTracks") ) {
      larlite::event_track* ev_tracks = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "stopmu3d" );

      // convert BMTrackCluster3D to larlite::track
      for ( int itrack=0; itrack<(int)data.stopmu_trackcluster_v.size(); itrack++ ) {
        const larlitecv::BMTrackCluster3D& track3d = data.stopmu_trackcluster_v.at(itrack);
        larlite::track lltrack = track3d.makeTrack();
        lltrack.set_track_id( itrack );
        ev_tracks->emplace_back( std::move(lltrack) );
      }
    }

    // output: stopmu-tagged pixels
    // use pixel2dclusters to fill out image
    if ( config.stopmu_write_cfg.get<bool>("WriteStopMuTaggedImage") ) {
      larcv::EventImage2D* stopmu_eventimgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stopmu" );
      for ( auto const& stopmu : data.stopmu_v ) {
        stopmu_eventimgs->Append( stopmu );
      }
    }

    // finally, store 2D pixels
    if ( config.stopmu_write_cfg.get<bool>("WriteStopMuPixels") ) {
      larcv::EventPixel2D* ev_stopmupixels = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "stopmupixels" );
      for ( size_t itrack=0; itrack<data.stopmu_trackcluster_v.size(); itrack++ ) {
        const larlitecv::BMTrackCluster3D& track3d = data.stopmu_trackcluster_v.at(itrack);
        std::vector< larcv::Pixel2DCluster > trackpixs_v = track3d.getTrackPixelsFromImages( input.img_v, input.badch_v,
											     config.thrumu_tracker_cfg.pixel_threshold,
											     config.thrumu_tracker_cfg.tag_neighborhood, 0.3 );
        for (size_t p=0; p<trackpixs_v.size(); p++) {
          ev_stopmupixels->Append( (larcv::PlaneID_t)p, trackpixs_v.at(p), data.stopmu_v.at(p).meta() );
        }
      }
    }
  }

  void WriteCROIPayload( const CROIPayload& data, const InputPayload& inputdata, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

    // CRITICAL OUTPUT: NO FLAGS

    // ROIs, StopMu Tracks and clusters, ThruMu Tracks and Clusters, Truth ROI, Bad channels
    larcv::EventROI* out_ev_roi = (larcv::EventROI*)dataco.get_larcv_data( larcv::kProductROI, "croi" );
    out_ev_roi->Set( data.croi_v );

    // store CROI pixels
    larcv::EventPixel2D* out_croi_pixels = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "croipixels" );
    for ( size_t itrack=0; itrack<data.flashdata_v.size(); itrack++) {
      auto const& flashdata = data.flashdata_v.at(itrack);
      if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 )
        continue;
      if ( data.flashdata_selected_v[itrack]==0 )
        continue;

      for (size_t p=0; p<flashdata.m_pixels.size(); p++) {
        out_croi_pixels->Append( (larcv::PlaneID_t)p, flashdata.m_pixels.at(p), inputdata.img_v.at(p).meta() );
      }
    }

    // untagged track
    if ( config.croi_write_cfg.get<bool>("WriteUntaggedTracks") ) { 
      larlite::event_track* evout_tracks_untagged = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "untagged3d" );

      for ( size_t itrack=0; itrack<data.flashdata_v.size(); itrack++) {
        auto const& flashdata = data.flashdata_v.at(itrack);
        if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 )
          continue;
        if ( flashdata.m_type==larlitecv::TaggerFlashMatchData::kUntagged ) {
          evout_tracks_untagged->push_back( flashdata.m_track3d );
        }
      }
    }

    // store untagged pixels
    if ( config.croi_write_cfg.get<bool>("WriteUntaggedPixels") ) {
      larcv::EventPixel2D* out_untagged_pixels = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "untaggedpixels" );
      for ( size_t itrack=0; itrack<data.flashdata_v.size(); itrack++) {
        auto const& flashdata = data.flashdata_v.at(itrack);
        if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 )
          continue;

        if ( flashdata.m_type==larlitecv::TaggerFlashMatchData::kUntagged ) {
          for (size_t p=0; p<flashdata.m_pixels.size(); p++) {
            out_untagged_pixels->Append( (larcv::PlaneID_t)p, flashdata.m_pixels.at(p), inputdata.img_v.at(p).meta() );
          }
        }
      }
    }

    // store reclustered pixelsx
    if ( config.croi_write_cfg.get<bool>("WriteReclusteredPixels")) {
      larcv::EventPixel2D* out_streclustered_pixels = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "streclusteredpixels" );
      for ( size_t itrack=0; itrack<data.stopthru_reclustered_pixels_v.size(); itrack++) {
	const std::vector< larcv::Pixel2DCluster >& pixel_v = data.stopthru_reclustered_pixels_v[itrack];
	for (size_t p=0; p<pixel_v.size(); p++) {
	  out_streclustered_pixels->Append( (larcv::PlaneID_t)p, pixel_v[p], inputdata.img_v[p].meta() );
	}
      }
    }

    // store reclustered tracks
    if ( config.croi_write_cfg.get<bool>("WriteReclusteredTracks")) {
      larlite::event_track* evout_tracks_streclustered = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "streclustered3d" );
      for ( size_t itrack=0; itrack<data.stopthru_reclustered_v.size(); itrack++) {
	larlite::track lltrack = larlitecv::T3D2LarliteTrack( data.stopthru_reclustered_v[itrack] );
	evout_tracks_streclustered->emplace_back( std::move(lltrack) );
      }

      // ---------------------------------------------------------------------------------------------------------
      larlite::event_track* evout_tracks_alltracks = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "all3d" );
      larlite::event_track* evout_tracks_mergedstopmu3d = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "mergedstopmu3d" );
      larlite::event_track* evout_tracks_mergedthrumu3d = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "mergedthrumu3d" );
      larlite::event_track* evout_tracks_mergeduntagged = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "mergeduntagged3d" );

      for ( size_t itrack=0; itrack<data.flashdata_v.size(); itrack++) {
        auto const& flashdata = data.flashdata_v.at(itrack);
        if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 )
          continue;
	evout_tracks_alltracks->push_back( flashdata.m_track3d );

	if ( flashdata.m_type==larlitecv::TaggerFlashMatchData::kStopMu )
	  evout_tracks_mergedstopmu3d->push_back( flashdata.m_track3d );
	else if ( flashdata.m_type==larlitecv::TaggerFlashMatchData::kThruMu )
	  evout_tracks_mergedthrumu3d->push_back( flashdata.m_track3d );
	else if ( flashdata.m_type==larlitecv::TaggerFlashMatchData::kUntagged )
	  evout_tracks_mergeduntagged->push_back( flashdata.m_track3d );
      }

      // ---------------------------------------------------------------------------------------------------------

      larcv::EventPixel2D* out_pixels_alltracks      = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "allpixels" );
      larcv::EventPixel2D* out_pixels_mergedstopmu3d = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "mergedstopmupixels" );
      larcv::EventPixel2D* out_pixels_mergedthrumu3d = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "mergedthrumupixels" );
      larcv::EventPixel2D* out_pixels_mergeduntagged = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "mergeduntaggedpixels" );
      for ( size_t itrack=0; itrack<data.flashdata_v.size(); itrack++) {
	const larlitecv::TaggerFlashMatchData& flashtrack = data.flashdata_v[itrack];
	const std::vector< larcv::Pixel2DCluster >& pixel_v = flashtrack.m_pixels;
	// all
	for (size_t p=0; p<pixel_v.size(); p++) {
	  out_pixels_alltracks->Append( (larcv::PlaneID_t)p, pixel_v[p], inputdata.img_v[p].meta() );

	  if ( flashtrack.m_type==larlitecv::TaggerFlashMatchData::kStopMu )
	    out_pixels_mergedstopmu3d->Append( (larcv::PlaneID_t)p, pixel_v[p], inputdata.img_v[p].meta() );
	  else if ( flashtrack.m_type==larlitecv::TaggerFlashMatchData::kThruMu )
	    out_pixels_mergedthrumu3d->Append( (larcv::PlaneID_t)p, pixel_v[p], inputdata.img_v[p].meta() );
	  else if ( flashtrack.m_type==larlitecv::TaggerFlashMatchData::kUntagged )
	    out_pixels_mergeduntagged->Append( (larcv::PlaneID_t)p, pixel_v[p], inputdata.img_v[p].meta() );
	}
      }

    }

    // selected tracks
    if ( config.croi_write_cfg.get<bool>("WriteSelectedTracks")) {
      larlite::event_track* evout_tracks_selected = (larlite::event_track*)dataco.get_larlite_data( larlite::data::kTrack, "croi3d" );

      for ( size_t itrack=0; itrack<data.flashdata_v.size(); itrack++) {
        auto const& flashdata = data.flashdata_v.at(itrack);
        if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 )
          continue;
        if ( data.flashdata_selected_v[itrack] ) {
          evout_tracks_selected->push_back( flashdata.m_track3d );
        }
      }
    }

    // combined tagged image
    if ( config.croi_write_cfg.get<bool>("WriteCombinedTaggedImage") ) {
      larcv::EventImage2D* evout_combined_v = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "combinedtags");
      larcv::EventImage2D* evout_thrumutagged_v = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "thrumutags");
      larcv::EventImage2D* evout_stopmutagged_v = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stopmutags");
      larcv::EventImage2D* evout_containedtagged_v = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "containedtags");
      for ( auto const& combined : data.combined_v ) {
	// split the tags
	larcv::Image2D thrumu( combined.meta() );
	larcv::Image2D stopmu( combined.meta() );
	larcv::Image2D contained( combined.meta() );
	thrumu.paint(0);
	stopmu.paint(0);
	contained.paint(0);
	for (int r=0; r<(int)combined.meta().rows(); r++) {
	  for (int c=0; c<(int)combined.meta().cols(); c++) {
	    float val = combined.pixel(r,c);
	    if ( val>5 && val<15 )
	      thrumu.set_pixel(r,c,255);
	    else if ( val>15  && val<25 )
	      stopmu.set_pixel(r,c,255);
	    else if ( val>25 && val<35 )
	      contained.set_pixel(r,c,255);
	  }
	}

	evout_thrumutagged_v->Emplace( std::move(thrumu) );
	evout_stopmutagged_v->Emplace( std::move(stopmu) );
	evout_containedtagged_v->Emplace( std::move(contained) );	
        evout_combined_v->Append( combined );
      }
    }

    // Store opflashes for all tracks
    if (config.croi_write_cfg.get<bool>("WriteTrackOpFlashes") ) {
      larlite::event_opflash* track_opflash_out = (larlite::event_opflash*)dataco.get_larlite_data( larlite::data::kOpFlash, "ophypo");

      for ( auto const& opflash : data.track_opflash_v) {
        track_opflash_out->push_back(opflash);
      }
    }

    // Store cut results
    if ( config.croi_write_cfg.get<bool>("WriteCutResults") ) {
      larlite::event_user* ev_user_info = (larlite::event_user*)dataco.get_larlite_data( larlite::data::kUserInfo, "croicutresults" );
      larlite::user_info cutinfo;

      if ( config.croi_selection_cfg.use_version==1 ) {

	for ( auto const& result : data.flashdata_passes_containment_v )
	  cutinfo.append( "containment", result );
	for ( auto const& fdwall : data.containment_dwall_v )
	  cutinfo.append( "containment_dwall", fdwall );
	
	for ( auto const& result : data.flashdata_passes_flashmatch_v )
	  cutinfo.append( "flashmatch", result );
	for ( auto const& intimechi2 : data.min_chi2_v )
	  cutinfo.append( "flashmatch_intimechi2", intimechi2 );
	
	for ( auto const& result : data.flashdata_passes_cosmicflash_ratio_v )
	  cutinfo.append( "cosmicratio", result );
	for ( auto const& dchi2 : data.cosmicflash_ratio_dchi_v ) {
	  cutinfo.append( "cosmicratio_dchi2", dchi2 );
	}
	
	for ( auto const& result : data.flashdata_passes_totpe_v )
	  cutinfo.append( "totalpe", result );
	for ( auto const& totpe_ratio : data.totpe_peratio_v ) {
	  cutinfo.append( "totalpe_peratio",totpe_ratio );
	}
      }
      else if ( config.croi_selection_cfg.use_version==2 ) {

	// cut result: one value per track
	for ( auto const& result : data.flashdata_passes_containment_v )
	  cutinfo.append( "containment", result );
	for ( auto const& result : data.v2_flashdata_passes_flashpos_v )
	  cutinfo.append( "flashpos", result );
	for ( auto const& result : data.flashdata_passes_cosmicflash_ratio_v )
	  cutinfo.append( "cosmicflash", result );

	// cut variables
	for ( auto const& result : data.containment_dwall_v )
	  cutinfo.append( "containment_dwall", result );	
	for ( auto const& result : data.cosmicflash_ratio_dchi_v )
	  cutinfo.append( "cosmicflash_chi2", result );	
	for ( auto const& result : data.v2_intime_meanz_v )
	  cutinfo.append( "flash_meanz", result );	
	for ( auto const& result : data.v2_intime_zfwhm_v )
	  cutinfo.append( "flash_zfwhm", result );	
	for ( auto const& result : data.v2_intime_pemax_v )
	  cutinfo.append( "flash_pemax", result );	
	for ( auto const& result : data.v2_track_zdiff_frac_v )
	  cutinfo.append( "track_zdiff_frac", result );	
	
      }
      else {
	throw std::runtime_error("WriteCROIPayload: unsupported version specified");
      }

      ev_user_info->emplace_back( std::move(cutinfo) );
    }
  }
}
