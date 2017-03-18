#include "PayloadWriteMethods.h"

#include <vector>

// larlite
#include "DataFormat/opflash.h"

// larcv
#include "DataFormat/DataFormatTypes.h"
#include "DataFormat/EventImage2D.h"

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


	void WriteThruMuPayload( const ThruMuPayload& data, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

    // side tagger -- real space hits
    if ( config.sidetagger_cfg.save_endpt_images ) {
    	if ( config.thrumu_write_cfg.get<bool>("WriteRealSpaceHitsImage") ) {
	      larcv::EventImage2D* realspace_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "realspacehits" );
  	    for ( auto const& img : data.realspacehit_image_v ) realspace_imgs->Append( img );
  	  }

      // boundary pixels
  	  if ( config.thrumu_write_cfg.get<bool>("WriteBoundaryPixelsImage") ) {
        larcv::EventImage2D* boundarypixels_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "boundarypixels" );
        for ( auto const& img : data.boundarypixel_image_v ) boundarypixels_imgs->Append( img );
      }
    }

    // Save All End-Points
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
      spacepoint_lists.push_back( &(data.used_spacepoint_v) );
      spacepoint_lists.push_back( &(data.unused_spacepoint_v) );
      for ( auto const& plist : spacepoint_lists ) {
      	for ( auto const& sp_v : *plist ) {
          int sptype = (int)sp_v.front().type;
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
        int sptype = (int)sp_v.at(0).type;
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
        const std::vector< larlitecv::BMTrackCluster2D >& trackcluster2d = track3d.plane_paths;
        for (int p=0; p<3; p++) {
          const larlitecv::BMTrackCluster2D& track = trackcluster2d.at(p);
          ev_tracks2d->Append( (larcv::PlaneID_t)p, track.pixelpath );
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
}