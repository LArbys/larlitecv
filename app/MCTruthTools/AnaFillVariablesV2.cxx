#include "AnaFillVariablesV2.h"

#include "TTree.h"

#include "LArUtil/LArProperties.h"

namespace larlitecv {

  void AnaFillVariablesV2::bindEventTree( TTree* tree ) {
    if ( tree==NULL )
      throw std::runtime_error( "[AnaFillVariables::bindEventTree] Error: need to provide non-NULL pointer" );
    m_ev_tree = tree;
    m_ev_tree->Branch( "croipixelarea",  ev_croipixelarea,   "croipixelarea[4]/I" );
    m_ev_tree->Branch( "maxpe",          &ev_maxpe,          "maxpe/F" );
    m_ev_tree->Branch( "fwhm",           &ev_fwhm,           "fwhm/F" );
    m_ev_tree->Branch( "meanz",          &ev_meanz,          "meanz/F" );
    m_ev_tree->Branch( "highestnufrac",  &ev_highestnufrac,  "ev_highestnufrac/F" );
    //m_ev_tree->Branch( "highestvtxfrac", &ev_highestvtxfrac, "ev_highestvtxfrac/F" ); // to do
  }

  void AnaFillVariablesV2::bindTrackTree( TTree* tree ) {
    if ( tree==NULL )
      throw std::runtime_error( "[AnaFillVariables::bindTrackTree] Error: need to provide non-NULL pointer" );
    m_track_tree = tree;
    m_track_tree->Branch( "track_highestnufrac", &track_highestnufrac, "track_highestnufrac/I" );
    m_track_tree->Branch( "track_maxpe_hypo",    &track_maxpe_hypo,    "track_maxpe_hypo/F" );
    m_track_tree->Branch( "track_zfracdiff",     &track_zfracdiff,     "track_zfracdiff/F" );
    m_track_tree->Branch( "track_chi2",          &track_chi2,          "track_chi2/F" );
    m_track_tree->Branch( "track_dwall",         &track_dwall,         "track_dwall/F" );
    m_track_tree->Branch( "track_nufrac",        &track_nufrac,        "track_nufrac/F" );
    m_track_tree->Branch( "track_flashdtick",    &track_flashdtick,    "track_flashdtick/F" );
    m_track_tree->Branch( "track_mctrackid",     &track_mctrackid,     "track_mctrackid/I" );
    m_track_tree->Branch( "track_matchable",     &track_matchable,     "track_matchable/I" );
    m_track_tree->Branch( "track_matched",       &track_matched,       "track_matched/I" );
    m_track_tree->Branch( "track_flashmcid",     &track_flashmcid,     "track_flashmcid/I" );
  }

  void AnaFillVariablesV2::clearEventVars() {
    for (int p=0; p<4; p++)
      ev_croipixelarea[p] = -1;
    ev_maxpe = -1;
    ev_fwhm  = -1;
    ev_meanz = -1;
    ev_highestnufrac = -1;
  }

  void AnaFillVariablesV2::clearTrackVars() {
    track_maxpe_hypo = -1;
    track_zfracdiff  = -1;
    track_chi2       = -1;
    track_dwall      = -1;
    track_nufrac     = -1;
    track_flashdtick = -1;
    track_mctrackid  = -1;
    track_flashmcid  = -1;
    track_matchable  = -1;
    track_matched    = -1;
  }

  void AnaFillVariablesV2::fillEventInfo( const larcv::EventImage2D& ev_img, const larcv::EventROI& ev_roi,
					  larlite::event_user* ev_user_info, const larcv::EventImage2D* ev_segment,
					  larcv::EventPixel2D* ev_allpixels_v,
					  const larlite::event_track* ev_track,
					  const CrossingPointAnaData_t& xingptdata ) {
    clearEventVars();
    clearTrackVars();
    
    fillCROIPixArea( ev_img, ev_roi );
    fillIntimeFlashInfo( ev_user_info );
    fillTrackInfo( ev_user_info, ev_track, xingptdata, ev_segment, ev_allpixels_v );
  }
  
  void AnaFillVariablesV2::fillCROIPixArea( const larcv::EventImage2D& ev_img, const larcv::EventROI& ev_roi ) {
    const std::vector<larcv::Image2D>& img_v = ev_img.Image2DArray();
    const std::vector<larcv::ROI>&     roi_v = ev_roi.ROIArray();
    
    if ( img_v.size()!=3 )
      throw std::runtime_error( "[AnaFillVariablesV2::fillCROIPixArea] Error: the assumed number of planes (3) does not match number of images" );
    
    // we want to count the pixels covered
    for (int p=0; p<4; p++)
      ev_croipixelarea[p] = 0;
    
    // we do something stupid: we make blank images and mark pixels.
    // This we should use intersection of boxes -- but we try to get away with it for now
    // we want to prevent reallocating memory, so we create a cache of images which we will
    //   set to zero at each call.
    if ( _area_counter_v.size()!=img_v.size() ) {
      _area_counter_v.clear();
      for ( auto const& img : img_v ) {
	larcv::Image2D blank( img.meta() );
	_area_counter_v.emplace_back( std::move(blank) );
      }
    }

    for ( auto& img : _area_counter_v )
      img.paint(0);
    
    
    for ( auto const& roi : roi_v ) {
      for ( auto& img : _area_counter_v ) {
      
	const larcv::ImageMeta& bbox = roi.BB( img.meta().plane() );
	// tag the row and col range
	size_t minr = bbox.row(bbox.max_y());
	size_t minc = bbox.col(bbox.min_x());
	
	for ( size_t r=minr; r<minr+bbox.rows(); r++)
	  for (size_t c=minc; c<minc+bbox.cols(); c++)
	    img.set_pixel(r,c,1);
	
      }
    }
    
    // now count
    for ( auto& img : _area_counter_v ) {
      const larcv::ImageMeta& meta = img.meta();
      for ( size_t r=0; r<meta.rows(); r++ ) {
	for ( size_t c=0; c<meta.cols(); c++) {
	  ev_croipixelarea[meta.plane()] += img.pixel(r,c);
	  ev_croipixelarea[3]            += img.pixel(r,c);
	}
      }
    }
    
    return;
  }
  
  void AnaFillVariablesV2::fillIntimeFlashInfo( larlite::event_user* ev_user_info ) {
    if ( ev_user_info==NULL ) {
      ev_meanz = -1.0;
      ev_fwhm  = -1.0;
      ev_maxpe = -1.0;
      return;
    }

    // debug
    std::cout << "number of user info objects: " << ev_user_info->size() << std::endl;
    larlite::user_info& info = ev_user_info->front();
    if ( info.exist_darray("flash_meanz") )
      ev_meanz     = info.get_darray( "flash_meanz" )->front();
    if ( info.exist_darray("flash_zfwhm") )
      ev_fwhm      = info.get_darray( "flash_zfwhm" )->front();
    if ( info.exist_darray("flash_pemax") )
      ev_maxpe     = info.get_darray( "flash_pemax" )->front();
  }
  
  void AnaFillVariablesV2::fillTrackInfo( larlite::event_user* ev_user_info ,
					  const larlite::event_track* ev_track,
					  const CrossingPointAnaData_t& xingptdata,
					  const larcv::EventImage2D* ev_segment,
					  larcv::EventPixel2D* ev_allpixels_v ) {

    std::cout << "[AnaFillVariablesV2::fillTrackInfo]";
    std::cout << " number of reco tracks: " << ev_allpixels_v->Pixel2DClusterArray(0).size() << std::endl;
      
    // fill reco variables
    larlite::user_info& info = ev_user_info->front();

    if ( !info.exist_darray("cosmicflash_chi2")
    	 || !info.exist_darray("containment_dwall")
    	 || !info.exist_darray("track_zdiff_frac") ) {
      std::cout << "Missing Double Array" << std::endl;
      return;
    }
    
    std::vector<double>* pcosmicflash_chi2_v = info.get_darray( "cosmicflash_chi2" );
    std::vector<double>* pcontainment_v      = info.get_darray( "containment_dwall" );
    std::vector<double>* ptrack_zdiff_frac_v = info.get_darray( "track_zdiff_frac" );

    int idx_highest_nufrac = -1;
    float highest_nufrac = -1;
    std::vector<double> nufrac_v( pcosmicflash_chi2_v->size(), -1.0 );    
    if ( ev_segment!=NULL && ev_allpixels_v!=NULL) {
      for (size_t itrack=0; itrack<pcosmicflash_chi2_v->size(); itrack++) {
	float nufrac = fillTrackNuFraction( ev_segment, itrack, ev_allpixels_v );
	if ( nufrac > highest_nufrac ) {
	  highest_nufrac     = nufrac;
	  idx_highest_nufrac = itrack;
	}
	nufrac_v[itrack] = nufrac;
      }
    }
    ev_highestnufrac = highest_nufrac;
    
    for (size_t itrack=0; itrack<pcosmicflash_chi2_v->size(); itrack++) {
      track_zfracdiff = ptrack_zdiff_frac_v->at(itrack);
      track_chi2      = pcosmicflash_chi2_v->at(itrack);
      track_dwall     = pcontainment_v->at(itrack);
      track_nufrac    = nufrac_v[itrack];
      if ( (int)itrack==idx_highest_nufrac )
	track_highestnufrac = 1;
      else
	track_highestnufrac = 0;

      // flash matching checks
      isTrackFlashMatched( itrack, ev_track, xingptdata, ev_user_info );
      
      m_track_tree->Fill();
    }
    
  }
  
  float AnaFillVariablesV2::fillTrackNuFraction( const larcv::EventImage2D* ev_segment,
						 const int itrack, larcv::EventPixel2D* ev_allpixels_v ) {

    larcv::EventPixel2D& evpixels = *ev_allpixels_v;
    const std::vector<larcv::Image2D>& seg_v = ev_segment->Image2DArray();
    float num_nu_pixels = 0;
    float num_tagged_nu = 0;
    for ( auto const& seg : seg_v ) {
      const larcv::ImageMeta& meta = seg.meta();
      const larcv::Pixel2DCluster& pixel_v = evpixels.Pixel2DClusterArray(meta.plane()).at(itrack);
      for ( auto const& pix : pixel_v ) {
	num_nu_pixels += 1.0;
	if ( seg.pixel( pix.Y(), pix.X() )>0 )
	  num_tagged_nu += 1.0;
      }
    }
    float nufrac = 0.;
    if ( num_nu_pixels>0 )
      nufrac = num_tagged_nu/num_nu_pixels;
    std::cout << "[itrack " << itrack << "] track_nufrac=" << nufrac << " ("  << num_tagged_nu << "/" << num_nu_pixels << ")" << std::endl;
    return nufrac;
  }

  
  void AnaFillVariablesV2::isTrackFlashMatched( const int itrack,
						const larlite::event_track* ev_track,
						const CrossingPointAnaData_t& xingptdata,
						larlite::event_user* ev_user_info ) {

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    
    // get the index of the cosmic flash
    larlite::user_info& info = ev_user_info->front();
    std::vector<int>* cosmicflash_index_v     = info.get_iarray( "cosmicflash_index" );
    std::vector<int>* cosmicflash_mctrackid_v = info.get_iarray( "cosmicflash_trackmctrackid" );
    std::vector<int>* cosmicflash_flashmcid_v = info.get_iarray( "cosmicflash_flashmctrackid" );

      
    int cosmicflash_idx = cosmicflash_index_v->at(itrack);
    track_mctrackid = cosmicflash_mctrackid_v->at(itrack);
    track_flashmcid = cosmicflash_flashmcid_v->at(itrack);
    if ( track_mctrackid==-1 ) {
      track_matchable = 0;
      track_matched   = 0;
      return;
    }
    
    // is it matchable. we check to see if there is any flash with the mctrackid of the track
    track_matchable = 0;
    for ( auto const& flashinfo : xingptdata.flashanainfo_v ) {
      if ( flashinfo.mctrackid==track_mctrackid ) {
	track_matchable = 1;
	break;
      }
    }
    if ( track_matchable==0 ) {
      track_matched   =  0;
      return;
    }

    // ok, did it match?
    if ( track_mctrackid==track_flashmcid )
      track_matched = 1;
    else
      track_matched = 0;
    
  }
  
}
