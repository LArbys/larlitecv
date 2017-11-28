#ifndef __TAGGER_CROI_TYPES_H__
#define __TAGGER_CROI_TYPES_H__

/* --------------------------------
//  ThruMuPayload
//  Simply the data output by the ThruMu tagger algos
//
// ------------------------------- */

#include "TaggerCROIVPayload.h"

#include <vector>

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/trigger.h"
#include "DataFormat/mctrack.h"

// larcv
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Image2D.h"

// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "ContainedROI/TaggerFlashMatchTypes.h"
#include "UntaggedClustering/ClusterGroupMatchingTypes.h"
#include "T3DMerge/T3DCluster.h"

namespace larlitecv {

  class InputPayload : public TaggerCROIVPayload {
  public:
  InputPayload() : TaggerCROIVPayload("Input") {};
    virtual ~InputPayload() {};

    std::vector<larcv::Image2D>            img_v;       //< input image
    std::vector<larcv::Image2D>            badch_v;     //< image marking bad channels
    std::vector<larcv::Image2D>            gapch_v;     //< image marking gap channels
    std::vector< larlite::event_opflash* > opflashes_v; //< container of opflashes to consider
    larlite::trigger*                      p_ev_trigger; //< event trigger data
    int run;
    int subrun;
    int event;
    int entry;

    virtual void saveSpace() {};
    void clear() {
      img_v.clear();
      badch_v.clear();
      gapch_v.clear();
      opflashes_v.clear(); // we are assuming, this doesn't own the objects (dangerous)
      p_ev_trigger = NULL;
    };

    // mc truth information
    const larlite::event_mctrack* p_ev_mctrack; //< pointer to mc track object if used
  

  };

  class ThruMuPayload : public TaggerCROIVPayload {
  public:

  ThruMuPayload() : TaggerCROIVPayload("ThruMu") { clear(); };
    virtual ~ThruMuPayload() {};

    std::vector< larcv::Image2D >     boundarypixel_image_v;
    std::vector< larcv::Image2D >     realspacehit_image_v;
    std::vector< larcv::Image2D >     tagged_v;

    // raw spacepoints from side and flash end point taggers
    std::vector< BoundarySpacePoint > side_spacepoint_v;
    std::vector< BoundarySpacePoint > anode_spacepoint_v;
    std::vector< BoundarySpacePoint > cathode_spacepoint_v;
    std::vector< BoundarySpacePoint > imgends_spacepoint_v;
    
    // filtered spacepoints
    std::vector< BoundarySpacePoint > side_filtered_v;
    std::vector< BoundarySpacePoint > anode_filtered_v;
    std::vector< BoundarySpacePoint > cathode_filtered_v;
    std::vector< BoundarySpacePoint > imgends_filtered_v;

    // track objects
    std::vector< BMTrackCluster3D >      trackcluster3d_v;
    std::vector< larlite::track >        track_v;
    std::vector< std::vector<larcv::Pixel2DCluster> > pixelcluster_v;

    std::vector< BoundarySpacePoint > used_spacepoint_v;
    std::vector< BoundarySpacePoint > unused_spacepoint_v;

    void saveSpace() {
      // clears out data unneeded downstream
      boundarypixel_image_v.clear();
      realspacehit_image_v.clear();
      side_spacepoint_v.clear();
      anode_spacepoint_v.clear();
      cathode_spacepoint_v.clear();
      imgends_spacepoint_v.clear();
      used_spacepoint_v.clear();
    };

    void clear() {
      boundarypixel_image_v.clear();
      realspacehit_image_v.clear();
      tagged_v.clear();

      side_spacepoint_v.clear();
      anode_spacepoint_v.clear();
      cathode_spacepoint_v.clear();
      imgends_spacepoint_v.clear();

      side_filtered_v.clear();
      anode_filtered_v.clear();
      cathode_filtered_v.clear();
      imgends_filtered_v.clear();
      
      trackcluster3d_v.clear();
      track_v.clear();
      pixelcluster_v.clear();
      used_spacepoint_v.clear();
      unused_spacepoint_v.clear();
    }

  };

  class StopMuPayload : public TaggerCROIVPayload {
  public:
  StopMuPayload() : TaggerCROIVPayload("StopMu") {};
    virtual ~StopMuPayload() {};

    larcv::EventPixel2D                                 stopmu_pixel_endpt_v;     // bad design. fix this garbage.
    std::vector< std::vector< const larcv::Pixel2D* > > stopmu_candidate_endpt_v; // refers to previous line. terrible!

    std::vector< larlitecv::BMTrackCluster3D >        stopmu_trackcluster_v;
    std::vector< larlite::track >                     track_v;
    std::vector< std::vector<larcv::Pixel2DCluster> > pixelcluster_v;

    std::vector< larcv::Image2D > stopmu_v;

    virtual void saveSpace() {
      stopmu_candidate_endpt_v.clear();
      stopmu_pixel_endpt_v.clear();
    };
    void clear() {
      stopmu_pixel_endpt_v.clear();
      stopmu_candidate_endpt_v.clear();
      stopmu_trackcluster_v.clear();
      track_v.clear();
      pixelcluster_v.clear();
      stopmu_v.clear();
    }
  };

  class CROIPayload : public TaggerCROIVPayload {
  public:
  CROIPayload() : TaggerCROIVPayload("CROI") {};
    virtual ~CROIPayload() {};

    // Save the information for the Clustered Tracks (for the notebook).

    // Before Reclustering
    std::vector<larlitecv::TaggerFlashMatchData> contained_tracks_v;

    // After Reclustering.
    std::vector< T3DCluster > reclustered_contained_v;

    std::vector< larlitecv::T3DCluster > stopthru_reclustered_v;
    std::vector< std::vector<larcv::Pixel2DCluster> > stopthru_reclustered_pixels_v;

    std::vector< larlitecv::PlaneClusterGroups > plane_groups_v;
    std::vector< larlitecv::ChargeVolume > vols_v;
    std::vector< larcv::Image2D > tagged_v;
    std::vector< larcv::Image2D > subimg_v;
    std::vector< larlitecv::TaggerFlashMatchData > flashdata_v;

    // Result Variables
    std::vector< int > flashdata_selected_v;

    std::vector< int > flashdata_passes_containment_v;         //v1,v2
    std::vector< int > flashdata_passes_cosmicflash_ratio_v;   //v1,v2
    std::vector< int > flashdata_passes_flashmatch_v;          //v1
    std::vector< int > v2_flashdata_passes_flashpos_v;         //v2    
    std::vector< int > flashdata_passes_totpe_v;               //v1

    std::vector<double> containment_dwall_v;      //< v1,v2: smallest distance of track to wall
    std::vector<double> min_chi2_v;               //< v1: min-chi2 score for each flash
    std::vector<double> totpe_peratio_v;          //< v1: fraction difference of total pe between hypothesis and in-time
    std::vector<double> cosmicflash_ratio_dchi_v; //< v1:dChi^2 between in-time flash and matched flash to track. v2: best chi^2 to out-of-time flashes from flash matching
    std::vector<int>    cosmicflash_bestindex_v;  //< v2: index of out-of-time cosmic flash that best matches track
    std::vector<int>    cosmicflash_mctrackid_v;  //< v2: if we make thrumu tracks using mc truth (and no reclustering) mctrackid is passed to track objects (-1 if not available)
    std::vector<int>    cosmicflash_bestflash_mctrackid_v; //< v2: mctrackid of flash objects paired to track (-1 if could not pair)
    int num_matchable_flashes;                    //< v2: number of tracks we could match to a flash
    int num_matched_flashes;                      //< v2: number of tracks correctly matched to a flash

    std::vector<double> v2_intime_meanz_v;        //< v2: pe-weighted mean z position of flash
    std::vector<double> v2_intime_zfwhm_v;        //< v2: max z-range using pmts with > 0.5 of max pe values
    std::vector<double> v2_intime_pemax_v;        //< v2: max pe
    std::vector<double> v2_track_zdiff_frac_v;    //< v2: z distance difference of track end to flash meanz in units of zfhwm

    std::vector< larcv::ROI > croi_v;
    std::vector< larcv::Image2D > combined_v;

    std::vector< larlite::opflash > track_opflash_v;

    virtual void saveSpace() {
    };

    void clear() {
      contained_tracks_v.clear();
      reclustered_contained_v.clear();
      stopthru_reclustered_v.clear();
      stopthru_reclustered_pixels_v.clear();
      plane_groups_v.clear();
      vols_v.clear();
      tagged_v.clear();
      subimg_v.clear();
      flashdata_v.clear();
      flashdata_selected_v.clear();
      flashdata_passes_containment_v.clear();
      flashdata_passes_cosmicflash_ratio_v.clear();
      flashdata_passes_flashmatch_v.clear();
      flashdata_passes_totpe_v.clear();
      containment_dwall_v.clear();
      min_chi2_v.clear();
      totpe_peratio_v.clear();
      cosmicflash_ratio_dchi_v.clear();
      croi_v.clear();
      combined_v.clear();
      track_opflash_v.clear();
      v2_flashdata_passes_flashpos_v.clear();  
      v2_intime_meanz_v.clear();        //< v2: pe-weighted mean z position of flash
      v2_intime_zfwhm_v.clear();        //< v2: max z-range using pmts with > 0.5 of max pe values
      v2_intime_pemax_v.clear();        //< v2: max pe
      v2_track_zdiff_frac_v.clear();    //< v2: z distance difference of track end to flash meanz in units of zfhwm
      num_matchable_flashes = 0;
      num_matched_flashes   = 0;
      
    }
  };
  
}

#endif
