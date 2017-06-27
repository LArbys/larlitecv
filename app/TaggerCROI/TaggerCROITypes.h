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
    };

  };

  class ThruMuPayload : public TaggerCROIVPayload {
  public:

  ThruMuPayload() : TaggerCROIVPayload("ThruMu") { clear(); };
    virtual ~ThruMuPayload() {};

    std::vector< larcv::Image2D >     boundarypixel_image_v;
    std::vector< larcv::Image2D >     realspacehit_image_v;
    std::vector< larcv::Image2D >     tagged_v;

    std::vector< BoundarySpacePoint > side_spacepoint_v;
    std::vector< BoundarySpacePoint > anode_spacepoint_v;
    std::vector< BoundarySpacePoint > cathode_spacepoint_v;
    std::vector< BoundarySpacePoint > imgends_spacepoint_v;

    std::vector< BoundarySpacePoint > side_filtered_v;
    std::vector< BoundarySpacePoint > anode_filtered_v;
    std::vector< BoundarySpacePoint > cathode_filtered_v;
    std::vector< BoundarySpacePoint > imgends_filtered_v;
    
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
    std::vector< int > flashdata_selected_v;

    std::vector< int > flashdata_passes_containment_v;
    std::vector< int > flashdata_passes_cosmicflash_ratio_v;
    std::vector< int > flashdata_passes_flashmatch_v;
    std::vector< int > flashdata_passes_totpe_v;

    std::vector<float> containment_dwall_v;
    std::vector<float> min_chi2_v; ///< min-chi2 score for each flash
    std::vector<float> totpe_peratio_v; //< fraction difference of total pe between hypothesis and in-time
    std::vector<float> cosmicflash_ratio_dchi_v; //< delta-Chi-squared between in-time flash and matched flash to track    

    std::vector< larcv::ROI > croi_v;
    std::vector< larcv::Image2D > combined_v;

    std::vector< larlite::opflash > track_opflash_v;

    virtual void saveSpace() {
    };

    void clean() {
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
    }
  };
  
}

#endif
