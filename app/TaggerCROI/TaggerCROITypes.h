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
#include "ThruMu/BoundarySpacePoint.h"
#include "ThruMu/BMTrackCluster3D.h"
#include "ContainedROI/TaggerFlashMatchTypes.h"
#include "UntaggedClustering/ClusterGroupMatchingTypes.h"

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

	};

  class ThruMuPayload : public TaggerCROIVPayload {
  public:

  	ThruMuPayload() : TaggerCROIVPayload("ThruMu") {};
  	virtual ~ThruMuPayload() {};

  	std::vector< larcv::Image2D >     boundarypixel_image_v;
  	std::vector< larcv::Image2D >     realspacehit_image_v;
  	std::vector< larcv::Image2D >     tagged_v;

  	std::vector< BoundarySpacePoint > side_spacepoint_v;
  	std::vector< BoundarySpacePoint > anode_spacepoint_v;
  	std::vector< BoundarySpacePoint > cathode_spacepoint_v;
  	std::vector< BoundarySpacePoint > imgends_spacepoint_v;

  	std::vector< BMTrackCluster3D >   trackcluster3d_v;

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

  };

  class StopMuPayload : public TaggerCROIVPayload {
  public:
  	StopMuPayload() : TaggerCROIVPayload("StopMu") {};
  	virtual ~StopMuPayload() {};

  	larcv::EventPixel2D                                 stopmu_pixel_endpt_v;     // bad design. fix this garbage.
  	std::vector< std::vector< const larcv::Pixel2D* > > stopmu_candidate_endpt_v; // refers to previous line. terrible!
  	std::vector< larlitecv::BMTrackCluster3D >          stopmu_trackcluster_v;
  	std::vector< larcv::Image2D >                       stopmu_v;

  	virtual void saveSpace() {
	  stopmu_candidate_endpt_v.clear();
	  stopmu_pixel_endpt_v.clear();
	};
  };

  class CROIPayload : public TaggerCROIVPayload {
  public:
  	CROIPayload() : TaggerCROIVPayload("CROI") {};
  	virtual ~CROIPayload() {};

  	std::vector< larlitecv::PlaneClusterGroups > plane_groups_v;
  	std::vector< larlitecv::ChargeVolume > vols_v;
    std::vector< larcv::Image2D > tagged_v;
    std::vector< larcv::Image2D > subimg_v;
    std::vector< larlitecv::TaggerFlashMatchData > flashdata_v;
    std::vector< int > flashdata_selected_v;

    std::vector< larcv::ROI > croi_v;
    std::vector< larcv::Image2D > combined_v;

    std::vector< larlite::opflash > track_opflash_v;

  	virtual void saveSpace() {
  	};
  };

}

#endif
