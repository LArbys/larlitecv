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
#include "DataFormat/Image2D.h"

// larlitecv
#include "ThruMu/BoundarySpacePoint.h"
#include "ThruMu/BMTrackCluster3D.h"

namespace larlitecv {

	class InputPayload : public TaggerCROIVPayload {
	public:
		InputPayload() : TaggerCROIVPayload("Input") {};
		virtual ~InputPayload() {};

		std::vector<larcv::Image2D>           img_v;       //< input image
		std::vector<larcv::Image2D>           badch_v;     //< image marking gap and bad channels
		std::vector< larlite::event_opflash* > opflashes_v; //< container of opflashes to consider

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
  	};

  };

}

#endif