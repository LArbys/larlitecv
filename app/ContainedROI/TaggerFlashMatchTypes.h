#ifndef __TAGGER_FLASH_MATCH_TYPES_H__
#define __TAGGER_FLASH_MATCH_TYPES_H__

#include <vector>

// LArCV
#include "DataFormat/Pixel2DCluster.h"

// LArLite
#include "DataFormat/track.h"

namespace larlitecv {

	class TaggerFlashMatchData {
	public:
		typedef enum { kThruMu=0, kStopMu, kUntagged } ClusterType_t;

		TaggerFlashMatchData( ClusterType_t type, const std::vector<larcv::Pixel2DCluster>& pixels, const larlite::track& track ) 
		 : m_type(type), m_pixels(pixels), m_track3d(track) {};
		virtual ~TaggerFlashMatchData() {};

	  ClusterType_t m_type;
		std::vector<larcv::Pixel2DCluster> m_pixels;
		larlite::track m_track3d;

	};
 
}

#endif