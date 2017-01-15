#ifndef __FLASHROISELECTION__
#define __FLASHROISELECTION__

/**
	// Class to Filter Out untagged ROIs using flash information

	*/

#include <vector>

// larcv
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "PMTWeights/PMTWireWeights.h"

// larlite
#include "DataFormat/opflash.h"

namespace larlitecv {
	
	class FlashROIMatchingConfig {
	public:
			FlashROIMatchingConfig();
			virtual ~FlashROIMatchingConfig() {};

			void setDefaults();

			std::vector<int> beam_tick_range;
			float us_per_tick;
			float pmtflash_thresh;
	};

	class FlashROIMatching {
	public:
		FlashROIMatching( const FlashROIMatchingConfig& config);

		virtual ~FlashROIMatching() {};

		FlashROIMatchingConfig m_config;
		larcv::pmtweights::PMTWireWeights m_pmtweights;


		// primary routine
		std::vector< larcv::ROI > SelectFlashConsistentROIs( const std::vector<larcv::ROI>& rois, const std::vector<larlite::event_opflash*>& opflashes_v,
			const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, 
			const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v );

		// supporting routines
    std::vector<larlite::opflash> SelectInTimeFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );
    void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );

	};


}

#endif