#ifndef __FLASHROISELECTION__
#define __FLASHROISELECTION__

/**
	// Class to Filter Out untagged ROIs using flash information

	*/

#include <vector>

#include "TTree.h"

// larcv
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Pixel2DCluster.h"
#include "PMTWeights/PMTWireWeights.h"

// larlite
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h" // SelectionTool

namespace larlitecv {
	
	class FlashROIMatchingConfig {
	public:
			FlashROIMatchingConfig();
			virtual ~FlashROIMatchingConfig() {};

			void setDefaults();

			std::vector<int> beam_tick_range;
			float us_per_tick;
			float pmtflash_thresh;
			bool store_calib_data;
	};

	class FlashROIMatching {

	public:
		FlashROIMatching( const FlashROIMatchingConfig& config);

		virtual ~FlashROIMatching() {};

		FlashROIMatchingConfig m_config;
		larcv::pmtweights::PMTWireWeights m_pmtweights;


		// primary routine
		std::vector< larcv::ROI > SelectFlashConsistentROIs( const std::vector<larlite::event_opflash*>& opflashes_v, const std::vector<larcv::Image2D>& img_v, 
			const std::vector< std::vector<larcv::Pixel2DCluster> >& untagged_clusters, const std::vector< larcv::ROI >& untagged_rois,
			const larcv::EventPixel2D& thrumu_clusters,
			const larcv::EventPixel2D& stopmu_clusters );

		// supporting routines
    std::vector<larlite::opflash> SelectInTimeFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );
    void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );

    // flash matching algorithm manager
    flashana::FlashMatchManager m_flash_matcher;

    // flash match data for calibration
    TTree* m_tree;
    int m_nuflag;
    int m_tagflag;
    float m_totalpe;
    float m_flashchi2;
    float m_flash_hypothesis[32];
    float m_measured[32];
    void writeCalibTree() { if ( m_config.store_calib_data) m_tree->Write(); }

	};


}

#endif
