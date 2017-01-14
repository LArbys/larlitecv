#include "FlashROIMatching.h"

namespace larlitecv {
	
	FlashROIMatchingConfig::FlashROIMatchingConfig() {
		setDefaults();
	}

	void FlashROIMatchingConfig::setDefaults() {
		beam_tick_range.resize(2);
		beam_tick_range[0] = 150;
		beam_tick_range[1] = 350;
		us_per_tick = 0.015625;
	}


	std::vector< larcv::ROI > FlashROIMatching::SelectFlashConsistentROIs( const std::vector<larcv::ROI>& rois, const std::vector<larlite::event_opflash*>& opflashes_v,
			const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, 
			const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v ) {

		// the vector we will fill
		std::vector<larcv::ROI> flash_matched_rois;

		// get all flashes in time with the beam
		std::vector<larlite::opflash> beam_flashes = SelectInTimeFlashes( opflashes_v );

		// for each flash we want to extract z-position and range, that is the first crude filter

		// filter all clusters consistent with this flash: untagged, thrumu, stopmu
		// 

		// 

		return flash_matched_rois;
	}

	std::vector<larlite::opflash> FlashROIMatching::SelectInTimeFlashes( const std::vector<larlite::event_opflash*>& opflashes_v ) {
		std::vector<larlite::opflash> beam_flashes;

		for ( auto const& ptr_ev_flash : opflashes_v ) {
			for ( auto const& opflash : *ptr_ev_flash ) {
				int tick = opflash.Time()/m_config.us_per_tick;
				if ( tick>=m_config.beam_tick_range[0] && tick <=m_config.beam_tick_range[1] ) {
					std::cout << "In-time flash found: " << opflash.Time() << "us from trigger. Tick=" << tick << std::endl;
					beam_flashes.push_back( opflash ); // our own copy
				}
			}
		}

		std::cout << "FlashROIMatching::SelectInTimeFlashes. Found " << beam_flashes.size() << " flashes." << std::endl;
		return beam_flashes;
	}

}