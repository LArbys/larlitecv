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
		pmtflash_thresh = 5.0;
	}

	FlashROIMatching::FlashROIMatching( const FlashROIMatchingConfig& config)
		: m_config(config), m_pmtweights("geoinfo.root") { 

	}



	std::vector< larcv::ROI > FlashROIMatching::SelectFlashConsistentROIs( const std::vector<larcv::ROI>& rois, const std::vector<larlite::event_opflash*>& opflashes_v,
			const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, 
			const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v ) {

		// the vector we will fill
		std::vector<larcv::ROI> flash_matched_rois;

		// get all flashes in time with the beam
		std::vector<larlite::opflash> beam_flashes = SelectInTimeFlashes( opflashes_v );

		// for each flash we want to extract z-position and range, that is the first crude filter
		std::vector<float> wire_means;
		std::vector<std::vector<float>> wire_ranges;
		for ( auto const& flash : beam_flashes ) {
			float wire_mean;
			std::vector<float> wire_range;
			GetFlashCenterAndRange( flash, wire_mean, wire_range );
			wire_means.push_back( wire_mean );
			wire_ranges.push_back( wire_range );
			std::cout << "flash position: wire_mean=" << wire_mean << " wire_range=[" << (int)wire_range[0] << "," << (int)wire_range[1] << "]" << std::endl;
		}

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

	void FlashROIMatching::GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange ) {

		zmean = 0.;
		zrange.resize(2,0.0);
    float qtot = 0.0;
    for (int ipmt=0; ipmt<32; ipmt++) {
    	zmean += flash.PE(ipmt)*m_pmtweights.pmtpos[ipmt][2];
    	qtot  += flash.PE( ipmt );
    }

    if (qtot>0)
    	zmean /= qtot;

    float mindist = 1000;
    float maxdist = 0;

    for (int ipmt=0; ipmt<32; ipmt++) {
      if (flash.PE(ipmt)>m_config.pmtflash_thresh) {
      	float dist = m_pmtweights.pmtpos[ipmt][2]-zmean;
      	if (dist<0 || mindist>dist)
      		mindist = dist;
      	else if ( dist>0 && maxdist<dist )
      		maxdist = dist;
      }
    }

    if (qtot<=0 )  {
    	// no information, so set wide range
    	zrange[0] = 0;
    	zrange[1] = 3455;
    }
    else {
			zrange[0] = (zmean+mindist-10.0)/0.3; // conversion from cm -> wire
			zrange[1] = (zmean+maxdist+10.0)/0.3; // conversion from cm -> wire
		}

		// bounds check
		if ( zrange[0]<0 ) zrange[0] = 0;
		if ( zrange[1]>=3455 ) zrange[1] = 3455;

	}

}