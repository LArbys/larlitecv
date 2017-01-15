#include "FlashROIMatching.h"

// larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "FhiclLite/FhiclLiteUtilFunc.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"
#include "OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "OpT0Finder/Base/FlashFilterFactory.h"
#include "OpT0Finder/Base/TPCFilterFactory.h"
#include "OpT0Finder/Base/FlashMatchFactory.h"
#include "OpT0Finder/Base/FlashHypothesisFactory.h"
#include "OpT0Finder/Base/FlashProhibitFactory.h"
#include "OpT0Finder/Base/CustomAlgoFactory.h"

namespace larlitecv {
	
	FlashROIMatchingConfig::FlashROIMatchingConfig(fcllite::PSet ps_) : ps(ps_) {
		setDefaults();
	}

	void FlashROIMatchingConfig::setDefaults() {
		beam_tick_range.resize(2);
		beam_tick_range[0] = 150;
		beam_tick_range[1] = 350;
		us_per_tick = 0.015625;
		pmtflash_thresh = 5.0;
		store_calib_data = true;
		gain_correction.resize(32,1.0);
		pixel_threshold = 10.0;
		flash_front_boundary = 3200.0;
		flash_back_boundary  = 3200.0 + 4600.0;
	}

	FlashROIMatchingConfig MakeFlashROIMatchingConfigFromFile( std::string fname ) {
		fcllite::PSet ps = fcllite::CreatePSetFromFile(fname);
		fcllite::PSet main = ps.get<fcllite::PSet>("ContainedROI");
		FlashROIMatchingConfig cfg(main.get<fcllite::PSet>("FlashROIMatchingConfig"));
		//flashana::PhotonLibHypothesisFactory plfact;
		//flashana::QLLMatchFactory qllfact;
		return cfg;
	}

	FlashROIMatching::FlashROIMatching( const FlashROIMatchingConfig& config)
		: m_config(config), m_pmtweights("geoinfo.root") { 

			m_flash_matcher.Configure( m_config.ps );


			if ( m_config.store_calib_data ) {
				m_tree = new TTree("flashroi","Flash-ROI Matching Tree");
				m_tree->Branch( "nuflag",     &m_nuflag,    "nuflag/I" );
				m_tree->Branch( "tagflag",    &m_tagflag,   "tagflag/I");
				m_tree->Branch( "totalpe",    &m_totalpe,   "totalpe/F");
				m_tree->Branch( "flashchi2",  &m_flashchi2, "flashchi2/F");
				m_tree->Branch( "hypothesis", m_flash_hypothesis, "hypothesis[32]/F");
				m_tree->Branch( "measured",   m_measured,   "measured[32]/F");
			}
			else {
				m_tree = NULL;
			}
	}



	std::vector< larcv::ROI > FlashROIMatching::SelectFlashConsistentROIs( const std::vector<larlite::event_opflash*>& opflashes_v, 
		const std::vector<larcv::Image2D>& img_v, 
		const std::vector< std::vector<larcv::Pixel2DCluster> >& untagged_clusters,  const std::vector< larcv::ROI >& untagged_rois,
		larcv::EventPixel2D& thrumu_clusters,
		larcv::EventPixel2D& stopmu_clusters ) {

		// the vector we will fill
		std::vector<larcv::ROI> flash_matched_rois;

		// get all flashes in time with the beam
		std::vector<larlite::opflash> beam_flashes = SelectInTimeFlashes( opflashes_v );

		// reset flash match manager
		m_flash_matcher.Reset();
		std::cout << "FlashROIMatching ready." << std::endl;
		m_flash_matcher.PrintConfig();

		// get geometry info
		const ::larutil::Geometry* geo = ::larutil::Geometry::GetME();

		// for each flash we want to extract z-position and range, that is the first crude filter
		std::vector<float> wire_means;
		std::vector<std::vector<float>> wire_ranges;
		int flash_id = 0;
		for ( auto const& flash : beam_flashes ) {
			// provide mean and basic width
			float wire_mean;
			std::vector<float> wire_range;
			GetFlashCenterAndRange( flash, wire_mean, wire_range );
			wire_means.push_back( wire_mean );
			wire_ranges.push_back( wire_range );
			std::cout << "flash position: wire_mean=" << wire_mean << " wire_range=[" << (int)wire_range[0] << "," << (int)wire_range[1] << "]" << std::endl;

			// also load flash info into flash manager
			flashana::Flash_t f;
			f.x = f.x_err = 0;
			f.y = flash.YCenter();
			f.z = flash.ZCenter();
			f.y_err = flash.YWidth();
			f.z_err = flash.ZWidth();
			f.pe_v.resize(geo->NOpDets());
			f.pe_err_v.resize(geo->NOpDets());
			for (unsigned int i = 0; i < f.pe_v.size(); i++) {
	  		unsigned int opdet = geo->OpDetFromOpChannel(i);
	  		f.pe_v[opdet] = flash.PE(i) / m_config.gain_correction[i];
	  		f.pe_err_v[opdet] = sqrt(flash.PE(i) / m_config.gain_correction[i]);
			}
			f.time = flash.Time();
			f.idx = flash_id;
			++flash_id;
			m_flash_matcher.Emplace(std::move(f));
		}

		// filter all clusters consistent with this flash: untagged, thrumu, stopmu
		for ( int iflash=0; iflash<(int)beam_flashes.size(); iflash++ ) {
			const larlite::opflash& flash = beam_flashes.at(iflash);
			float wire_mean = wire_means.at(iflash);
			const std::vector<float>& wire_range = wire_ranges.at(iflash);

			// we check compatibility of all clusters (thrumu/stopmu/untagged)
			// note: untagged only for now to see how we do.  but we will have to fix thrumu/stopmu that have flash-tagged ends
			int roi_idx = 0;
			for ( auto const& untagged_roi : untagged_rois ) {
				const larcv::ImageMeta& yplane_bb = untagged_roi.BB().at(2);
				if ( (wire_range[0]<=yplane_bb.min_x() && yplane_bb.min_x()<=wire_range[1])
					|| (wire_range[0]<=yplane_bb.max_x() && yplane_bb.max_x()<=wire_range[1]) ) {
					flash_matched_rois.push_back( untagged_roi );
					std::cout << "flash matched roi: wire-range=[" << wire_range[0] << "-" << wire_range[1] << "] " 
						<< " y-plane roi=[" << yplane_bb.min_x() << "-" << yplane_bb.max_x() << "] "
						<< std::endl;
					// we make a qcluster object with this cluster
					flashana::QCluster_t qcluster = MakeTPCFlashObject( untagged_clusters.at(roi_idx).at(2), img_v.at(2) );
					m_flash_matcher.Emplace(std::move(qcluster));
				}
				roi_idx++;
			}// end of loop over untagged clusters/rois

			// thrumu clusters
			int thrumu_idx = 0;
			for ( const auto& pixel_cluster : thrumu_clusters.Pixel2DClusterArray((larcv::PlaneID_t)2) ) {
				std::cout << "thrumu cluster #" << thrumu_idx << " numpixels=" << pixel_cluster.size() << std::endl;
				if ( ClusterWithinRange(wire_range,pixel_cluster,img_v.at(2))) {
					std::cout << "thrumu cluster #" << thrumu_idx << " is compatible with flash." << std::endl;
					flashana::QCluster_t qcluster = MakeTPCFlashObject( pixel_cluster, img_v.at(2) );
					m_flash_matcher.Emplace(std::move(qcluster));
				}
				thrumu_idx++;
			}

			// stopmu clusters
			int stopmu_idx = 0;
			for ( const auto& pixel_cluster : stopmu_clusters.Pixel2DClusterArray((larcv::PlaneID_t)2) ) {
				std::cout << "stopmu cluster #" << stopmu_idx << " numpixels=" << pixel_cluster.size() << std::endl;
				if ( ClusterWithinRange(wire_range,pixel_cluster,img_v.at(2))) {
					std::cout << "stopmu cluster #" << stopmu_idx << " is compatible with flash." << std::endl;
					flashana::QCluster_t qcluster = MakeTPCFlashObject( pixel_cluster, img_v.at(2) );
					m_flash_matcher.Emplace(std::move(qcluster));
				}
				stopmu_idx++;
			}
		}

		std::cout << "Loaded " << m_flash_matcher.FlashArray().size() << " flashes "
			<< " and " << m_flash_matcher.QClusterArray().size() << " clusters" << std::endl;

		// run flash matching code
		m_results = m_flash_matcher.Match();
		std::cout << "Flash matching results" << std::endl;
		for ( auto &result : m_results ) {
			std::cout << "tpc_id=" << result.tpc_id << " flash_id="  << result.flash_id 
				<< " score=" << result.score
				<< std::endl;
		}

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

	bool FlashROIMatching::ClusterWithinRange( const std::vector<float>& wire_range, const larcv::Pixel2DCluster& pixels, const larcv::Image2D& img ) {
		const larcv::ImageMeta& meta = img.meta();
		bool within_range = false;
		for ( auto& pixel : pixels ) {
			if ( pixel.Y()<0 || pixel.Y()>=img.meta().rows()) continue;
			if ( pixel.X()<0 || pixel.X()>=img.meta().cols()) continue;
			if ( img.pixel(pixel.Y(),pixel.X())<m_config.pixel_threshold ) continue;
			float tick = meta.pos_y( pixel.Y() );
			float wire = meta.pos_x( pixel.X() );
			if ( wire_range[0]<=wire && wire <=wire_range[1] )
				within_range =  true;
			if ( tick<m_config.flash_front_boundary || tick>m_config.flash_back_boundary )
				return false;
		}
		return within_range;
	}

	flashana::QCluster_t FlashROIMatching::MakeTPCFlashObject( const larcv::Pixel2DCluster& pixels, const larcv::Image2D& img ) {
		flashana::QCluster_t qcluster;
		float cm_per_tick = 0.5*::larutil::LArProperties::GetME()->DriftVelocity();  // [usec/tick]*[cm/usec]
		for ( auto &pixel : pixels ) {
			flashana::QPoint_t pt;
			if ( pixel.Y()<img.meta().min_y() || pixel.Y()>=img.meta().max_y()) continue;
			if ( pixel.X()<0 || pixel.X()>=img.meta().max_x()) continue;
			if ( img.pixel(pixel.Y(),pixel.X())<m_config.pixel_threshold ) continue;
			int tick = img.meta().pos_y( pixel.Y() );
			int wire = img.meta().pos_x( pixel.X() );
		  pt.y = 0;
	  	pt.x = tick * cm_per_tick;
	  	pt.z = wire * 0.3;
	  	pt.q = img.pixel(pixel.Y(),pixel.X()) / 30. * 8000. * 0.3;
	  	qcluster.emplace_back(std::move(pt));
	  }
	  return qcluster;
	}

}