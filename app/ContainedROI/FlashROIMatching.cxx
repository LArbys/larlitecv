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
	
	flashana::PhotonLibHypothesisFactory plfact;
	flashana::QLLMatchFactory qllfact;

	FlashROIMatchingConfig::FlashROIMatchingConfig(fcllite::PSet ps_) : ps(ps_) {
		setDefaults();
	}

	void FlashROIMatchingConfig::setDefaults() {
		beam_tick_range.resize(2);
		beam_tick_range[0] = 150;
		beam_tick_range[1] = 400;
		us_per_tick = 0.015625;
		pmtflash_thresh = 5.0;
		store_calib_data = true;
		gain_correction.resize(32,1.0);
		pixel_threshold = 10.0;
		flash_front_boundary = 3200.0;
		flash_back_boundary  = 3200.0 + 4600.0;
		min_qcluster_size = 5;
		verbosity = 0;
		qneighborhood = 5;
		trigger_tick = 3200;
	}

	FlashROIMatchingConfig MakeFlashROIMatchingConfigFromFile( std::string fname ) {
		fcllite::PSet ps = fcllite::CreatePSetFromFile(fname);
		fcllite::PSet main = ps.get<fcllite::PSet>("ContainedROI");
		fcllite::PSet flashmatchcfg = main.get<fcllite::PSet>("FlashROIMatchingConfig");
		FlashROIMatchingConfig cfg(flashmatchcfg);
		cfg.store_calib_data = flashmatchcfg.get<bool>("StoreCalibData");
		cfg.verbosity = flashmatchcfg.get<int>("Verbosity",0);
		std::vector<int> input_beam_range = flashmatchcfg.get<std::vector<int> >("BeamTickRange");
		if ( input_beam_range.size()!=2 )
			throw std::runtime_error("Invalid BeamTickRange size. Should be 2 (value for lower and upper bound)");
		cfg.beam_tick_range[0] = input_beam_range[0];
		cfg.beam_tick_range[1] = input_beam_range[1];
		//cfg.verbosity = 0;
		return cfg;
	}

	FlashROIMatching::FlashROIMatching( const FlashROIMatchingConfig& config)
		: m_config(config), m_pmtweights("geoinfo.root") { 

			ptr_segimg = NULL;
			SetVerbosity( config.verbosity );

			m_flash_matcher.Configure( m_config.ps );
			m_flash_matcher.GetAlgo((flashana::Algorithm_t)2)->set_verbosity( (flashana::msg::Level_t)0 );
			//m_flash_matcher.set_verbosity( (flashana::msg::Level_t)0 );

			if ( m_config.store_calib_data ) {
				m_file = new TFile("output_flashmatch_data.root", "recreate");
				m_tree = new TTree("flashroi","Flash-ROI Matching Tree");
				m_tree->Branch( "entry",      &m_entry,     "entry/I" );
				m_tree->Branch( "run",	      &m_run,       "run/I" );
				m_tree->Branch( "subrun",     &m_subrun,    "subrun/I" );
				m_tree->Branch( "event",     	&m_event,     "event/I" );
				m_tree->Branch( "nuflag",     &m_nuflag,    "nuflag/I" );
				m_tree->Branch( "tagflag",    &m_tagflag,   "tagflag/I");
				m_tree->Branch( "clustertype",&m_clustertype,"clustertype/I");
				m_tree->Branch( "passcuts",   &m_passcuts,  "passcuts/I");
				m_tree->Branch( "score",      &m_score,     "score/F");
				m_tree->Branch( "totalpe",    &m_totalpe,   "totalpe/F");
				m_tree->Branch( "predictpe",  &m_predictpe, "predictpe/F");
				m_tree->Branch( "rawtotpe",   &m_rawtotpe,  "rawtotpe/F");
				m_tree->Branch( "flashchi2",  &m_flashchi2, "flashchi2/F");
				m_tree->Branch( "fracnupixels", &m_fraction_nupixels, "fractionupixels/F" );
				m_tree->Branch( "flashrange", m_flash_range,"flash_range[2]/F" );
				m_tree->Branch( "hypothesis", m_flash_hypothesis, "hypothesis[32]/F");
				m_tree->Branch( "raw_hypothesis", m_raw_hypothesis, "raw_hypothesis[32]/F");
				m_tree->Branch( "measured",   m_measured,   "measured[32]/F");
			}
			else {
				m_tree = NULL;
				m_file = NULL;
			}
	}



	std::vector< larcv::ROI > FlashROIMatching::SelectFlashConsistentROIs( const std::vector<larlite::event_opflash*>& opflashes_v, 
		const std::vector<larcv::Image2D>& img_v, 
		const std::vector< std::vector<larcv::Pixel2DCluster> >& untagged_clusters,  const std::vector< larcv::ROI >& untagged_rois,
		larcv::EventPixel2D& thrumu_clusters,
		larcv::EventPixel2D& stopmu_clusters ) {

		if ( m_verbosity>0 ) {
			std::cout << "=========================================================================" << std::endl;
			std::cout << "======== FlashROIMatching: Selecting Candidate Clusters =================" << std::endl;
		}

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
		if ( m_verbosity>0 ) {
			std::cout << "=========================================================================" << std::endl;
			std::cout << "======== FlashROIMatching: Selecting Candidate Clusters =================" << std::endl;
			std::cout << "In-Time Flashes" << std::endl;
		}

		for ( auto const& flash : beam_flashes ) {
			// provide mean and basic width
			float wire_mean;
			std::vector<float> wire_range;
			GetFlashCenterAndRange( flash, wire_mean, wire_range );
			wire_means.push_back( wire_mean );
			wire_ranges.push_back( wire_range );
			if ( m_verbosity>0 )
				std::cout << "Flash #" << flash_id << ": wire_mean=" << wire_mean << " wire_range=[" << (int)wire_range[0] << "," << (int)wire_range[1] << "]" << std::endl;

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

		std::vector< FlashROIMatching::CandidateFlashMatchedROI_t > candidates = GenerateFlashTPCObjects( beam_flashes, wire_ranges, 
			untagged_clusters, untagged_rois, thrumu_clusters, stopmu_clusters, img_v );


		// count number of pixels around neutrino interaction
		int nsegpixels = 0;
		std::cout << "segment ptr=" << ptr_segimg << std::endl;
		if ( ptr_segimg!=NULL && m_config.store_calib_data ) {
			nsegpixels = NumOfSegPixels( img_v.at(2) );
			std::cout << "number of segment pixels=" << nsegpixels << std::endl;
		}

		// run flash matching code	
		m_results = GenerateFittedFlashHypotheses( candidates );
		std::cout << "Flash matching results: " << m_results.size() << std::endl;

		m_results_unfitted = GenerateUnfittedFlashHypotheses( candidates );

		// make the selection
		std::set< std::pair<int,int> > cluster_indices;
		SelectCandidateROIs( candidates, m_results_unfitted, cluster_indices );

		// the vector we will fill
		std::vector<larcv::ROI> flash_matched_rois;
		for ( auto& candidate_pair : cluster_indices ) {
			std::vector< larcv::Pixel2DCluster > plane_clusters;
			for ( size_t p=0; p<3; p++ ) {
				if ( candidate_pair.first==0 )
					plane_clusters.push_back(  untagged_clusters.at( candidate_pair.second ).at(p)  );
				else if ( candidate_pair.first==1 )
					plane_clusters.push_back(  thrumu_clusters.Pixel2DClusterArray(p).at( candidate_pair.second ) );
				else if ( candidate_pair.first==2 )
					plane_clusters.push_back(  stopmu_clusters.Pixel2DClusterArray(p).at( candidate_pair.second ) );
			}

			larcv::ROI roi = GenerateClusterROI( plane_clusters, img_v );
			flash_matched_rois.emplace_back( std::move(roi));
		}


		// Save output to ROOT file
		if( m_config.store_calib_data) {
			for ( auto &result : m_results ) {
				std::cout << "tpc_id=" << result.tpc_id << " flash_id="  << result.flash_id 
					<< " score=" << result.score
					<< std::endl;
				const CandidateFlashMatchedROI_t& candidate = candidates.at( result.tpc_id );
				const flashana::FlashMatch_t& unfitted_flashmatch = m_results_unfitted.at( result.tpc_id );

				const larcv::Pixel2DCluster& expanded_cluster = candidate.expanded_clusters.at(0);
				int noverlap = NumOfSegPixelOverlap( expanded_cluster, img_v.at(2).meta() );

				//// generate un-optimized flash
				//const flashana::QLLMatch* matchalgo = (flashana::QLLMatch*)m_flash_matcher.GetAlgo( flashana::kFlashMatch );
				//flashana::Flash_t raw_hypothesis = matchalgo->GetEstimate( m_flash_matcher.QClusterArray().at(result.tpc_id ) );

				// save information about the flash
				m_fraction_nupixels = float(noverlap)/float(nsegpixels);
				m_score = result.score;
				m_totalpe = 0;
				m_predictpe = 0;
				m_rawtotpe = 0;
				m_passcuts = 0;
				for (int ipmt=0; ipmt<32; ipmt++) {
					m_flash_hypothesis[ipmt] = result.hypothesis.at(ipmt);
					m_raw_hypothesis[ipmt] = unfitted_flashmatch.hypothesis.at(ipmt);
					m_measured[ipmt] = m_flash_matcher.FlashArray().at(result.flash_id).pe_v.at(ipmt);
					m_totalpe += m_measured[ipmt];
					m_predictpe += m_flash_hypothesis[ipmt];
					m_rawtotpe += m_raw_hypothesis[ipmt];
				}
				m_clustertype = candidate.source_type;
				m_flash_range[0] = wire_ranges.at( result.flash_id ).at(0);
				m_flash_range[1] = wire_ranges.at( result.flash_id ).at(1);
				m_flashchi2 = unfitted_flashmatch.score;

				std::pair<int,int> candidate_key( candidate.source_type, candidate.source_idx );
				auto it_passing_candidates = cluster_indices.find( candidate_key );
				if ( it_passing_candidates!=cluster_indices.end() )
					m_passcuts = 1;

				m_tree->Fill();
			}
		}

		return flash_matched_rois;
	}

	std::vector<larlite::opflash> FlashROIMatching::SelectInTimeFlashes( const std::vector<larlite::event_opflash*>& opflashes_v ) {
		std::vector<larlite::opflash> beam_flashes;
		for ( auto const& ptr_ev_flash : opflashes_v ) {
			for ( auto const& opflash : *ptr_ev_flash ) {
				int tick = opflash.Time()/m_config.us_per_tick;
				if ( tick>=m_config.beam_tick_range[0] && tick <=m_config.beam_tick_range[1] ) {
					if ( m_verbosity>0 )
						std::cout << "In-time flash found: " << opflash.Time() << "us from trigger. Tick=" << tick << std::endl;
					beam_flashes.push_back( opflash ); // our own copy
				}
				else {
					if ( m_verbosity>0)
						std::cout << "Rejected flash: " << opflash.Time() << "us from trigger. Tick=" << tick << std::endl;
				}
			}
		}

		if ( m_verbosity>0 )
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

	flashana::QCluster_t FlashROIMatching::MakeTPCFlashObject( const larcv::Pixel2DCluster& pixels, const larcv::Image2D& img, 
		larcv::Pixel2DCluster& expanded_cluster ) {
		// around each pixel we look around a neighborhood for pixels with charge.
		// we make a QPoint for each pixel. We take care not fo double count

		larcv::Image2D qtracker( img.meta() );
		qtracker.paint(0.0);

		flashana::QCluster_t qcluster;
		float cm_per_tick = 0.5*::larutil::LArProperties::GetME()->DriftVelocity();  // [usec/tick]*[cm/usec]
		for ( auto &pixel : pixels ) {
			flashana::QPoint_t pt;
			int row = pixel.Y();
			int col = pixel.X();
			for (int dr=-m_config.qneighborhood; dr<=m_config.qneighborhood; dr++) {
				int r = row+dr;
				if ( r<0 || r>=(int)img.meta().rows() ) continue;
				for ( int dc=-m_config.qneighborhood; dc<=m_config.qneighborhood; dc++) {
					int c = col+dc;
					if ( c<0 || c>=(int)img.meta().cols() ) continue;
					if ( img.pixel(r,c)<m_config.pixel_threshold) continue;

					if ( qtracker.pixel(r,c)>0 ) continue; // already counted

					int tick = img.meta().pos_y( r );
					int wire = img.meta().pos_x( c );
				  pt.y = 0;
	  			pt.x = (tick-m_config.trigger_tick) * cm_per_tick;
	  			if ( pt.x<0 ) continue; // not physical
	  			pt.z = wire * 0.3;
	  			pt.q = img.pixel(r,c) / 30. * 8000. * 0.3;
	  			qcluster.emplace_back(std::move(pt));
	  			qtracker.set_pixel( r,c, 1.0 );
	  			expanded_cluster.push_back( larcv::Pixel2D(c,r) );
				}//loop over c neighborhood
			}//loop over r neighborhood
	  }//loop over track pixels
	  return qcluster;
	}

	int FlashROIMatching::NumOfSegPixels( const larcv::Image2D& img ) {
		int nsegpixels = 0;
		if ( ptr_segimg==NULL ) return 0;

		const larcv::ImageMeta& meta = ptr_segimg->meta();
		for (size_t r=0; r<meta.rows(); r++) {
			for ( size_t c=0; c<meta.cols(); c++ ) {
				if ( ptr_segimg->pixel(r,c)>0 ) {
					float tick = meta.pos_y(r);
					float wire = meta.pos_x(c);
					int imgrow = img.meta().row( tick );
					int imgcol = img.meta().col( wire );
					if ( img.pixel(imgrow,imgcol)>m_config.pixel_threshold )
						nsegpixels++;
				}
			}
		}
		return nsegpixels;
	}

	int FlashROIMatching::NumOfSegPixelOverlap( const larcv::Pixel2DCluster& pixels, const larcv::ImageMeta& meta ) {
		int noverlap = 0;
		if ( ptr_segimg==NULL ) return 0;

		const larcv::ImageMeta& segmeta = ptr_segimg->meta();
		for ( auto const& pixel : pixels ) {
			float tick = meta.pos_y( pixel.Y() );
			float wire = meta.pos_x( pixel.X() );
			if ( tick<=segmeta.min_y() || tick>=segmeta.max_y() ) continue;
			if ( wire<=segmeta.min_x() || wire>=segmeta.max_x() ) continue;
			try {
				int segrow = segmeta.row( tick );
				int segcol = segmeta.col( wire );
				if ( ptr_segimg->pixel(segrow,segcol)>0 ) noverlap++;
			}
			catch (...) {} // move on
		}
		return noverlap;
	}

	std::vector<FlashROIMatching::CandidateFlashMatchedROI_t> FlashROIMatching::GenerateFlashTPCObjects( const std::vector<larlite::opflash>& beam_flashes, 
		const std::vector< std::vector<float> >& flash_wire_ranges,
		const std::vector< std::vector<larcv::Pixel2DCluster> >& untagged_clusters, const std::vector< larcv::ROI >& untagged_rois,
		larcv::EventPixel2D& thrumu_clusters,
		larcv::EventPixel2D& stopmu_clusters,
		const std::vector<larcv::Image2D>&  img_v ) {
		// This routine loops through all of the 'track' clusters. 
		// For all those that are compatible with the different flashes, we 
		//  1) generate a new Pixel2DCluster on each plane, this time making sure to collect as much charge as possible
		//  2) use this information to generate a QCluster_t object that will be used to generate various flash hypotheses
		//  3) 

		std::vector<FlashROIMatching::CandidateFlashMatchedROI_t> output; // container of the data we're trying to collect and organize

		// filter all clusters consistent with this flash: untagged, thrumu, stopmu
		int qcluster_idx = 0;

		// we check compatibility of all clusters (thrumu/stopmu/untagged)
		// note: untagged only for now to see how we do.  but we will have to fix thrumu/stopmu that have flash-tagged ends

		// UNTAGGED ROIs
		for ( int roi_idx=0; roi_idx<(int)untagged_rois.size(); roi_idx++ ) {
			const larcv::ROI& untagged_roi = untagged_rois.at(roi_idx);
			const std::vector<larcv::Pixel2DCluster>& untagged_cluster = untagged_clusters.at(roi_idx);
			const larcv::ImageMeta& yplane_bb = untagged_roi.BB().at(2);

			if ( m_verbosity>0 )
				std::cout << "Untagged cluster #" << roi_idx 
					<< ": wire-span [" << yplane_bb.min_x() << ", " << yplane_bb.max_x() << "]"
					<< "  tick-span=[" << yplane_bb.min_y() << ","  << yplane_bb.max_y() << "]";

			// check flash compatibility
			std::vector<int> compatible_flashes;
			for ( int iflash=0; iflash<(int)beam_flashes.size(); iflash++ ) {
				const std::vector<float>& wire_range = flash_wire_ranges.at(iflash);
				if ( (wire_range[0]<=yplane_bb.min_x() && yplane_bb.min_x()<=wire_range[1])
					|| (wire_range[0]<=yplane_bb.max_x() && yplane_bb.max_x()<=wire_range[1]) ) {

					compatible_flashes.push_back(iflash);
				}
			}

			if ( m_verbosity>0 )
				std::cout << " Num compatible Flashes=" << compatible_flashes.size() << std::endl;
			if ( compatible_flashes.size()==0 ) {
				continue; // no matching flashes
			}

			// compatible ROI
			CandidateFlashMatchedROI_t candidate;
			candidate.flash_idxs = compatible_flashes;

			// we make a qcluster object with this cluster
			larcv::Pixel2DCluster expanded_cluster;
			candidate.qcluster = MakeTPCFlashObject( untagged_cluster.at(2), img_v.at(2), expanded_cluster );
			if ( (int)candidate.qcluster.size()>m_config.min_qcluster_size ) {
				// we will keep this one
				candidate.qcluster_idx = qcluster_idx;
				candidate.source_type = 0; // untagged roi
				candidate.source_idx = roi_idx;
				candidate.expanded_clusters.emplace_back( std::move(expanded_cluster) );
				qcluster_idx++;				
				output.emplace_back( std::move(candidate) );
			}// end of loop over untagged clusters/rois
		}

		// thrumu clusters
		for ( int thrumu_idx=0; thrumu_idx<(int)thrumu_clusters.Pixel2DClusterArray( (larcv::PlaneID_t)2 ).size(); thrumu_idx++ ) {
			const larcv::Pixel2DCluster& pixel_cluster = thrumu_clusters.Pixel2DClusterArray((larcv::PlaneID_t)2).at(thrumu_idx);
			//std::cout << "thrumu cluster #" << thrumu_idx << " numpixels=" << pixel_cluster.size() << std::endl;
			if ( m_verbosity>0 ) {
				// this is a lot for some output...
				std::vector< larcv::Pixel2DCluster > temp;
				temp.push_back( pixel_cluster );
				larcv::ROI thrumu_roi = GenerateClusterROI( temp, img_v );
				std::cout << "Thrumu cluster #" << thrumu_idx 
					<< ": wire-span [" << thrumu_roi.BB(0).min_x() << ", " << thrumu_roi.BB(0).max_x() << "]"
					<< "  tick-span=[" << thrumu_roi.BB(0).min_y() << "," << thrumu_roi.BB(0).max_y() << "]";
			}


			// check flash compatibility
			std::vector<int> compatible_flashes;
			for ( int iflash=0; iflash<(int)beam_flashes.size(); iflash++ ) {
				const std::vector<float>& wire_range = flash_wire_ranges.at(iflash);
				if ( ClusterWithinRange(wire_range,pixel_cluster,img_v.at(2))) {
					compatible_flashes.push_back(iflash);
				}
			}

			if ( m_verbosity>0 )
				std::cout << " Num compatible Flashes=" << compatible_flashes.size() << std::endl;


			if ( compatible_flashes.size()==0 ) continue; // no matching flashes
			//std::cout << "thrumu cluster #" << thrumu_idx << " is compatible with flash." << std::endl;

			// compatible ROI
			CandidateFlashMatchedROI_t candidate;
			candidate.flash_idxs = compatible_flashes;

			// we make a qcluster object with this cluster
			larcv::Pixel2DCluster expanded_cluster;			
			candidate.qcluster = MakeTPCFlashObject( pixel_cluster, img_v.at(2), expanded_cluster );
			if ( (int)candidate.qcluster.size()>m_config.min_qcluster_size ) {
				// we will keep this one
				candidate.qcluster_idx = qcluster_idx;
				candidate.source_type = 1; // thrumu
				candidate.source_idx = thrumu_idx;
				candidate.expanded_clusters.emplace_back( std::move(expanded_cluster) );
				qcluster_idx++;
				output.emplace_back( std::move(candidate) );
			}
		}//end of thrumu loop

		// stopmu clusters
		for ( int stopmu_idx=0; stopmu_idx<(int)stopmu_clusters.Pixel2DClusterArray((larcv::PlaneID_t)2).size(); stopmu_idx++ ) {
			const larcv::Pixel2DCluster& pixel_cluster = stopmu_clusters.Pixel2DClusterArray((larcv::PlaneID_t)2).at(stopmu_idx);
			if ( pixel_cluster.size()<1 ) {
				if ( m_verbosity>0 )
					std::cout << "StopMu cluster #" << stopmu_idx << " too small." << std::endl;
				continue;
			}
			if ( m_verbosity>0 ) {
				// this is a lot for some output...
				std::vector< larcv::Pixel2DCluster > temp;
				temp.push_back( pixel_cluster );
				larcv::ROI stopmu_roi = GenerateClusterROI( temp, img_v );
				std::cout << "StopMu cluster #" << stopmu_idx 
					<< ": wire-span [" << stopmu_roi.BB(0).min_x() << ", " << stopmu_roi.BB(0).max_x() << "]"
					<< "  tick-span=[" << stopmu_roi.BB(0).min_y() << ","  << stopmu_roi.BB(0).max_y() << "]";
			}

			// check flash compatibility
			std::vector<int> compatible_flashes;
			for ( int iflash=0; iflash<(int)beam_flashes.size(); iflash++ ) {
				const std::vector<float>& wire_range = flash_wire_ranges.at(iflash);
				if ( ClusterWithinRange(wire_range,pixel_cluster,img_v.at(2))) {
					compatible_flashes.push_back(iflash);
				}
			}

			if ( m_verbosity>0 )
				std::cout << " Num compatible Flashes=" << compatible_flashes.size() << std::endl;

			if ( compatible_flashes.size()==0 ) continue; // no matching flashes
			//std::cout << "stopmu cluster #" << stopmu_idx << " is compatible with flash." << std::endl;

			// compatible ROI
			CandidateFlashMatchedROI_t candidate;
			candidate.flash_idxs = compatible_flashes;

			// we make a qcluster object with this cluster
			larcv::Pixel2DCluster expanded_cluster;			
			candidate.qcluster = MakeTPCFlashObject( pixel_cluster, img_v.at(2), expanded_cluster );
			if ( (int)candidate.qcluster.size()>m_config.min_qcluster_size ) {
				// we will keep this one
				candidate.qcluster_idx = qcluster_idx;
				candidate.source_type = 2; // stopmu
				candidate.source_idx = stopmu_idx;
				candidate.expanded_clusters.emplace_back( std::move(expanded_cluster) );
				qcluster_idx++;
				output.emplace_back( std::move(candidate) );
			}
		}//end of stopmu loop

		std::cout << "ID'd " << output.size() << " in-time flash compatible TPC 'track' clusters using loose z-position cut." << std::endl;

		return output;

	}

	
	std::vector<flashana::FlashMatch_t> FlashROIMatching::GenerateFittedFlashHypotheses( std::vector< CandidateFlashMatchedROI_t> & candidates ) {
		// Use FlashMatchManager (QLLMatch+PhotonLibraryHypothesis) to generate flash hypotheses
		for ( auto& candidate : candidates ) {
			m_flash_matcher.Add( candidate.qcluster );
		}
		return m_flash_matcher.Match();
	}

	
	std::vector< flashana::FlashMatch_t > FlashROIMatching::GenerateUnfittedFlashHypotheses( std::vector<CandidateFlashMatchedROI_t>& candidates ) {
		std::vector< flashana::FlashMatch_t > output;
		const flashana::QLLMatch* matchalgo = (flashana::QLLMatch*)m_flash_matcher.GetAlgo( flashana::kFlashMatch );		

		for ( auto& candidate : candidates ) {
			// first we loop over flashes
			std::vector<flashana::FlashMatch_t> unfitted_matches;
			int bestmatch = 0;
			double best_chi2 = 1e6;
			for ( int flash_idx=0; flash_idx<(int)candidate.flash_idxs.size(); flash_idx++ ) {
				const flashana::Flash_t& theflash = m_flash_matcher.FlashArray().at(flash_idx);
				flashana::Flash_t unfitted_hypothesis = matchalgo->GetEstimate( candidate.qcluster );

				double chi2 = 0.;
				int ndf = 0;
				for ( int ipmt=0; ipmt<(int)unfitted_hypothesis.pe_v.size(); ipmt++ ) {
					double diff = unfitted_hypothesis.pe_v.at(ipmt)-theflash.pe_v.at(ipmt);
					double err  = sqrt( unfitted_hypothesis.pe_v.at(ipmt)+theflash.pe_v.at(ipmt) );
					if ( err!=0.0 ) {
						chi2 += (diff*diff)/(err*err);
						ndf++;
					}
				}
				if ( ndf-1>0 )
					chi2 /= float(ndf-1);
				else
					chi2 = -1.0;

				if ( (flash_idx==0) || best_chi2>chi2 ) {
					bestmatch = flash_idx;
					best_chi2 = chi2;
				}

				flashana::FlashMatch_t unfit_match( candidate.qcluster_idx, flash_idx, chi2 );
				unfit_match.hypothesis = unfitted_hypothesis.pe_v;
				unfitted_matches.emplace_back( std::move(unfit_match) );
			}

			output.emplace_back( std::move(unfitted_matches.at(bestmatch)) );
		}

		return output;
	}

	void FlashROIMatching::SelectCandidateROIs( std::vector<CandidateFlashMatchedROI_t>& candidates, 
		std::vector< flashana::FlashMatch_t>& matches, std::set< std::pair<int,int> >& cluster_indices ) {
		// here we decide which ROIs to pass on
		// We pass on
		//  1) uncontained clusters if
		//       a) ratio of raw/measured total pe is > 0.2
		//       b) chi2<100
		//  2) stopmu flash
		//       a) if entering point is anode/cathode and (peratio>0.2 && chi2<40)
		//  2) thrumu ?

		for ( auto& candidate : candidates ) {
			// get flash result
			const flashana::FlashMatch_t& match = matches.at( candidate.qcluster_idx );
			// get chi2 and ratio
			float chi2 = match.score; // (if using unfitted flash match)
			float totpe_measured = 0.;
			float totpe_hypothesis = 0.;
			for ( int ipmt=0; ipmt<32; ipmt++ ) {
				totpe_measured += m_flash_matcher.FlashArray().at( match.flash_id ).pe_v.at(ipmt);
				totpe_hypothesis += match.hypothesis.at(ipmt);
			}
			float peratio = totpe_hypothesis/totpe_measured;

			bool accept = false;
			if ( candidate.source_type==0 ) {
				// untagged cluster object
				if ( peratio>0.1 && chi2<100 ) {
					accept = true;
				}
			}
			else if ( candidate.source_type==1 ) {
				// through-going muons
				if ( peratio>0.5 && peratio<1.5 && chi2 < 15 )
					accept = true;
			}
			else if ( candidate.source_type==2 ) {
				// stop mu
				// what was the end point that created this muon?
				if ( peratio>0.5 && peratio<1.5 && chi2 < 15 )
					accept = true;				
			}

			if ( accept ) { 
				int idx_source_type = candidate.source_type;
				int idx_source = candidate.source_idx;
				cluster_indices.insert( std::make_pair<int,int>( std::move(idx_source_type), std::move(idx_source) ) );
			}
		}
	}//end of selection function

	larcv::ROI FlashROIMatching::GenerateClusterROI( const std::vector<larcv::Pixel2DCluster>& cluster, const std::vector<larcv::Image2D>& img_v ) {

		larcv::ROI roi;

		for (size_t p=0; p<cluster.size(); p++ ) {
			const larcv::ImageMeta& meta = img_v.at(p).meta();
			const larcv::Pixel2DCluster& plane_cluster = cluster.at(p);
			std::vector<int> extrema(4); // low col, high col, low row, high row
			extrema.at(0) = plane_cluster.front().X();
			extrema.at(1) = plane_cluster.front().X();
			extrema.at(2) = plane_cluster.front().Y();
			extrema.at(3) = plane_cluster.front().Y();			

			for ( auto const& pixel : plane_cluster ) {
			  if ( extrema[0]>(int)pixel.X() ) extrema[0] = (int)pixel.X();
			  if ( extrema[1]<(int)pixel.X() ) extrema[1] = (int)pixel.X();
			  if ( extrema[2]>(int)pixel.Y() ) extrema[2] = (int)pixel.Y();
			  if ( extrema[3]<(int)pixel.Y() ) extrema[3] = (int)pixel.Y();
			}

			int drow = extrema[3]-extrema[2];
			int dcol = extrema[1]-extrema[0];

			// origin in upper corner
			float height = drow*meta.pixel_height();
			float width  = dcol*meta.pixel_width();
			float origin_x = meta.pos_x( extrema[0] );
			float origin_y = meta.pos_y( extrema[2] );

			larcv::ImageMeta bb( width, height, drow, dcol, origin_x, origin_y, (larcv::PlaneID_t)p );
			//std::cout << " defining ROI with origin(" << origin_x << "," << origin_y << ") BB: " << bb.dump() << std::endl;
			roi.AppendBB( bb );
		}

		return roi;
	}

}
