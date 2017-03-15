#include "TaggerFlashMatchAlgo.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"

namespace larlitecv {

	TaggerFlashMatchAlgo::TaggerFlashMatchAlgo()
	 : m_pmtweights("geoinfo.root") {}

	TaggerFlashMatchAlgo::TaggerFlashMatchAlgo( TaggerFlashMatchAlgoConfig& config )
	 : m_config(config), m_pmtweights("geoinfo.root") {
	 	setVerbosity( m_config.verbosity );
	}

	std::vector<larcv::ROI> TaggerFlashMatchAlgo::FindFlashMatchedContainedROIs( const std::vector<TaggerFlashMatchData>& inputdata, 
	  const std::vector<larlite::event_opflash*>& opflashes_v ) {

		// get the in-beam flashes
		std::vector<flashana::Flash_t> data_flashana = GetInTimeFlashana( opflashes_v );

		// merge objects [later if needed]

		// choose contained candidates
		std::vector<int> passes_containment( inputdata.size(), 0 );
		ChooseContainedCandidates( inputdata, passes_containment ); // simply check then 3D bounding box

		// flash-match candidates
		std::vector<int> passes_flashmatch( inputdata.size(), 0 );
		//std::vector<flashana::QCluster_t> qclusters_v = GenerateQClusters( inputdata );
		//ChooseInTimeFlashMatchedCandidates( qclusters_v, data_flashana, passes_flashmatch );

		std::vector<larcv::ROI> roi_v;

		for ( size_t i=0; i<inputdata.size(); i++) {
			if ( passes_containment[i] || passes_flashmatch[i] ) {
				// Make ROI
			}
		}

		// Merge ROI

		return roi_v;
	}

	flashana::QCluster_t TaggerFlashMatchAlgo::GenerateQCluster( const TaggerFlashMatchData& data ) {
		const larlite::track& track = data.m_track3d;
 
		flashana::QCluster_t qcluster;

		for ( int i=0; i<(int)track.NumberTrajectoryPoints()-1; i++ ) {

			const TVector3& pt     = track.LocationAtPoint(i);
			const TVector3& nextpt = track.LocationAtPoint(i+1);
			const TVector3& ptdir  = track.DirectionAtPoint(i);
			float dist = 0.;
			for ( int v=0; v<3; v++ ) {
				float dx = nextpt[v]-pt[v];
				dist += dx*dx;
			}
			dist = sqrt(dist);

			int nsteps = dist/m_config.qcluster_stepsize;
			if ( fabs( nsteps*m_config.qcluster_stepsize - dist )>0.01 ) {
				nsteps+=1;
			}

			float step = dist/float(nsteps);

			for (int istep=0; istep<nsteps; istep++) {
				float pos[3] = {0};
				for (int v=0; v<3; v++)
					pos[v] = pt[v] + step*ptdir[v];
				flashana::QPoint_t qpt( pos[0], pos[1], pos[2], step*m_config.MeV_per_cm );
				qcluster.emplace_back( std::move(qpt) );
			}
		}

		return qcluster;
	}

	std::vector< flashana::QCluster_t > TaggerFlashMatchAlgo::GenerateQClusters( const std::vector<TaggerFlashMatchData>& inputdata ) {
		std::vector< flashana::QCluster_t> qcluster_v;
		for ( auto const& data : inputdata ) {
			flashana::QCluster_t qcluster = GenerateQCluster( data );
			qcluster_v.emplace_back( std::move(qcluster) );
		}
		return qcluster_v;
	}

	flashana::Flash_t TaggerFlashMatchAlgo::GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster ) {
    const flashana::QLLMatch* matchalgo = (flashana::QLLMatch*)m_flash_matcher.GetAlgo( flashana::kFlashMatch );
    flashana::Flash_t unfitted_hypothesis = matchalgo->GetEstimate( qcluster );
    return unfitted_hypothesis;
	}

	std::vector<larlite::opflash> TaggerFlashMatchAlgo::SelectInTimeOpFlashes( const std::vector<larlite::event_opflash*>& opflashes_v ) {

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
			std::cout << "TaggerFlashMatchAlgo::SelectInTimeFlashes. Found " << beam_flashes.size() << " flashes." << std::endl;

		return beam_flashes;
	}

	void TaggerFlashMatchAlgo::GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange ) {

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

	flashana::Flash_t TaggerFlashMatchAlgo::MakeDataFlash( const larlite::opflash& flash ) {

    // get geometry info
		const larutil::Geometry* geo = ::larutil::Geometry::GetME();

		// provide mean and basic width
		float wire_mean;
		std::vector<float> wire_range;
		GetFlashCenterAndRange( flash, wire_mean, wire_range );
		if ( m_verbosity>0 )
			std::cout << "Flash: wire_mean=" << wire_mean << " wire_range=[" << (int)wire_range[0] << "," << (int)wire_range[1] << "]" << std::endl;

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
		f.idx = 0;

		return f;
	}

	std::vector<flashana::Flash_t> TaggerFlashMatchAlgo::MakeDataFlashes( const std::vector<larlite::opflash>& opflash_v ) {
		std::vector<flashana::Flash_t> flashes_v;
		for ( auto const& opflash : opflash_v ) {
			flashana::Flash_t flash = MakeDataFlash( opflash );
			flashes_v.emplace_back( std::move(flash) );
		}
		return flashes_v;
	}

  std::vector<flashana::Flash_t> TaggerFlashMatchAlgo::GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v ) {

  	std::vector<larlite::opflash> intime_opflashes_v = SelectInTimeOpFlashes( opflashes_v );

  	std::vector<flashana::Flash_t> intime_flashana_v = MakeDataFlashes( intime_opflashes_v );

  	return intime_flashana_v;
  }

  bool TaggerFlashMatchAlgo::IsClusterContained( const TaggerFlashMatchData& data ) {
    // simple selection
    // we make a cut on fiducial volume as a function of x

  	// first we need the bounding box of the points
  	float bb[3][2] = {0}; // extremes in each dimension
  	float extrema[3][2][3] = {-1.0e6 }; // each extrema point stored here. (dim,min/max,xyz)
  	for ( int v=0; v<3; v++) {
   		bb[v][0] = 1.0e6;
   		bb[v][1] = -1.0e6;
  	}

  	const larlite::track& track = data.m_track3d;
  	for ( size_t i=0; i<track.NumberTrajectoryPoints(); i++ ) {
  		const TVector3& xyz = track.LocationAtPoint(i);
  		for (int v=0; v<3; v++) {
  			// minvalue
  			if ( bb[v][0]>xyz[v] ) {
  				bb[v][0] = xyz[v];
  				for (int j=0; j<3; j++)
	  				extrema[v][0][j] = xyz[j];
  			}
  			// maxvalue
  			if ( bb[v][1]<xyz[v] ) {
  				bb[v][1] = xyz[v];
  				for (int j=0; j<3; j++)
	  				extrema[v][1][j] = xyz[j];
  			}
  		}
  	}



  	// we adjust for space charge effects
  	std::vector<double> delta_ymin = m_sce.GetPosOffsets( extrema[1][0][0], -118.0, extrema[1][0][2] );
  	std::vector<double> delta_ymax = m_sce.GetPosOffsets( extrema[1][1][0],  118.0, extrema[1][1][2] );
  	std::vector<double> delta_zmin = m_sce.GetPosOffsets( extrema[2][0][0], extrema[2][0][1], 0.0    );
  	std::vector<double> delta_zmax = m_sce.GetPosOffsets( extrema[2][1][0], extrema[2][1][1], 1037.0 );

  	if ( m_verbosity>=2 ) {
    	std::cout << "bounds: x=[" << bb[0][0] << "," << bb[0][1] << "] "
    	 << " y=[" << bb[1][0] << "," << bb[1][1] << "] "
  	   << " z=[" << bb[2][0] << "," << bb[2][1] << "] "
  	   << " dy=[" << delta_ymin[1] << "," << delta_ymax[1] << "] "
  	   << " dz=[" << delta_zmin[2] << "," << delta_zmax[2] << "] ";
  	}

  	// x extrema, not a function of the other dimensions
  	if ( bb[0][0]<m_config.FVCutX[0] || bb[0][1]>m_config.FVCutX[1] )
  		return false;
  	if ( (bb[1][0]-delta_ymin[1]) < m_config.FVCutY[0] || (bb[1][1]-delta_ymax[1]) > m_config.FVCutY[1] )
  		return false;
  	if ( (bb[2][0]-delta_zmin[2]) < m_config.FVCutZ[0] || (bb[2][1]-delta_zmax[2]) > m_config.FVCutZ[1] )
  		return false;

  	// we made it!
  	return true;
  }

  void TaggerFlashMatchAlgo::ChooseContainedCandidates( const std::vector<TaggerFlashMatchData>& inputdata, std::vector<int>& passes_containment ) {
    if ( passes_containment.size()!=inputdata.size() ) {
    	passes_containment.resize( inputdata.size(), 0 );
    }

    if ( m_verbosity>=2 ) {
    	std::cout << "[TaggerFlashMatchAlgo::ChooseContainedCandidates]" << std::endl;
    }
    for ( size_t i=0; i<inputdata.size(); i++ ) {

    	if ( m_verbosity>=2 ) {
    	  std::cout << " Candidate #" << i << ", ";
    	  if ( inputdata.at(i).m_type==TaggerFlashMatchData::kThruMu )
    	  	std::cout << "ThruMu";
    	  else if ( inputdata.at(i).m_type==TaggerFlashMatchData::kStopMu )
    	  	std::cout << "StopMu";
    	  else if ( inputdata.at(i).m_type==TaggerFlashMatchData::kUntagged )
    	  	std::cout << "Untagged";
    	  std::cout << ": ";
    	}

    	passes_containment[i] = ( IsClusterContained( inputdata.at(i) ) ) ? 1 : 0;

    	if ( m_verbosity>=2 ) {
      	if ( passes_containment[i] )
      		std::cout << " {contained}" << std::endl;
      	else
    	  	std::cout << "{uncontained}" << std::endl;
    	}
    }
  }

  float TaggerFlashMatchAlgo::InTimeFlashComparison( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster ) {
  	flashana::Flash_t flash_hypo = GenerateUnfittedFlashHypothesis( qcluster );

  	float smallest_chi2 = -1.0;
  	for ( auto const& intime_flash : intime_flashes_v ) {
  		int ndof = 0;
  		float chi2 = 0.;
  		for (size_t i=0; i<intime_flash.pe_v.size(); i++) {
  			if ( flash_hypo.pe_v.at(i)>0 ) {
    			float dpe = (intime_flash.pe_v.at(i) - flash_hypo.pe_v.at(i))/flash_hypo.pe_v.at(i);
  				ndof++;
  				chi2 += dpe;
  			}
  		}
  		chi2 /= ndof;
  		if ( smallest_chi2<0 || chi2<smallest_chi2 )
  			smallest_chi2 = chi2;
  	}
  	return smallest_chi2;
  }

  void TaggerFlashMatchAlgo::ChooseInTimeFlashMatchedCandidates( const std::vector<flashana::QCluster_t>& inputdata, 
  	const std::vector<flashana::Flash_t>& intime_flashes, std::vector<int>& passes_flashmatch ) {

  	if ( passes_flashmatch.size()!=inputdata.size() ) {
  		passes_flashmatch.resize( inputdata.size(), 0 );
  	}

  	for ( size_t q = 0; q<inputdata.size(); q++ ) {
  		const flashana::QCluster_t& qcluster = inputdata.at(q);
  		float chi2 = InTimeFlashComparison( intime_flashes, qcluster );
  		passes_flashmatch[q] = ( chi2 < m_config.flashmatch_chi2_cut ) ? 1 : 0;
  	}
  }

}