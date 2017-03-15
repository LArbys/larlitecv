#include "TaggerFlashMatchAlgo.h"

namespace larlitecv {

	TaggerFlashMatchAlgo::TaggerFlashMatchAlgo( TaggerFlashMatchAlgoConfig& config )
	 : m_config(config) {
	}

	std::vector<larcv::ROI> TaggerFlashMatchAlgo::FindFlashMatchedContainedROIs( const std::vector<TaggerFlashMatchData>& inputdata, 
	  const std::vector<larlite::event_opflash*>& opflashes_v ) {

		// get the in-beam flashes
		std::vector<flashana::Flash_t> data_flashana = GetInTimeFlashana( opflashes_v );

		// merge objects [later if needed]

		// choose contained candidates
		std::vector<int> passes_containment( inputdata.resize(), 0 );
		ChooseContainedCandidates( inputdata, passes_containment ); // simply check then 3D bounding box

		// flash-match candidates
		std::vector<flashana::QCluster_t> qclusters_v = GenerateQClusters( inputdata );
	}

	flashana::QCluster_t TaggerFlashMatchAlgo::GenerateQCluster( const TaggerFlashMatchData& data ) {
		const larlite::track& track = data.m_track3d;
 
		flashana::QCluster_t qcluster;

		for ( int i=0; i<track.NumberOfTrajectorPoints()-1; i++ ) {

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

	std::vector<flashana::QCluster_t> TaggerFlashMatchAlgo::GenerateQClusters( const TaggerFlashMatchData& inputdata ) {
		std::vector<flashana::QCluster_t> qcluster_v;
		for ( auto const& data : inputdata ) {
			flashana::QCluster_t qcluster = GenerateQCluster( data );
			qcluster_v.emplace_back( std::move(qcluster) );
		}
		return qcluster_v;
	}

	flashana::Flash_t TaggerFlashMatchAlgo::GenerateUnfittedFlashHypothesis( const larlite::QCluster_t& qcluster ) {
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
		const ::larutil::Geometry* geo = ::larutil::Geometry::GetME();

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

		return f;
	}

  std::vector<flashana::Flash_t> TaggerFlashMatchAlgo::GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v ) {

  	std::vector<larlite::opflash> intime_opflashes_v = SelectInTimeOpFlashes( opflashes_v );

  	std::vector<flashana::Flash_t> intime_flashana_v = MakeDataFlash( intime_opflashes_v );

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
  	for ( size_t i=0; i<track.NumberOfTrajectorPoints(); i++ ) {
  		const TVector3& xyz = track.LocationAtPoint(i);
  		for (int v=0; v<3; v++) {
  			// minvalue
  			if ( bb[v][0]>xyz[v] ) {
  				bb[v][0] = xyz[v];
  				for (int j=0; j<3; j++)
	  				extrema[v][0][i] = xyz[j];
  			}
  			// maxvalue
  			if ( bb[v][1]<xyz[v] ) {
  				bb[v][1] = xyz[v];
  				for (int j=0; j<3; j++)
	  				extrema[v][1][i] = xyz[j];
  			}
  		}
  	}

  	// x extrema, not a function of the other dimensions
  	if ( bb[0][0]<m_config.FVCutX[0] || bb[0][1]>m_config.FVCutX[1] )
  		return false;

  	// yz extrema a function of x
  	if ( bb[1][0]<m_config.FVCutCurveYtop.Eval( extrema[1][0][0] ) || b[1][1]>m_config.FVCutCurveYbottom.Eval( extrema[1][1][0] ) )
  		return false;
  	if ( bb[2][0]<m_config.FVCutCurveZtop.Eval( extrema[2][0][0] ) || b[2][1]>m_config.FVCutCurveZbottom.Eval( extrema[2][1][0] ) )
  		return false;

  	// we made it!
  	return true;
  }

  void TaggerFlashMatchAlgo::ChooseContainedCandidates( const TaggerFlashMatchData& inputdata, std::vector<int>& passes_containment ) {
    
  }

}