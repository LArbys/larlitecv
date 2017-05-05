#include "TaggerFlashMatchAlgo.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"

// larlitecv
#include "FlashMatchMetricMethods.h"

namespace larlitecv {

  TaggerFlashMatchAlgo::TaggerFlashMatchAlgo()
    : m_pmtweights("geoinfo.root") {}

  TaggerFlashMatchAlgo::TaggerFlashMatchAlgo( TaggerFlashMatchAlgoConfig& config )
    : m_config(config), m_pmtweights("geoinfo.root") {
    setVerbosity( m_config.verbosity );
    m_flash_matcher.Configure( m_config.m_flashmatch_config );

    const larutil::Geometry* geo = ::larutil::Geometry::GetME();
    for ( size_t opch=0; opch<32; opch++) {
      m_opch_from_opdet.insert( std::make_pair<int,int>( geo->OpDetFromOpChannel(opch), opch ) );
    }
  }

  
  std::vector<larcv::ROI> TaggerFlashMatchAlgo::FindFlashMatchedContainedROIs( const std::vector<TaggerFlashMatchData>& inputdata,
									       const std::vector<larlite::event_opflash*>& opflashes_v, std::vector<int>& flashdata_selected ) {

    // the output
    std::vector<larcv::ROI> roi_v;

    // clear state
    m_opflash_hypos.clear();
    m_min_chi2.clear();

    if (flashdata_selected.size()!=inputdata.size()) {
      flashdata_selected.resize( inputdata.size(), 0 );
    }

    // get the in-beam flashes
    std::vector<flashana::Flash_t> data_flashana = GetInTimeFlashana( opflashes_v );
    if ( data_flashana.size()==0 ) {
      return roi_v;
    }


    // merge objects [later if needed]

    // choose contained candidates
    m_passes_containment.resize( inputdata.size(), 0 );
    m_passes_flashmatch.resize( inputdata.size(), 0 );
    m_passes_totpe.resize( inputdata.size(), 0 );    
    m_passes_cosmicflash_ratio.resize( inputdata.size(), 0 );

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

      m_passes_containment[i] = ( IsClusterContained( inputdata.at(i) ) ) ? 1 : 0;

      if ( m_verbosity>=2 ) {
      	if ( m_passes_containment[i] )
          std::cout << " {contained}" << std::endl;
      	else
          std::cout << "{uncontained}" << std::endl;
      }

      // flash-match candidates
      flashana::QCluster_t qcluster = GenerateQCluster( inputdata.at(i) );

      // in-time flash matching
      float totpe_data = 0;
      float totpe_hypo = 0;
      m_passes_flashmatch[i] = ( DoesQClusterMatchInTimeFlash( data_flashana, qcluster, totpe_data, totpe_hypo )  ) ? 1 : 0;

      // totpe
      m_passes_totpe[i] = ( DoesTotalPEMatch( totpe_data, totpe_hypo ) ) ? 1 : 0;

      // cosmic versus in-time flash log-likelihood ratio
      float dchi2=0;
      m_passes_cosmicflash_ratio[i] = ( DoesQClusterMatchInTimeBetterThanCosmic( data_flashana, qcluster, inputdata[i], i, dchi2 ) ) ? 1 : 0;

      if ( m_verbosity>=2 ) {
        std::cout << "  ";

      	if ( m_passes_containment[i] )
          std::cout << " {contained}";
      	else
          std::cout << "{uncontained}";

        if ( m_passes_flashmatch[i]==1 ) {
          std::cout << " {in-time flash-matched}";
        }
        else {
          std::cout << " {not-matched}";
        }

	if ( m_passes_totpe[i]==1 )
	  std::cout << " {totpe matched}";
	else
	  std::cout << " {totpe fails}";

	if ( !m_passes_cosmicflash_ratio[i] )
	  std::cout << " {fails cosmic-flash match: dchi2=" << dchi2 << "}";
	else if ( dchi2!=0 )
	  std::cout << " {passes cosmic-flash match: dchi2=" << dchi2 << "}";
	else
	  std::cout << " { no flash end }";
	  
	
	if ( m_passes_flashmatch[i] && m_passes_containment[i] && m_passes_cosmicflash_ratio[i] && m_passes_totpe[i] )
          std::cout << " **PASSES**";
        std::cout << std::endl;
      }
    }//end of input data loop

    for ( size_t i=0; i<inputdata.size(); i++) {
      if ( m_passes_containment[i] && m_passes_flashmatch[i] && m_passes_cosmicflash_ratio[i] && m_passes_totpe[i] ) {
        // Make ROI
        flashdata_selected[i] = 1;
      }
    }

    return roi_v;
  }

  bool TaggerFlashMatchAlgo::didTrackPassContainmentCut( int itrack ) {
    if ( itrack<0 || itrack>=(int)m_passes_containment.size() )
      return false;
    return m_passes_containment[itrack];
  }

  bool TaggerFlashMatchAlgo::didTrackPassFlashMatchCut( int itrack ) {
    if ( itrack<0 || itrack>=(int)m_passes_flashmatch.size() )
      return false;
    return m_passes_flashmatch[itrack];
  }
  
  bool TaggerFlashMatchAlgo::didTrackPassCosmicFlashCut( int itrack ) {
    if ( itrack<0 || itrack>=(int)m_passes_cosmicflash_ratio.size() )
      return false;
    return m_passes_cosmicflash_ratio[itrack];
  }

  flashana::QCluster_t TaggerFlashMatchAlgo::GenerateQCluster( const TaggerFlashMatchData& data ) {
    const larlite::track& track = data.m_track3d;
    return GenerateQCluster(track);
  }
  
  flashana::QCluster_t TaggerFlashMatchAlgo::GenerateQCluster( const larlite::track& track ) {

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
        if ( opflash.TotalPE()<m_config.flashpe_thresh )
          continue;
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
      std::cout << "Flash: pe=" << flash.TotalPE() << " wire_mean=" << wire_mean << " wire_range=[" << (int)wire_range[0] << "," << (int)wire_range[1] << "]" << std::endl;

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

  larlite::opflash TaggerFlashMatchAlgo::MakeOpFlashFromFlash(const flashana::Flash_t& inFlash) {

    const larutil::Geometry* geo = ::larutil::Geometry::GetME();

    std::vector<double> PEperOpDet( inFlash.pe_v.size(), 0.0);
    for (unsigned int femch = 0; femch < inFlash.pe_v.size(); femch++) {
      unsigned int opdet = geo->OpDetFromOpChannel(femch);
      PEperOpDet[femch] = inFlash.pe_v[opdet];
    }

    larlite::opflash flash( inFlash.time, 0, 0, 0, PEperOpDet, false, false, 1.0, inFlash.y, inFlash.y_err, inFlash.z, inFlash.z_err );

    return flash;
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

  std::vector< std::vector<float> > TaggerFlashMatchAlgo::GetAABoundingBox( const larlite::track& track )  {

    std::vector< std::vector<float> > bb(3);
    for (int v=0; v<3; v++) {
      bb.at(v).resize(2,0.0);
    }

    float extrema[3][2][3] = {-1.0e6 }; // each extrema point stored here. (dim,min/max,xyz)
    for ( int v=0; v<3; v++) {
      bb[v][0] = 1.0e6;
      bb[v][1] = -1.0e6;
    }

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

    return bb;
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
    //std::vector<double> delta_ymin = m_sce.GetPosOffsets( extrema[1][0][0], -118.0, extrema[1][0][2] );
    //std::vector<double> delta_ymax = m_sce.GetPosOffsets( extrema[1][1][0],  118.0, extrema[1][1][2] );
    //std::vector<double> delta_zmin = m_sce.GetPosOffsets( extrema[2][0][0], extrema[2][0][1], 0.0    );
    //std::vector<double> delta_zmax = m_sce.GetPosOffsets( extrema[2][1][0], extrema[2][1][1], 1037.0 );
    std::vector<double> delta_ymin(3,0);
    std::vector<double> delta_ymax(3,0);
    std::vector<double> delta_zmin(3,0);
    std::vector<double> delta_zmax(3,0);


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

  float TaggerFlashMatchAlgo::InTimeFlashComparison( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster, float& best_totpe_data, float& tot_pe_hypo ) {
    flashana::Flash_t flash_hypo = GenerateUnfittedFlashHypothesis( qcluster );
    larlite::opflash opflash_hypo = MakeOpFlashFromFlash( flash_hypo );
    larlite::opflash gaus2d_hypo = larlitecv::GetGaus2DPrediction( opflash_hypo );
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();

    float smallest_chi2 = -1.0;
    tot_pe_hypo = 0.;
    std::vector< float > expectation(32,0);
    for ( auto const& intime_flash : intime_flashes_v ) {
      float ll = 0.;
      float tot_pe = 0.;
      for (size_t i=0; i<intime_flash.pe_v.size(); i++) {
        float observed = intime_flash.pe_v.at(i);
        float expected = flash_hypo.pe_v.at(i)*m_config.fudge_factor;
        tot_pe += observed;
        tot_pe_hypo += expected;
        if ( observed>0 && expected==0 ) {
	  if ( !m_config.use_gaus2d )
	    expected = 1.0e-3;
	  else
	    expected = gaus2d_hypo.PE( geo->OpDetFromOpChannel(i) )*m_config.fudge_factor;
	}
	expectation[i] = expected;
	
        float dpe = (expected-observed);
        if ( observed>0 )
          dpe += observed*( log( observed ) - log( expected ) );
        ll += 2.0*dpe;
      }
      ll /= 32.0;
      if ( smallest_chi2<0 || ll<smallest_chi2 ) {
        smallest_chi2 = ll;
	best_totpe_data = tot_pe;
      }
      std::cout << "  [observed] ";
      for (int ich=0; ich<32; ich++ ) {
        std::cout << std::setw(5) << (int)intime_flash.pe_v.at(  geo->OpDetFromOpChannel(ich) );
      }
      std::cout << " TOT=" << tot_pe << " LL=" << ll << std::endl;
    }
    tot_pe_hypo /= float( intime_flashes_v.size() );
    
    std::cout << "  [expected] ";
    for ( int ich=0; ich<32; ich++ ) {
      //std::cout << std::setw(5) << (int)(flash_hypo.pe_v.at(  geo->OpDetFromOpChannel(ich) )*m_config.fudge_factor);
      //std::cout << std::setw(5) << (int)(flash_hypo.pe_v.at( ich )*m_config.fudge_factor);
      std::cout << std::setw(5) << (int)expectation[ich];
    }
    std::cout << " TOT=" << tot_pe_hypo << " BestLL=" << smallest_chi2 << std::endl;
    
    // store flash info and chi2 for diagnostic info
    // we sneakily use this info again. also store it to study cut performance
    m_opflash_hypos.emplace_back(std::move(opflash_hypo));
    m_min_chi2.push_back( smallest_chi2 );

    return smallest_chi2;
  }
  
  bool TaggerFlashMatchAlgo::DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster, float& totpe_data, float& totpe_hypo ) {

    float chi2 = InTimeFlashComparison( intime_flashes_v, qcluster, totpe_data, totpe_hypo );

    // our cut needs to be position dependent, so we need a position. we get the q-weighted mean and the intervals.
    float totw = 0.;
    float mean[3] = {0};
    float min_x = 1.0e6;
    float max_x = -1.0e6;
    for ( auto const& qpt : qcluster ) {
      mean[0] += qpt.x*qpt.q;
      mean[1] += qpt.y*qpt.q;
      mean[2] += qpt.z*qpt.q;
      totw += qpt.q;
      if ( qpt.x < min_x ) min_x = qpt.x;
      if ( qpt.x > max_x ) max_x = qpt.x;
    }
    if ( totw==0 ) {
      return false;
    }
    else {
      for (int i=0; i<3; i++) mean[i] /= totw;
    }

    //  still a dumb cut
    if ( mean[0] < 100.0 && chi2< m_config.flashmatch_chi2_cut*1.5 )
      return true;
    else if ( mean[0]>100.0 && chi2<m_config.flashmatch_chi2_cut )
      return true;
    else
      return false;
    
  }

  bool TaggerFlashMatchAlgo::DoesTotalPEMatch( float totpe_data, float totpe_hypo )  {
    float frac_diff = fabs(totpe_data-totpe_hypo)/totpe_hypo;
    if ( frac_diff > m_config.totpe_sigma_cut )
      return false;
    return true;
  }
  
  bool TaggerFlashMatchAlgo::DoesQClusterMatchInTimeBetterThanCosmic( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster,
								      const larlitecv::TaggerFlashMatchData& taggertrack, const int trackidx, float& dchi2 ) {
    if ( !taggertrack.hasStartFlash() && !taggertrack.hasEndFlash() ) {
      // No flash associated to this track.
      return true;
    }

    // check that the matching flash is not the intime flash
    
    float intime_min_chi2 = m_min_chi2[trackidx]; // we cheat and use this stored value instead of recalculating

    float totpe_cosmicflash = 0.;
    float totpe_hypoflash = 0;
    flashana::Flash_t flash_hypo = GenerateUnfittedFlashHypothesis( qcluster );
    larlite::opflash opflash_hypo = MakeOpFlashFromFlash( flash_hypo );
    
    float cosmicflash_start_chi2 = 1.0e6;    
    if ( taggertrack.hasStartFlash() ) {
      std::vector< larlite::opflash > cosmicflash_start_v;
      float start_tick = taggertrack.m_pstart_flash->Time()/m_config.us_per_tick;
      float start_ly = m_config.fudge_factor;
      if ( start_tick<m_config.beam_tick_range[0] || start_tick>m_config.beam_tick_range[1] )
	start_ly = m_config.fudge_factor_cosmic;
      cosmicflash_start_v.push_back( *(taggertrack.m_pstart_flash ) );
      cosmicflash_start_chi2 = larlitecv::CalculateFlashMatchChi2( cosmicflash_start_v, opflash_hypo, totpe_cosmicflash, totpe_hypoflash, start_ly, m_config.use_gaus2d, false );      
    }
    
    float cosmicflash_end_chi2 = 1.0e6;        
    if ( taggertrack.hasEndFlash() ) {
      std::vector< larlite::opflash > cosmicflash_end_v;      
      float end_ly = m_config.fudge_factor;      
      float end_tick = taggertrack.m_pend_flash->Time()/m_config.us_per_tick;
      if ( end_tick<m_config.beam_tick_range[0] || end_tick>m_config.beam_tick_range[1] )
	end_ly = m_config.fudge_factor_cosmic; // use cosmic disc.
      cosmicflash_end_v.push_back( *(taggertrack.m_pend_flash ) );
      cosmicflash_end_chi2 = larlitecv::CalculateFlashMatchChi2( cosmicflash_end_v, opflash_hypo, totpe_cosmicflash, totpe_hypoflash, end_ly, m_config.use_gaus2d, false );
    }
    
    float cosmicflash_chi2 = ( cosmicflash_start_chi2<cosmicflash_end_chi2 ) ? cosmicflash_start_chi2 : cosmicflash_end_chi2;
    
    dchi2 = intime_min_chi2 - cosmicflash_chi2;
    
    if ( cosmicflash_chi2 < intime_min_chi2 ) {
      // matches cosmic flash better
      return false;
    }
    else {
      return true;
    }
    
    return true; // never gets here
  }
  
}
