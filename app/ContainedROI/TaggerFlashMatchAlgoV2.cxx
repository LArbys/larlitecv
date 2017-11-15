
#include "TaggerFlashMatchAlgoV2.h"

#include <sstream>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"

namespace larlitecv {

  TaggerFlashMatchAlgoV2::TaggerFlashMatchAlgoV2() {
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << " Default constructor only declared for ROOT dict. purposes. Should not be used." << std::endl;
    throw std::runtime_error(msg.str());
  }

  TaggerFlashMatchAlgoV2::TaggerFlashMatchAlgoV2( TaggerFlashMatchAlgoConfig& config )
    : m_config(config),m_genflashmatch(config.genflashmatch_cfg) {
    // set verbosity
    setVerbosity( m_config.verbosity );
    // opdet-opch map
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();
    for ( size_t opch=0; opch<32; opch++) {
      m_opch_from_opdet.insert( std::make_pair<int,int>( geo->OpDetFromOpChannel(opch), opch ) );
    }
  }

  std::vector<larcv::ROI> TaggerFlashMatchAlgoV2::FindFlashMatchedContainedROIs( const std::vector<TaggerFlashMatchData>& inputdata,
									       const std::vector<larlite::event_opflash*>& opflashes_v, std::vector<int>& flashdata_selected ) {

    // the output
    std::vector<larcv::ROI> roi_v;

    // geometry
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();    

    // clear record of decisions
    if (flashdata_selected.size()!=inputdata.size()) {
      flashdata_selected.resize( inputdata.size(), 0 );
    }

    // get the in-beam flashes
    std::vector<flashana::Flash_t> data_flashana       = m_genflashmatch.GetInTimeFlashana( opflashes_v );
    std::vector<flashana::Flash_t> cosmicdata_flashana = m_genflashmatch.GetCosmicFlashana( opflashes_v ); 
    if ( data_flashana.size()==0 ) {
      return roi_v;
    }

    // merge objects [later if needed]

    // prepare result containers
    m_passes_containment.resize( inputdata.size(), 0 );
    m_passes_flashmatch.resize( inputdata.size(), 0 );
    m_passes_totpe.resize( inputdata.size(), 0 );
    m_passes_cosmicflash_ratio.resize( inputdata.size(), 0 );

    m_opflash_hypos.clear();
    m_containment_dwall.clear();
    m_cosmicflash_ratio_dchi.resize( inputdata.size(), -1 );
    m_min_chi2.resize( inputdata.size(), -1 );
    m_totpe_peratio.resize( inputdata.size(),- 1 );
    m_cosmicflash_ratio_dchi.resize( inputdata.size(), -1 );
    m_genflashmatch.getFlashMatchManager().Reset();


    std::cout << "TPC List ------ -------------------------------------" << std::endl;
    for ( size_t itrack=0; itrack<inputdata.size(); itrack++ ) {
      const TaggerFlashMatchData& track = inputdata[itrack];
      float xmin = track.m_track3d.Vertex()[0];
      float xmax = track.m_track3d.End()[0];
      std::cout << "  #" << itrack << " xmin=" << xmin << " xmax=" << xmax << std::endl;
    }
    std::cout << "-------------------------------------------------------" << std::endl;
    
    std::cout << "Flash List ------ -------------------------------------" << std::endl;
    for ( size_t iflash=0; iflash<data_flashana.size(); iflash++) {
      float usec = data_flashana[iflash].time;
      std::cout << "  #" << iflash<< ": idx= " << data_flashana[iflash].idx << " pe=" << data_flashana[iflash].TotalPE() << " time=" << data_flashana[iflash].time << std::endl;
    }
    for ( size_t iflash=0; iflash<cosmicdata_flashana.size(); iflash++) {
      float usec = cosmicdata_flashana[iflash].time;
      std::cout << "  #" << iflash<< ": idx= " << cosmicdata_flashana[iflash].idx << " pe=" << cosmicdata_flashana[iflash].TotalPE() << " time=" << cosmicdata_flashana[iflash].time << std::endl;
    }
    std::cout << "-------------------------------------------------------" << std::endl;

    // setup the qclusters
    setupQClusters( inputdata );
    
    // in-time tests: 
    // (1) get z-mean and z-FWHM bound for flash
    // (2) keep those tracks within bounds

    std::vector<FlashRange_t> flashrange_v = getFlashRange( data_flashana );
    
    
    std::vector<flashana::FlashMatch_t> results = m_genflashmatch.getFlashMatchManager().Match();
    // print results and make flash-hypothesis
    m_intime_bestflash_hypo_v.resize(inputdata.size());
    m_intime_bestflash_chi2_v.resize(inputdata.size(),-1.0);
    for (auto& flash : m_intime_bestflash_hypo_v)
      flash.pe_v.resize(32,0.);
    std::cout << "In-time FlashMatch Result -----------------------------------" << std::endl;
    for ( auto const& flashmatch : results ) {
      flashana::Flash_t& flashhypo = m_intime_bestflash_hypo_v[flashmatch.tpc_id];
      flashhypo.pe_v.resize(32);
      float x_offset = flashmatch.tpc_point.x;
      flashana::QCluster_t best_cluster( m_qcluster_v[flashmatch.tpc_id] );
      double min_x = 1e9;
      for (size_t i = 0; i < best_cluster.size(); ++i) {
	auto const &pt = best_cluster[i];
	if ( pt.y<-116.5 || pt.y>116.5 || pt.z<0.0 || pt.z>1036.8 )
	  continue;
	if (pt.x < min_x) { min_x = pt.x; }
      }
      for ( auto& pt : best_cluster ) {
	pt.x += x_offset-min_x;
      }
      ((flashana::BaseFlashMatch*)m_genflashmatch.getFlashMatchManager().GetAlgo( flashana::kFlashMatch ))->FillEstimate( best_cluster, flashhypo );
      std::cout << "  score: " << flashmatch.score << "  flashid=" << flashmatch.flash_id << " tpcid=" << flashmatch.tpc_id
		<< " x-offset=" << x_offset << " x-min=" << min_x << " dx=" << x_offset-min_x
		<< std::endl;
      m_intime_bestflash_chi2_v[flashmatch.tpc_id] = -2.0*log10(flashmatch.score);
    }
    std::cout << "-------------------------------------------------------------" << std::endl;

    // --------------------------------------------------------------------------------
    // out-of-time flash tests: load in-time flash(es) and non-extended flashes
    bool extend_tracks = true;
    setupFlashMatchInterface( data_flashana, cosmicdata_flashana, extend_tracks );
    std::vector<flashana::FlashMatch_t> results_outoftime = m_genflashmatch.getFlashMatchManager().Match();
    // print results and make flash-hypothesis
    m_cosmic_bestflash_hypo_v.resize(inputdata.size());
    m_cosmic_bestflash_idx_v.resize(inputdata.size());
    m_cosmic_bestflash_chi2_v.resize(inputdata.size(),-1);    
    for (auto& flash : m_cosmic_bestflash_hypo_v)
      flash.pe_v.resize(32,0.);    
    std::cout << "Out-of-time FlashMatch Result -----------------------------------" << std::endl;
    for ( auto const& flashmatch : results_outoftime ) {
      flashana::Flash_t& flashhypo = m_cosmic_bestflash_hypo_v[flashmatch.tpc_id];
      flashhypo.pe_v.resize(32);      
      float x_offset = flashmatch.tpc_point.x;
      flashana::QCluster_t best_cluster( m_qcluster_extended_v[flashmatch.tpc_id] );
      double min_x = 1e9;
      for (size_t i = 0; i < best_cluster.size(); ++i) {
	auto const &pt = best_cluster[i];
	if ( pt.y<-116.5 || pt.y>116.5 || pt.z<0.0 || pt.z>1036.8 )
	  continue;	
	if (pt.x < min_x) { min_x = pt.x; }
      }      
      for ( auto& pt : best_cluster ) {
	pt.x += x_offset-min_x;
      }
      ((flashana::BaseFlashMatch*)m_genflashmatch.getFlashMatchManager().GetAlgo( flashana::kFlashMatch ))->FillEstimate( best_cluster, flashhypo );

      // cosmic disc correction from chris
      for (size_t ich=0; ich<flashhypo.pe_v.size(); ich++ ) {
	float pe = flashhypo.pe_v.at(ich);
	if ( pe < 60.0 )
	  pe *= 0.424;
	else if ( pe > 60.0 )
	  pe *= 0.354;
	flashhypo.pe_v[ich] = pe;
      }
      m_cosmic_bestflash_idx_v[flashmatch.tpc_id] = flashmatch.flash_id;
      std::cout << "  score: " << flashmatch.score << "  flashid=" << flashmatch.flash_id << " tpcid=" << flashmatch.tpc_id
		<< " x-offset=" << x_offset << " x-min=" << min_x << " dx=" << x_offset-min_x
		<< std::endl;
		
      m_cosmic_bestflash_chi2_v[flashmatch.tpc_id] = -2.0*log10(flashmatch.score);
    }
    std::cout << "-----------------------------------------------------------------" << std::endl;

    //m_genflashmatch.getFlashMatchManager().PrintFullResult();
    
    // choose contained candidates    
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
      
      // in-time flash matching
      float totpe_data = 0;
      float totpe_hypo = 0;
      float min_chi2   = -1;
      m_passes_flashmatch[i] = ( DoesQClusterMatchInTimeFlash( data_flashana, inputdata[i], totpe_data, totpe_hypo, min_chi2 )  ) ? 1 : 0;

      // totpe
      m_passes_totpe[i] = ( m_genflashmatch.DoesTotalPEMatch( totpe_data, totpe_hypo ) ) ? 1 : 0;

      flashana::QCluster_t qcluster;
      m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qcluster, inputdata[i].m_track3d, 100.0, 15.0 );
      
      larlite::opflash ophypo = m_genflashmatch.MakeOpFlashFromFlash( m_intime_bestflash_hypo_v[i] );
      m_opflash_hypos.emplace_back( std::move(ophypo) );

      // store values
      
      // intime chi2
      m_min_chi2[i] = m_intime_bestflash_chi2_v[i];
      // intime chi2 - cosmic chi2
      if ( m_min_chi2[i]<0 ) {
	m_cosmicflash_ratio_dchi[i] = -10 - m_cosmic_bestflash_chi2_v[i];
      }
      else {
	m_cosmicflash_ratio_dchi[i] = m_min_chi2[i]-m_cosmic_bestflash_chi2_v[i];
      }
      // intime pe ratio
      if ( m_min_chi2[i]<0 ) {
	m_totpe_peratio[i] = -1;
      }
      else {
	m_totpe_peratio[i] = data_flashana.front().TotalPE() / m_intime_bestflash_hypo_v[i].TotalPE();
      }
	
      if ( m_verbosity>=2  ) {
	// data in-time flash
	std::cout << "  [observed] ";
	for (int ich=0; ich<32; ich++ ) {
	  std::cout << std::setw(5) << (int)data_flashana.front().pe_v.at(  ich );
	}
	std::cout << " TOT=" << data_flashana.front().TotalPE() << " CHI2=" << "XXX" << std::endl;
    
	std::cout << "  [in-time ] ";
	for ( int ich=0; ich<32; ich++ ) {
	  //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at(  geo->OpDetFromOpChannel(ich) )*m_config.fudge_factor);
	  //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at( ich )*m_config.fudge_factor);
	  std::cout << std::setw(5) << (int)m_intime_bestflash_hypo_v[i].pe_v[ich];
	}
	std::cout << " TOT=" << m_intime_bestflash_hypo_v[i].TotalPE()
		  << " BestCHI2=" << m_intime_bestflash_chi2_v[i] << std::endl;

	std::cout << "  [bestdata] ";
	for ( int ich=0; ich<32; ich++ ) {
	  //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at(  geo->OpDetFromOpChannel(ich) )*m_config.fudge_factor);
	  //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at( ich )*m_config.fudge_factor);
	  std::cout << std::setw(5) << (int)cosmicdata_flashana[m_cosmic_bestflash_idx_v[i]].pe_v[ich];
	}
	std::cout << " TOT=" << cosmicdata_flashana[m_cosmic_bestflash_idx_v[i]].TotalPE()  << std::endl;

	std::cout << "  [besthypo] ";
	for ( int ich=0; ich<32; ich++ ) {
	  //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at(  geo->OpDetFromOpChannel(ich) )*m_config.fudge_factor);
	  //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at( ich )*m_config.fudge_factor);
	  std::cout << std::setw(5) << (int)m_cosmic_bestflash_hypo_v[i].pe_v[ich];
	}
	std::cout << " TOT=" << m_cosmic_bestflash_hypo_v[i].TotalPE()
		  << " BestCHI2=" << m_cosmic_bestflash_chi2_v[i] << std::endl;
	
      }
	
      
      // cosmic versus in-time flash log-likelihood ratio
      float dchi2=0;
      m_passes_cosmicflash_ratio[i] = ( DoesQClusterMatchInTimeBetterThanCosmic( cosmicdata_flashana, inputdata[i], min_chi2, dchi2 ) ) ? 1 : 0;
      m_cosmicflash_ratio_dchi.push_back( dchi2 );

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

      //m_qcluster_v.emplace_back( std::move(qcluster) );
      //m_opflash_hypos.emplace_back( std::move(ophypo) );
      
    }//end of input data loop

    for ( size_t i=0; i<inputdata.size(); i++) {
      if ( m_passes_containment[i] && m_passes_flashmatch[i] && m_passes_cosmicflash_ratio[i] && m_passes_totpe[i] ) {
        // Make ROI
        flashdata_selected[i] = 1;
      }
    }


    // std::cout << "QCLUSTER List ------ -------------------------------------" << std::endl;
    // for ( size_t itrack=0; itrack<inputdata.size(); itrack++ ) {
    //   const flashana::QCluster_t& qcluster = m_qcluster_v[itrack];
    //   std::cout << "  #" << itrack << ": " << qcluster.size() << " points" << std::endl;
    //   for ( auto const& qpt : qcluster ) {
    // 	std::cout << "      (" << qpt.x << "," << qpt.y << "," << qpt.z << ") q=" << qpt.q << std::endl;
    //   }
    // }
    // std::cout << "-------------------------------------------------------" << std::endl;

    
    return roi_v;
  }

  bool TaggerFlashMatchAlgoV2::didTrackPassContainmentCut( int itrack ) {
    if ( itrack<0 || itrack>=(int)m_passes_containment.size() )
      return false;
    return m_passes_containment[itrack];
  }

  bool TaggerFlashMatchAlgoV2::didTrackPassFlashMatchCut( int itrack ) {
    if ( itrack<0 || itrack>=(int)m_passes_flashmatch.size() )
      return false;
    return m_passes_flashmatch[itrack];
  }

  bool TaggerFlashMatchAlgoV2::didTrackPassCosmicFlashCut( int itrack ) {
    if ( itrack<0 || itrack>=(int)m_passes_cosmicflash_ratio.size() )
      return false;
    return m_passes_cosmicflash_ratio[itrack];
  }

  flashana::QCluster_t TaggerFlashMatchAlgoV2::GenerateQCluster( const TaggerFlashMatchData& data ) {
    // we use the generalflashmatch algo qcluster maker, which will extend the ends if needed.
    flashana::QCluster_t qclust;
    const larlite::track& track = data.m_track3d;
    m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qclust, track, 1000.0, 10.0 );
    return qclust;
  }

  bool TaggerFlashMatchAlgoV2::IsClusterContained( const TaggerFlashMatchData& data ) {
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

    float dwall[6];
    dwall[0] = fabs(bb[0][0]-m_config.FVCutX[0]);
    dwall[1] = fabs(bb[0][1]-m_config.FVCutX[1]);
    dwall[2] = fabs(bb[1][0]-m_config.FVCutY[0]);
    dwall[3] = fabs(bb[1][1]-m_config.FVCutY[1]);    
    dwall[4] = fabs(bb[2][0]-m_config.FVCutZ[0]);
    dwall[5] = fabs(bb[2][1]-m_config.FVCutZ[1]);

    
    // x extrema, not a function of the other dimensions
    if ( bb[0][0]<m_config.FVCutX[0] || bb[0][1]>m_config.FVCutX[1] ) {
      if ( dwall[0]<dwall[1] )
	m_containment_dwall.push_back( dwall[0] );
      else
	m_containment_dwall.push_back( dwall[1] );
      return false;
    }
    if ( (bb[1][0]-delta_ymin[1]) < m_config.FVCutY[0] || (bb[1][1]-delta_ymax[1]) > m_config.FVCutY[1] ) {
      if ( dwall[2]<dwall[3] )
	m_containment_dwall.push_back( dwall[2] );
      else
	m_containment_dwall.push_back( dwall[3] );
      return false;
    }
    if ( (bb[2][0]-delta_zmin[2]) < m_config.FVCutZ[0] || (bb[2][1]-delta_zmax[2]) > m_config.FVCutZ[1] ) {
      if ( dwall[4]<dwall[5] )
	m_containment_dwall.push_back( dwall[4] );
      else
	m_containment_dwall.push_back( dwall[5] );      
      return false;
    }

    float mindwall = 1e6;
    for ( int i=0; i<6; i++) {
      if ( mindwall>dwall[i] )
	mindwall = dwall[i];
    }

    m_containment_dwall.push_back( mindwall );
    
    // we made it!
    return true;
  }

  void TaggerFlashMatchAlgoV2::ChooseContainedCandidates( const std::vector<TaggerFlashMatchData>& inputdata, std::vector<int>& passes_containment ) {
    if ( passes_containment.size()!=inputdata.size() ) {
      passes_containment.resize( inputdata.size(), 0 );
    }

    if ( m_verbosity>=2 ) {
      std::cout << "[TaggerFlashMatchAlgoV2::ChooseContainedCandidates]" << std::endl;
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

  bool TaggerFlashMatchAlgoV2::DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v,
							   const larlitecv::TaggerFlashMatchData& taggertrack,
							   float& totpe_data, float& totpe_hypo, float& min_chi2 ) {
    flashana::QCluster_t qcluster;
    // generate qcluster using genflash. Do not extend as this is the in-time hypothesis
    m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qcluster, taggertrack.m_track3d, 1000.0, 10.0 );
    
    bool matches = m_genflashmatch.DoesQClusterMatchInTimeFlash( intime_flashes_v, qcluster, totpe_data, totpe_hypo, min_chi2 );
    
    return matches;
  }
  
  bool TaggerFlashMatchAlgoV2::DoesQClusterMatchInTimeBetterThanCosmic( const std::vector<flashana::Flash_t>& cosmic_flashes_v,
								      const larlitecv::TaggerFlashMatchData& taggertrack,
								      const float intime_min_chi2, float& dchi2 ) {
    if ( !taggertrack.hasStartFlash() && !taggertrack.hasEndFlash() ) {
      // No flash associated to this track.
      return true;
    }

    flashana::QCluster_t qcluster;
    // generate qcluster using genflash.
    m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qcluster, taggertrack.m_track3d, 1000.0, 10.0 );

    // check that the matching flash is not the intime flash

    float totpe_cosmicflash = 0.;
    float totpe_hypoflash = 0;
    dchi2 = 0;
    int bestflash_idx = 0;

    m_genflashmatch.GetFlashMinChi2( cosmic_flashes_v, qcluster, 1, totpe_cosmicflash, totpe_hypoflash, dchi2, bestflash_idx );

    dchi2 -= intime_min_chi2;
    
    // float cosmicflash_start_chi2 = 1.0e6;
    // if ( taggertrack.hasStartFlash() ) {
    //   std::vector< larlite::opflash > cosmicflash_start_v;
    //   float start_tick = taggertrack.m_pstart_flash->Time()/m_config.us_per_tick;
    //   float start_ly = m_config.fudge_factor;
    //   if ( start_tick<m_config.beam_tick_range[0] || start_tick>m_config.beam_tick_range[1] )
    //     start_ly = m_config.fudge_factor_cosmic;
    //   cosmicflash_start_v.push_back( *(taggertrack.m_pstart_flash ) );
    //   cosmicflash_start_chi2 = larlitecv::CalculateFlashMatchChi2( cosmicflash_start_v, opflash_hypo, totpe_cosmicflash, totpe_hypoflash, start_ly, m_config.use_gaus2d, false );
    // }

    // float cosmicflash_end_chi2 = 1.0e6;
    // if ( taggertrack.hasEndFlash() ) {
    //   std::vector< larlite::opflash > cosmicflash_end_v;
    //   float end_ly = m_config.fudge_factor;
    //   float end_tick = taggertrack.m_pend_flash->Time()/m_config.us_per_tick;
    //   if ( end_tick<m_config.beam_tick_range[0] || end_tick>m_config.beam_tick_range[1] )
    //     end_ly = m_config.fudge_factor_cosmic; // use cosmic disc.
    //   cosmicflash_end_v.push_back( *(taggertrack.m_pend_flash ) );
    //   cosmicflash_end_chi2 = larlitecv::CalculateFlashMatchChi2( cosmicflash_end_v, opflash_hypo, totpe_cosmicflash, totpe_hypoflash, end_ly, m_config.use_gaus2d, false );
    // }

    // float cosmicflash_chi2 = ( cosmicflash_start_chi2<cosmicflash_end_chi2 ) ? cosmicflash_start_chi2 : cosmicflash_end_chi2;

    // dchi2 = intime_min_chi2 - cosmicflash_chi2;
    // m_cosmicflash_ratio_dchi.push_back( dchi2 );

    if ( dchi2<0 ) {
      // matches cosmic flash better
      return false;
    }
    else {
      return true;
    }

    // // our cut needs to be position dependent, so we need a position. we get the q-weighted mean and the intervals.                                                                                   
    // float totw = 0.;
    // float mean[3] = {0};
    // float min_x = 1.0e6;
    // float max_x = -1.0e6;
    // for ( auto const& qpt : qcluster ) {
    //   mean[0] += qpt.x*qpt.q;
    //   mean[1] += qpt.y*qpt.q;
    //   mean[2] += qpt.z*qpt.q;
    //   totw += qpt.q;
    //   if ( qpt.x < min_x ) min_x = qpt.x;
    //   if ( qpt.x > max_x ) max_x = qpt.x;
    // }
    // if ( totw==0 ) {
    //   return false;
    // }
    // else {
    //   for (int i=0; i<3; i++) mean[i] /= totw;
    // }

    // //  still a dumb cut                                                                                                                                                                               
    // if ( mean[0] < 100.0 && min_chi2< m_config.flashmatch_chi2_cut*1.5 )
    //   return true;
    // else if ( mean[0]>100.0 && min_chi2<m_config.flashmatch_chi2_cut )
    //   return true;
    // else
    //   return false;
    

    return true; // never gets here
  }

  void TaggerFlashMatchAlgoV2::setupQClusters( const std::vector<TaggerFlashMatchData>& taggertracks_v ) {

    m_qcluster_v.clear();
    m_qcluster_extended_v.clear();
    
    // load qclusters, no extension for in-time tests
    for ( auto const& taggertrack : taggertracks_v ) {
      flashana::QCluster_t qcluster;
      m_genflashmatch.ExpandQClusterStartingWithLarliteTrack( qcluster, taggertrack.m_track3d, 0.0, false, false );
      m_qcluster_v.emplace_back( std::move(qcluster) );
    }
    
    // load qclusters, extension for out-of-time tests
    for ( auto const& taggertrack : taggertracks_v ) {
      
      // we extend if near the side boundarys
      //  or if there is a flash end
      
      auto const& tvecstart = taggertrack.m_track3d.Vertex();
      auto const& tvecend   = taggertrack.m_track3d.End();
      
      bool extend_start = false;
      if ( tvecstart.Y()<-107.0 || tvecstart.Y()>107.0 || tvecstart.Z()<10.0 || tvecstart.Z()>1026.0 )
	extend_start = true;
      
      bool extend_end   = false;      
      if ( tvecend.Y()<-107.0 || tvecend.Y()>107.0 || tvecend.Z()<10.0 || tvecend.Z()>1026.0 )
	extend_end = true;      
      
      flashana::QCluster_t qcluster;
      //m_genflashmatch.ExpandQClusterStartingWithLarliteTrack( qcluster, taggertrack.m_track3d, 1000.0, extend_start, extend_end );
      m_genflashmatch.ExpandQClusterStartingWithLarliteTrack( qcluster, taggertrack.m_track3d, 1000.0, extend_start, extend_end );      
      m_qcluster_extended_v.emplace_back( std::move(qcluster) );
    }
  }
  
  void TaggerFlashMatchAlgoV2::setupFlashMatchInterface( std::vector<flashana::Flash_t>& data_flashana,
						       std::vector<flashana::Flash_t>& cosmicdata_flashana,
						       bool use_extended ) {
    // we fill up the qcluster and flash objects into the flashmatchmanager located inside the
    //  the flash interface class: generalflashmatchalgo
    m_genflashmatch.getFlashMatchManager().Reset();
    
    if ( !use_extended ) {
      // in-time flash and non-extended tracks
      for ( auto& flash : data_flashana ) {
	m_genflashmatch.getFlashMatchManager().Add( flash, false );
      }
      for ( auto& qcluster : m_qcluster_v ) {
	m_genflashmatch.getFlashMatchManager().Add( qcluster );
      }
    }
    else {
      // cosmic-disc fash
      for ( auto& flash : cosmicdata_flashana ) {
	m_genflashmatch.getFlashMatchManager().Add( flash, true );
      }
      for ( auto& qcluster : m_qcluster_extended_v ) {
	m_genflashmatch.getFlashMatchManager().Add( qcluster );
      }
    }
  }

  std::vector<TaggerFlashMatchAlgoV2::FlashRange_t> TaggerFlashMatchAlgoV2::getFlashRange( const std::vector<flashana::Flash_t>& intime_flash_t ) {
    const larutil::Geometry* geo = larutil::Geometry::GetME();
    std::vector<FlashRange_t> flashrange_v(intime_flash_t.size());
    int iflash = -1;
    for ( auto const& flash : intime_flash_t ) {
      iflash++;
      FlashRange_t& range = flashrange_v[iflash];
      range.maxq  = 0;
      range.maxch = 0;
      range.meanz = 0;
      range.zfwhm[0] = 1050.0; // min
      range.zfwhm[1] = 0;      // max

      float totpe = 0;

      // fill maxq, maxch, meanz
      for ( size_t ifemch=0; ifemch<flash.pe_v.size(); ifemch++ ) {
	float pe     = flash.pe_v[ifemch];
	std::vector<double> pos(3);
	geo->GetOpChannelPosition( ifemch, pos );
	range.meanz += pos[2]*pe;
	totpe       += pe;
	if ( pe>range.maxq ) {
	  range.maxq  = pe;
	  range.maxch = ifemch;
	}
      }
      range.meanz /= totpe;

      // get range
      for ( size_t ifemch=0; ifemch<flash.pe_v.size(); ifemch++ ) {
	float pe     = flash.pe_v[ifemch];
	if ( pe < 0.5*range.maxq )
	  continue;
	std::vector<double> pos(3);
	geo->GetOpChannelPosition( ifemch, pos );

	if ( pos[2]<range.zfwhm[0] )
	  range.zfwhm[0] = pos[2];
	if ( pos[2]>range.zfwhm[1] )
	  range.zfwhm[1] = pos[2];
      }
      
    }

    return flashrange_v;
  }

  
}
