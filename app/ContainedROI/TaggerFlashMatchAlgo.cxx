#include "TaggerFlashMatchAlgo.h"

#include <sstream>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"

namespace larlitecv {

  TaggerFlashMatchAlgo::TaggerFlashMatchAlgo() {
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << " Default constructor only declared for ROOT dict. purposes. Should not be used." << std::endl;
    throw std::runtime_error(msg.str());
  }

  TaggerFlashMatchAlgo::TaggerFlashMatchAlgo( TaggerFlashMatchAlgoConfig& config )
    : m_config(config),m_genflashmatch(config.genflashmatch_cfg) {
    // set verbosity
    setVerbosity( m_config.verbosity );
    // opdet-opch map
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

    m_min_chi2.clear();
    m_opflash_hypos.clear();
    m_totpe_peratio.clear();
    m_cosmicflash_ratio_dchi.clear();
    m_min_chi2.reserve( inputdata.size() );
    m_totpe_peratio.reserve( inputdata.size() );
    m_cosmicflash_ratio_dchi.reserve( inputdata.size() );
    m_genflashmatch.getFlashMatchManager().Reset();    

    // we fill the flash-match code with qclusters--with extensions along with the flashes
    setupFlashMatchInterface( data_flashana, cosmicdata_flashana, inputdata );

    // we do the many-to-many match: will use this later when comparing in-time to cosmic flashes
    m_genflashmatch.getFlashMatchManager().Match();

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

      larlite::opflash ophypo = m_genflashmatch.make_flash_hypothesis( inputdata[i].m_track3d );
      m_opflash_hypos.emplace_back( std::move(ophypo) );
      
      // in-time flash matching
      float totpe_data = 0;
      float totpe_hypo = 0;
      float min_chi2   = -1;
      m_passes_flashmatch[i] = ( DoesQClusterMatchInTimeFlash( data_flashana, inputdata[i], totpe_data, totpe_hypo, min_chi2 )  ) ? 1 : 0;

      // totpe
      m_passes_totpe[i] = ( m_genflashmatch.DoesTotalPEMatch( totpe_data, totpe_hypo ) ) ? 1 : 0;

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
    // we use the generalflashmatch algo qcluster maker, which will extend the ends if needed.
    flashana::QCluster_t qclust;
    const larlite::track& track = data.m_track3d;
    m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qclust, track, 1000.0, 10.0 );
    return qclust;
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

  bool TaggerFlashMatchAlgo::DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v,
							   const larlitecv::TaggerFlashMatchData& taggertrack,
							   float& totpe_data, float& totpe_hypo, float& min_chi2 ) {
    flashana::QCluster_t qcluster;
    // generate qcluster using genflash. Do not extend as this is the in-time hypothesis
    m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qcluster, taggertrack.m_track3d, 1000.0, 10.0 );
    
    bool matches = m_genflashmatch.DoesQClusterMatchInTimeFlash( intime_flashes_v, qcluster, totpe_data, totpe_hypo, min_chi2 );
    
    return matches;
  }
  
  bool TaggerFlashMatchAlgo::DoesQClusterMatchInTimeBetterThanCosmic( const std::vector<flashana::Flash_t>& cosmic_flashes_v,
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

  void TaggerFlashMatchAlgo::setupFlashMatchInterface( std::vector<flashana::Flash_t>& data_flashana,
						       std::vector<flashana::Flash_t>& cosmicdata_flashana,
						       const std::vector<TaggerFlashMatchData>& taggertracks_v ) {
    // we fill up the qcluster and flash objects into the flashmatchmanager located inside the
    //  the flash interface class: generalflashmatchalgo

    // in-time flash
    for ( auto& flash : data_flashana ) {
      m_genflashmatch.getFlashMatchManager().Add( flash, false );
    }

    // cosmic-disc fash
    for ( auto& flash : cosmicdata_flashana ) {
      m_genflashmatch.getFlashMatchManager().Add( flash, true );
    }

    // load qclusters
    for ( auto const& taggertrack : taggertracks_v ) {
      flashana::QCluster_t qcluster;
      m_genflashmatch.ExpandQClusterNearBoundaryFromLarliteTrack( qcluster, taggertrack.m_track3d, 10000.0, 15.0 );
      m_genflashmatch.getFlashMatchManager().Emplace( std::move(qcluster) );
    }
    
  }

}
