#ifndef __GENERAL_FLASHMATCH_ALGO_H__
#define __GENERAL_FLASHMATCH_ALGO_H__

/* ============================================================================
 * GeneralFlashMatchAlgo
 * Methods to compare a flash of light with a track using the larlite flash matching infrastructure. 
 *
 * This is meant to introduce flash matching into the tagger's operations. 
 *
 * author(s): Chris Barnes (barnchri@umich.edu)
 *
 * ============================================================================*/


#include <vector>
#include <cmath>

// larlite 
#include "OpT0Finder/Algorithms/QLLMatch.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h"
#include "OpT0Finder/Algorithms/LightPath.h"
#include "DataFormat/mctrack.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/TimeService.h" // For the G4 to electronics time.

// larcv
#include "Base/PSet.h"
#include "DataFormat/ROI.h"
#include "PMTWeights/PMTWireWeights.h"

// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "GeneralFlashMatchAlgoConfig.h"

namespace larlitecv {

  // define the class
  class GeneralFlashMatchAlgo {
  public:
    GeneralFlashMatchAlgo();
    GeneralFlashMatchAlgo( const GeneralFlashMatchAlgoConfig& config );
    virtual ~GeneralFlashMatchAlgo() {};

    // This might need more information based on the needs of converting a flash to an opflash.
    std::vector < int > non_anodecathode_flash_idx( const std::vector < larlite::opflash* >& opflash_v,
						    std::vector < int > anode_flash_idx_v, std::vector < int > cathode_flash_idx_v );

    std::vector < larlite::opflash* > generate_single_denomination_flash_list( std::vector< larlite::opflash > full_opflash_v, std::vector < int > idx_v);
    
    std::vector<larlite::opflash> SelectInTimeOpFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );

    bool DoesTotalPEMatch( float totpe_data, float totpe_hypo );
    
    std::vector<flashana::Flash_t> GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v );

    std::vector<flashana::Flash_t> GetCosmicFlashana( const std::vector<larlite::event_opflash*>& opflashes_v );

    bool DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster, float& totpe_data, float& totpe_hypo, float& min_chi2 );


    void GetFlashMinChi2( const std::vector<flashana::Flash_t>& flashana_v, const flashana::QCluster_t& qcluster, const int in_or_out_of_time,
			  float& totpe_data, float& totpe_hypo, float& min_chi2, int& flash_idx );

    std::vector< BoundaryFlashIndex > generate_boundaryflashindex_vector_for_event( std::vector< larlite::event_opflash* > opflashes_v );

    std::vector< larlite::opflash > generate_single_opflash_vector_for_event( std::vector< larlite::event_opflash* > opflashes_v );

    std::vector< int > generate_single_opflash_idx_vector_for_event( std::vector< larlite::event_opflash* > opflashes_v );

    std::vector< larlite::track >  generate_tracks_between_passes( const std::vector< larlitecv::BMTrackCluster3D > trackcluster3d_v);
    
    void ExpandQClusterStartingWithLarliteTrack( flashana::QCluster_t& qcluster, const larlite::track& larlite_track,
						 double extension, bool extend_start, bool extend_end );

    void ExpandQClusterNearBoundaryFromLarliteTrack( flashana::QCluster_t& qcluster, const larlite::track& larlite_track,
						     const double extension, const double boundarydist );

    bool isNearActiveVolumeEdge( ::larlite::geoalgo::Vector pt, double d );
    
    void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );

    flashana::Flash_t GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster );

    larlite::opflash MakeOpFlashFromFlash( const flashana::Flash_t& Flash );

    larlite::opflash GetGaus2DPrediction( const larlite::opflash& flash );
    
    flashana::Flash_t MakeDataFlash( const larlite::opflash& flash );

    std::vector< flashana::Flash_t > MakeDataFlashes( std::vector< larlite::opflash > opflash_v );

    larlite::opflash make_flash_hypothesis(const larlite::track input_track);

    std::vector<larlite::opflash> make_flash_hypothesis_vector( const std::vector<larlite::track>& input_track_v );
    
    float generate_chi2_in_track_flash_comparison(const flashana::QCluster_t qcluster, const larlite::opflash data_opflash, float& totpe_data_flash, float& totpe_hypo_flash, int flash_prod_idx);

    const GeneralFlashMatchAlgoConfig& getConfig() { return m_config; };

    void setVerbosity( int v ) { m_verbosity = v; };

    static GeneralFlashMatchAlgo* GetME();
    static GeneralFlashMatchAlgo* GetME( const larlitecv::GeneralFlashMatchAlgoConfig& config );

    flashana::FlashMatchManager& getFlashMatchManager() { return m_flash_matcher; }; //< convenience function

  protected:

    // CONTAINMENT CUT                                                                                                                                                                              
    std::vector<int> m_passes_containment;
    std::vector<float> m_containment_dwall;

    // IN-TIME FLASH CUT                                                                                                                                                                              
    std::vector<int> m_passes_flashmatch;
    std::vector<larlite::opflash> m_opflash_hypos; ///< flash hypotheses for each track                                                                                                               
    std::vector<float> m_min_chi2; ///< min-chi2 score for each flash                                                                                                                                 

    // TOTAL PE                                                                                                                                                                                     
    std::vector<int> m_passes_totpe;
    std::vector<float> m_totpe_peratio; //< fraction difference of total pe between hypothesis and in-time                                                                                             

    std::vector<int> m_passes_cosmicflash_ratio;
    std::vector<float> m_cosmicflash_ratio_dchi; //< delta-Chi-squared between in-time flash and matched flash to track                                                                                 


  private:
    const GeneralFlashMatchAlgoConfig m_config;   //< configuration class
    larcv::pmtweights::PMTWireWeights m_pmtweights;
    flashana::FlashMatchManager m_flash_matcher; //< Generates Flash Hypothesis producer
    //SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
    int m_verbosity;
    std::map<int,int> m_opch_from_opdet;

    static GeneralFlashMatchAlgo* _global_instance;

  };


} 

#endif

