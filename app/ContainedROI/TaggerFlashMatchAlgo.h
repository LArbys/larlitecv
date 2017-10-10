#ifndef __TAGGER_FLASH_MATCH_ALGO_H__
#define __TAGGER_FLASH_MATCH_ALGO_H__

/* ============================================================================
 * TaggerFlashMatchAlgo
 * 
 * Takes tagger-track objects and uses flash matching to determine ROIs
 *  around the tracks for downstream analysis.
 * 
 *
 * This is meant to introduce flash matching into the tagger's operations. 
 *
 * author(s): Chris Barnes (barnchri@umich.edu)
 *
 * ============================================================================*/


// stdlib
#include <map>

// LArCV
#include "Base/PSet.h"
#include "LArCV/build/include/DataFormat/ROI.h"
#include "PMTWeights/PMTWireWeights.h"

// larlite
#include "DataFormat/opflash.h"

// larlitecv
#include "FlashMatchInterface/GeneralFlashMatchAlgo.h"
#include "TaggerFlashMatchAlgoConfig.h"
#include "TaggerFlashMatchTypes.h"


namespace larlitecv {

  class TaggerFlashMatchAlgo {

    TaggerFlashMatchAlgo(); //< default declared for ROOT dict. Do not use.

  public:
    TaggerFlashMatchAlgo( TaggerFlashMatchAlgoConfig& config ); // Use this constructor.
    virtual ~TaggerFlashMatchAlgo() {};


    // Primary Method
    std::vector<larcv::ROI> FindFlashMatchedContainedROIs( const std::vector<TaggerFlashMatchData>& inputdata,
							   const std::vector<larlite::event_opflash*>& opflashes_v, std::vector<int>& flashdata_selected );
    
    void setVerbosity( int verbosity ) { m_verbosity = verbosity; };
   
    bool DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v,
				       const larlitecv::TaggerFlashMatchData& taggertrack,
				       float& totpe_data, float& totpe_hypo, float& min_chi2 );

    bool IsClusterContained( const TaggerFlashMatchData& data );   
   
    void ChooseContainedCandidates( const std::vector<TaggerFlashMatchData>& inputdata, std::vector<int>& passes_containment );
   
    bool DoesQClusterMatchInTimeBetterThanCosmic( const std::vector<flashana::Flash_t>& cosmic_flashes_v,
						  const larlitecv::TaggerFlashMatchData& taggertrack,
						  const float intime_min_chi2, float& dchi2 );
    
   
    bool didTrackPassContainmentCut( int itrack );
    bool didTrackPassFlashMatchCut( int itrack );
    bool didTrackPassCosmicFlashCut( int itrack );

    flashana::QCluster_t GenerateQCluster( const TaggerFlashMatchData& data );
   
    // selection results
    const std::vector<int>& getContainmentCutResults() { return m_passes_containment; };
    const std::vector<int>& getFlashMatchCutResults() { return m_passes_flashmatch; };
    const std::vector<int>& getCosmicRatioCutResults() { return m_passes_cosmicflash_ratio; };
    const std::vector<int>& getTotalPECutResults() { return m_passes_totpe; };
   
    // selection variables
    const std::vector<float>& getContainmentCutValues() { return m_containment_dwall; };
    const std::vector<float>& getInTimeChi2Values() { return m_min_chi2; };
    const std::vector<float>& getTotPEratioValues() { return m_totpe_peratio; };
    const std::vector<float>& getCosmicDeltaChi2Values() { return m_cosmicflash_ratio_dchi; };
   
    const TaggerFlashMatchAlgoConfig& getConfig() { return m_config; };

    const std::vector<larlite::opflash>& getOpFlashHypotheses() const { return m_opflash_hypos; };
   
  protected:


    void setupFlashMatchInterface( std::vector<flashana::Flash_t>& data_flashana,
				   std::vector<flashana::Flash_t>& cosmicdata_flashana,
				   const std::vector<TaggerFlashMatchData>& taggertracks_v );

    
    const TaggerFlashMatchAlgoConfig m_config;   //< configuration class
   
    GeneralFlashMatchAlgo m_genflashmatch;   //< interface to flash match tools
   
    //SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
    int m_verbosity;
    std::map<int,int> m_opch_from_opdet;
   
    // -----------------------------------------------------------------------------------------
    // We save the cut results and cut variables of this algo. We'll save it to a larlite::user_info product if requested
    // This is useful for study of the cuts
   
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
   
  };
  
}


#endif
