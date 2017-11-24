#ifndef __TAGGER_FLASH_MATCH_ALGO_V2_H__
#define __TAGGER_FLASH_MATCH_ALGO_V2_H__

/* ============================================================================
 * TaggerFlashMatchAlgoV2
 * 
 * Takes tagger-track objects and uses flash matching to determine ROIs
 *  around the tracks for downstream analysis.
 *
 *  This version makes a selection only on position of flash and track.
 *  Also filters out tracks that match to cosmic discriminator flashes.
 *  Version made in response to BNB vs. EXTBNB disagreement.
 * 
 *
 * This is meant to introduce flash matching into the tagger's operations. 
 *
 * author(s): Chris Barnes (barnchri@umich.edu), Taritree (twongj01@tufts.edu)
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
#include "MCTruthTools/crossingPointsAnaMethods.h"


namespace larlitecv {

  class TaggerFlashMatchAlgoV2 {

    TaggerFlashMatchAlgoV2(); //< default declared for ROOT dict. Do not use.

  public:
    TaggerFlashMatchAlgoV2( TaggerFlashMatchAlgoConfig& config ); // Use this constructor.
    virtual ~TaggerFlashMatchAlgoV2() {};


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
    
   
    /* bool didTrackPassContainmentCut( int itrack ); */
    /* bool didTrackPassFlashMatchCut( int itrack ); */
    /* bool didTrackPassCosmicFlashCut( int itrack ); */

    flashana::QCluster_t GenerateQCluster( const TaggerFlashMatchData& data );

    struct FlashRange_t {
      int  maxch;
      float maxq;
      float meanz;
      float zfwhm[2];
    };
    std::vector<FlashRange_t> getFlashRange( const std::vector<flashana::Flash_t>& intime_flash_t );

    void matchFlashAndTrackPosition( const std::vector<FlashRange_t>& rangeinfo_v, const std::vector<TaggerFlashMatchData>& inputdata,
				     std::vector<int>& passes_intimepos, std::vector<double>& trackend_zdiff_frac );
   
    // selection results
    const std::vector<int>& getIntimeFlashPosResults()   const { return m_passes_intimepos;  };    
    const std::vector<int>& getContainmentCutResults()   const { return m_passes_containment; };
    const std::vector<int>& getCosmicFlashMatchResults() const { return m_passes_cosmic_flashmatch; };
    const std::vector<int>& getFinalResults()            const { return m_passes_finalresult; };
   
    // selection variables
    const std::vector<double>& getIntimeMeanZ()          const { return m_intime_meanz; };
    const std::vector<double>& getIntimeZFWHM()          const { return m_intime_zfwhm; };
    const std::vector<double>& getIntimePEMax()          const { return m_intime_pemax; };
    const std::vector<double>& getIntimeZDiffFrac()      const { return m_trackend_zdiff_frac; };    
    const std::vector<double>& getContainmentCutValues() const { return m_containment_dwall; };
    const std::vector<double>& getCosmicFlashMatchChi2Values()  const { return m_cosmic_bestflash_chi2_v; };
    const std::vector<int>&    getCosmicMatchedFlashIndices()   const { return m_cosmic_bestflash_idx_v; };
    const std::vector<int>&    getCosmicBestFlashMCTrackIDs()   const { return m_cosmic_bestflash_mctrackid_v; };    
    const std::vector<int>&    getCosmicMCTrackIDs()     const { return m_cosmic_mctrackid_v; };
    const int                  getMatchableFlashes()     const { return m_num_matchable_flashes; };
    const int                  getMatchedFlashes()       const { return m_num_matched_flashes; };

    //const std::vector<larlite::opflash>& getOpFlashHypotheses() const { return m_opflash_hypos; };
    
    const TaggerFlashMatchAlgoConfig& getConfig() { return m_config; };

    void provideTruthCrossingData( const CrossingPointAnaData_t& truthxingdata ) { ptruthxingdata=&truthxingdata; };
    
  protected:

    void setupQClusters( const std::vector<TaggerFlashMatchData>& taggertracks_v );
    
    void setupFlashMatchInterface( std::vector<flashana::Flash_t>& data_flashana,
				   std::vector<flashana::Flash_t>& cosmicdata_flashana,
				   bool use_extended );

    
    const TaggerFlashMatchAlgoConfig m_config;   //< configuration class
   
    GeneralFlashMatchAlgo m_genflashmatch;   //< interface to flash match tools
   
    //SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
    int m_verbosity;
    std::map<int,int> m_opch_from_opdet;

    // MC Truth variables -- to study flash matching accuracy
    const CrossingPointAnaData_t* ptruthxingdata;
    std::vector<int> m_cosmicflash_mctrackid;
   
    // -----------------------------------------------------------------------------------------
    // We save the cut results and cut variables of this algo. We'll save it to a larlite::user_info product if requested
    // This is useful for study of the cuts

    void clearResults( int n_intime_flashes, int n_tracks );

    // IN-TIME POSITION CUT
    std::vector<int>     m_passes_intimepos;
    std::vector<double>  m_intime_meanz;        ///< pe-weighted z position
    std::vector<double>  m_intime_zfwhm;        ///< z range that includes tubes above half max
    std::vector<double>  m_intime_pemax;        ///< pe max value
    std::vector<double>  m_trackend_zdiff_frac; ///< passes cut: (meanz-trackend)/zfwhm
    
    // CONTAINMENT CUT
    std::vector<int>    m_passes_containment;
    std::vector<double> m_containment_dwall;
   
    // COSMIC-FLASH MATCH
    std::vector<int>                  m_passes_cosmic_flashmatch;
    std::vector<double>               m_cosmic_bestflash_chi2_v; //< delta-Chi-squared between in-time flash and matched flash to track
    std::vector<flashana::QCluster_t> m_qcluster_v;              // w/o extension -- for in-time tests
    std::vector<flashana::QCluster_t> m_qcluster_extended_v;     // w/ extension  -- for cosmic-disc tests
    std::vector<flashana::Flash_t>    m_cosmic_bestflash_hypo_v;
    std::vector<int>                  m_cosmic_bestflash_idx_v;
    std::vector<int>                  m_cosmic_bestflash_mctrackid_v; //< mctrack id of paired flash
    std::vector<int>                  m_cosmic_mctrackid_v;           //< mctrack id of track (only if loaded MC track)
    int                               m_num_matchable_flashes;
    int                               m_num_matched_flashes;


    // FINAL RESULT
    std::vector<int> m_passes_finalresult;
   
  };
  
}


#endif
