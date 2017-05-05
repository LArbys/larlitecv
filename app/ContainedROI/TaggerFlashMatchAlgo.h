#ifndef __TAGGER_FLASH_MATCH_ALGO_H__
#define __TAGGER_FLASH_MATCH_ALGO_H__

#include <map>

// LArCV
#include "Base/PSet.h"
#include "DataFormat/ROI.h"
#include "PMTWeights/PMTWireWeights.h"

// larlite
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h" // SelectionTool
#include "OpT0Finder/Base/OpT0FinderTypes.h"
#include "FhiclLite/PSet.h"

// larlitecv
#include "TaggerFlashMatchAlgoConfig.h"
#include "TaggerFlashMatchTypes.h"
//#include "SCE/SpaceChargeMicroBooNE.h"


namespace larlitecv {

  class TaggerFlashMatchAlgo {

  	TaggerFlashMatchAlgo();
  public:
  	TaggerFlashMatchAlgo( TaggerFlashMatchAlgoConfig& config );
  	virtual ~TaggerFlashMatchAlgo() {};

    // collect the clusters. what are their representation?
    // I have vector<Pixel2D> and larlite::spacepoints
    // for both thrumu, stopmu, and contained clusters

    // we probably should do some spatial-merging?

    // we use all our 3D info to build a flash hypothesis

    // we make decision based on
    // (1) flash-hypothesis matching, this is a loose cut chi-2 cut.
    // (2) containment, where we are being cautious. containment includes out-of-tpc consideration

   std::vector<larcv::ROI> FindFlashMatchedContainedROIs( const std::vector<TaggerFlashMatchData>& inputdata,
		       const std::vector<larlite::event_opflash*>& opflashes_v, std::vector<int>& flashdata_selected );

   void setVerbosity( int verbosity ) { m_verbosity = verbosity; };

   flashana::QCluster_t GenerateQCluster( const TaggerFlashMatchData& data );
   flashana::QCluster_t GenerateQCluster( const larlite::track& track );
   flashana::Flash_t    GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster );
   flashana::Flash_t    MakeDataFlash( const larlite::opflash& flash );
   void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );
   std::vector<flashana::Flash_t>    GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v );
   std::vector<larlite::opflash>     SelectInTimeOpFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );
   std::vector<flashana::Flash_t>    MakeDataFlashes( const std::vector<larlite::opflash>& opflashes );
   std::vector<flashana::QCluster_t> GenerateQClusters( const std::vector<TaggerFlashMatchData>& inputdata );
   void ChooseContainedCandidates( const std::vector<TaggerFlashMatchData>& inputdata, std::vector<int>& passes_containment );
   bool IsClusterContained( const TaggerFlashMatchData& data );
   float InTimeFlashComparison( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster, float& totpe_data, float& totpe_hypo );

   larlite::opflash MakeOpFlashFromFlash(const flashana::Flash_t& inFlash);
   void clearFlashInfo() {
     m_opflash_hypos.clear();
     m_min_chi2.clear();
   };
   const std::vector<larlite::opflash>& getOpFlashHypotheses() { return m_opflash_hypos; };

    static std::vector< std::vector<float> > GetAABoundingBox( const larlite::track& track );

    bool DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster, float& totpe_data, float& totpe_hypo );

    bool DoesQClusterMatchInTimeBetterThanCosmic( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster,
						  const larlitecv::TaggerFlashMatchData& taggertrack, const int trackidx, float& dchi2 );
    bool DoesTotalPEMatch( float totpe_data, float totpe_hypo );
    
    bool didTrackPassContainmentCut( int itrack );
    bool didTrackPassFlashMatchCut( int itrack );
    bool didTrackPassCosmicFlashCut( int itrack );

    const std::vector<int>& getContainmentCutResults() { return m_passes_containment; };
    const std::vector<int>& getFlashMatchCutResults() { return m_passes_flashmatch; };
    const std::vector<int>& getCosmicRatioCutResults() { return m_passes_cosmicflash_ratio; };
    const std::vector<int>& getTotalPECutResults() { return m_passes_totpe; };    
    
    const TaggerFlashMatchAlgoConfig& getConfig() { return m_config; };

  protected:

    const TaggerFlashMatchAlgoConfig m_config;   //< configuration class
    larcv::pmtweights::PMTWireWeights m_pmtweights;
    flashana::FlashMatchManager m_flash_matcher; //< Generates Flash Hypothesis producer
    //SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
    int m_verbosity;
    std::map<int,int> m_opch_from_opdet;
    std::vector<larlite::opflash> m_opflash_hypos; ///< flash hypotheses for each track
    std::vector<float> m_min_chi2; ///< min-chi2 score for each flash

    std::vector<int> m_passes_containment;
    std::vector<int> m_passes_flashmatch;
    std::vector<int> m_passes_totpe;    
    std::vector<int> m_passes_cosmicflash_ratio;
    

  };

}


#endif
