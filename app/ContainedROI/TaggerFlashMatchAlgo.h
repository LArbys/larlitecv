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
#include "SCE/SpaceChargeMicroBooNE.h"


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
	    const std::vector<larlite::event_opflash*>& opflashes_v );

	  void setVerbosity( int verbosity ) { m_verbosity = verbosity; };


  	flashana::QCluster_t GenerateQCluster( const TaggerFlashMatchData& data );
	  flashana::Flash_t    GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster );
	  flashana::Flash_t    MakeDataFlash( const larlite::opflash& flash );
	  void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );
    std::vector<flashana::Flash_t>    GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v );	  
	  std::vector<larlite::opflash>     SelectInTimeOpFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );
	  std::vector<flashana::Flash_t>    MakeDataFlashes( const std::vector<larlite::opflash>& opflashes );
	  std::vector<flashana::QCluster_t> GenerateQClusters( const std::vector<TaggerFlashMatchData>& inputdata );	  
	  void ChooseContainedCandidates( const std::vector<TaggerFlashMatchData>& inputdata, std::vector<int>& passes_containment );
	  bool IsClusterContained( const TaggerFlashMatchData& data );
    float InTimeFlashComparison( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster );
    void ChooseInTimeFlashMatchedCandidates( const std::vector<flashana::QCluster_t>& inputdata, 
    	const std::vector<flashana::Flash_t>& intime_flashes, std::vector<int>& passes_flashmatch );


  protected:

  	const TaggerFlashMatchAlgoConfig m_config;   //< configuration class
    larcv::pmtweights::PMTWireWeights m_pmtweights;  	
		flashana::FlashMatchManager m_flash_matcher; //< Generates Flash Hypothesis producer
		SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
		int m_verbosity;
		std::map<int,int> m_opch_from_opdet;



  };

}


#endif
