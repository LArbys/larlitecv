#ifndef __TAGGER_FLASH_MATCH_ALGO_H__
#define __TAGGER_FLASH_MATCH_ALGO_H__

// LArCV
#include "Base/PSet.h"

// larlite
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h" // SelectionTool
#include "OpT0Finder/Base/OpT0FinderTypes.h"
#include "FhiclLite/PSet.h"

// larlitecv
#include "TaggerFlashMatchAlgoConfig.h"


namespace larlitecv {

  class TaggerFlashMatchAlgo {

  	TaggerFlashMatchAlgo() {};
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


  	flashana::QCluster_t GenerateQCluster( const TaggerFlashMatchData& data );
	  flashana::Flash_t    GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster );
	  flashana::Flash_t    MakeDataFlash( const larlite::opflash& flash );
	  void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );
	  std::vector<larlite::opflash>     SelectInTimeOpFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );
    std::vector<flashana::Flash_t>    GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v );
	  std::vector<flashana::QCluster_t> GenerateQClusters( const TaggerFlashMatchData& inputdata );
	  void ChooseContainedCandidates( const TaggerFlashMatchData& inputdata, std::vector<int>& passes_containment );
	  bool IsClusterContained( const TaggerFlashMatchData& data );

  protected:

  	const TaggerFlashMatchAlgoConfig m_config;   //< configuration class
		flashana::FlashMatchManager m_flash_matcher; //< Generates Flash Hypothesis producer


  };

}


#endif
