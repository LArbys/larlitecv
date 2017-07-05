/* Function to compare a flash of light with a track using the flash matching infrastructure. 

   This is meant to introduce flash matching into the tagger's operations. */
#ifndef __GENERAL_FLASHMATCH_ALGO_H__
#define __GENERAL_FLASHMATCH_ALGO_H__

#include <vector>
#include <cmath>

// larlite 
#include "OpT0Finder/Algorithms/QLLMatch.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h"

// larcv
#include "Base/PSet.h"

// larlitecv
#include "TaggerCROI/TaggerCROITypes.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "GeneralFlashMatchAlgoConfig.h"

namespace larlitecv {

  // define the class
  class GeneralFlashMatchAlgo {

        GeneralFlashMatchAlgo();
	
  public:
        GeneralFlashMatchAlgo( GeneralFlashMatchAlgoConfig& config );
        virtual ~GeneralFlashMatchAlgo() {};

	// This might need more information based on the needs of converting a flash to an opflash.
	std::vector < int > GeneralFlashMatchAlgo::non_anodecathode_flash_idx( const std::vector < larlite::opflash* >& opflash_v, std::vector < int > anode_flash_idx_v, std::vector < int > cathode_flash_idx_v );
	std::vector < opflash* > generate_single_denomination_flash_list(const std::vector< larlite::opflash > full_opflash_v, std::vector < int > idx_v);
	flashana::QCluster_t GenerateQCluster( const larlite::track& track );
	void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );
	flashana::Flash_t GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster );
	larlite::opflash MakeOpFlashFromFlash(const flashana::Flash_t& Flash);
	flashana::Flash_t MakeDataFlash( const larlite::opflash& flash );
	larlite::opflash make_flash_hypothesis(const larlite::track input_track);
	std::vector<larlite::opflash> GeneralFlashMatchAlgo::make_flash_hypothesis_vector( const std::vector<larlite::opflash>& input_track_v );
	float generate_chi2_in_track_flash_comparison(const larlite::opflash opflash_hypo, const larlite::opflash data_opflash);

	const TaggerFlashMatchAlgoConfig& getConfig() { return m_config; };

  private:
	const TaggerFlashMatchAlgoConfig m_config;   //< configuration class
	larcv::pmtweights::PMTWireWeights m_pmtweights;
	flashana::FlashMatchManager m_flash_matcher; //< Generates Flash Hypothesis producer
	//SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
	int m_verbosity;
	std::map<int,int> m_opch_from_opdet;

  };


} 

#endif

