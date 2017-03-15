#ifndef __TAGGER_FLASH_MATCH_ALGO_CONFIG_H__
#define __TAGGER_FLASH_MATCH_ALGO_CONFIG_H__

#include <vector>

// LArCV
#include "Base/PSet.h"

// LArLite
#include "FhiclLite/PSet.h"

namespace larlitecv {

	class TaggerFlashMatchAlgoConfig {
  public:
  	TaggerFlashMatchAlgoConfig();
  	virtual ~TaggerFlashMatchAlgoConfig() {};

  	float qcluster_stepsize;
  	float MeV_per_cm;
  	float fudge_factor;
  	std::vector<int> beam_tick_range;
  	larlite::PSet flashmatch_config;
  	std::vector<float> FVCutX;
  	TGraph FVCutCurveY;
  	TGraph FVCutCurveZ;

  	static TaggerFlashMatchAlgoConfig MakeTaggerFlashMatchAlgoConfigFromPSet( const larcv::PSet& pset );
	};

}

#endif