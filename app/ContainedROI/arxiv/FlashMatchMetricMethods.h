#ifndef __FLASH_MATCH_METRIC_METHODS_H__
#define __FLASH_MATCH_METRIC_METHODS_H__

#include <vector>

// LArLite
#include "DataFormat/opflash.h"
#include "DataFormat/track.h"

#include "TaggerFlashMatchAlgo.h"

namespace larlitecv {

  float CalculateFlashMatchChi2( const std::vector<larlite::opflash>& flash_data_v, const larlite::opflash& flash_hypothesis,
				 float& totpe_data, float& totpe_hypo, const float fudge_factor, const bool use_gaus2d, const bool verbose=false );

  float CalculateShapeOnlyUnbinnedLL( const std::vector<larlite::opflash>& flash_data_v, const larlite::opflash& flash_hypothesis,
    float& totpe_data, float& totpe_hypo, const bool verbose );

  float ScanFlashMatchChi2( const float dz, const float dx, const float dy, const float stepsize,
			    const std::vector<larlite::opflash>& flash_data_v, const larlite::track& taggertrack,
			    TaggerFlashMatchAlgo& algo );
  
  larlite::opflash GetGaus2DPrediction( const larlite::opflash& flash );  
  

}

#endif
