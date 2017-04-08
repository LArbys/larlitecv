#ifndef __FLASH_MATCH_METRIC_METHODS_H__
#define __FLASH_MATCH_METRIC_METHODS_H__

#include <vector>

// LArLite
#include "DataFormat/opflash.h"

namespace larlitecv {

  float CalculateFlashMatchChi2( const std::vector<larlite::opflash>& flash_data_v, const larlite::opflash& flash_hypothesis,
    float& totpe_data, float& totpe_hypo, const float fudge_factor, const bool verbose=false );

  float CalculateShapeOnlyUnbinnedLL( const std::vector<larlite::opflash>& flash_data_v, const larlite::opflash& flash_hypothesis,
    float& totpe_data, float& totpe_hypo, const bool verbose );

}

#endif
