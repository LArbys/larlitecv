#ifndef __MATCH_TAGGER_DATA_2_FLASH_H__
#define __MATCH_TAGGER_DATA_2_FLASH_H__

// std lib
#include <vector>

// larlite
#include "DataFormat/opflash.h"

// larlitecv
#include "TaggerFlashMatchTypes.h"
#include "TaggerTypes/BoundarySpacePoint.h"

namespace larlitecv {

  void MatchTaggerData2Flash( std::vector< TaggerFlashMatchData >& taggerdata_v, const std::vector< larlite::event_opflash* >& opflash_v,
			      const std::vector< BoundarySpacePoint >& anode_spacepts, const std::vector< BoundarySpacePoint >& cathode_spacepts, const float max_dist );

}

#endif
