#ifndef __FOXTROTLEAD_H__
#define __FOXTROTLEAD_H__

/* =======================================================
 * FoxTrotLead
 * 
 * abstract interface class for a function that chooses the
 * next step for the FoxTrotAlgo.
 *
 * =======================================================*/

#include <vector>
#include "Segment3DAlgoTypes.h"
#include "FoxTrotTrackerAlgoTypes.h"
#include "FoxTrotTrackerAlgoConfig.h"

namespace larlitecv {

  class FoxTrotLead {
  public:
    FoxTrotLead() {};
    virtual ~FoxTrotLead() {};

    // implement this
    virtual bool chooseBestSegment( const FoxStep& current, const std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past,
				    const FoxTrotTrackerAlgoConfig& config, int& best_seg_idx ) = 0;
    
  };

}

#endif
