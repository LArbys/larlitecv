#ifndef __FOX_TROT_TRACKER_ALGO_CONFIG_H__
#define __FOX_TROT_TRACKER_ALGO_CONFIG_H__

#include "Base/PSet.h"

namespace larlitecv {

  class FoxTrotTrackerAlgoConfig {

  public:
    FoxTrotTrackerAlgoConfig();
    virtual ~FoxTrotTrackerAlgoConfig() {};

    float step_size;
    int num_step_attempts;
    std::vector<float> pixel_thresholds;
    int min_hit_width;
    int hit_neighborhood;
    float segment_frac_w_charge;
    float radius_reduction_factor;
    float min_cosine;
    int max_steps;
    int verbosity;

    static FoxTrotTrackerAlgoConfig makeFromPSet( const larcv::PSet& pset );

  };

}


#endif

