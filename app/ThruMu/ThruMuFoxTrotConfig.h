#ifndef __THRUMU_FOX_TROT_CONFIG_H__
#define __THRUMU_FOX_TROT_CONFIG_H__

#include "ChargeSegmentAlgos/FoxTrotTrackerAlgoConfig.h"

namespace larlitecv {

  class ThruMuFoxTrotConfig {
  public:
    ThruMuFoxTrotConfig();
    virtual ~ThruMuFoxTrotConfig() {};

    float endpoint_radius;
    int maxiters;
    int hit_neighborhood;
    bool use_thrumu_lead;
    std::vector<float> pixel_thresholds;

    FoxTrotTrackerAlgoConfig foxtrotalgo_cfg;

  };
}
#endif
