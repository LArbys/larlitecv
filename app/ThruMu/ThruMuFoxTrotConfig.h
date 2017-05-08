#ifndef __THRUMU_FOX_TROT_CONFIG_H__
#define __THRUMU_FOX_TROT_CONFIG_H__

#include "ChargeSegmentAlgos/FoxTrotTrackerAlgoConfig.h"

namespace larlitecv {

  class ThruMuFoxTrotConfig {
  public:
    ThruMuFoxTrotConfig();
    virtual ~ThruMuFoxTrotConfig() {};

    
    FoxTrotTrackerAlgoConfig foxtrotalgo_cfg;

  };
}
#endif
