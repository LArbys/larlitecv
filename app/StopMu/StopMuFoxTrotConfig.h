#ifndef __STOPMU_FOXTROT_CONFIG_H__
#define __STOPMU_FOXTROT_CONFIG_H__

// larcv
#include "Base/PSet.h"

#include "FoxTrotTrackerAlgoConfig.h"

namespace larlitecv {

  class StopMuFoxTrotConfig {
  public:
    StopMuFoxTrotConfig();
    virtual ~StopMuFoxTrotConfig() {};

    FoxTrotTrackerAlgoConfig foxtrotalgo_cfg;
    int min_num_steps;

    static StopMuFoxTrotConfig makeFromPSet( const larcv::PSet& ps );
    
  };

}


#endif