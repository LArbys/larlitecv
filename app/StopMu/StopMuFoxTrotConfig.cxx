#include "StopMuFoxTrotConfig.h"

namespace larlitecv {

  StopMuFoxTrotConfig::StopMuFoxTrotConfig() {
    min_num_steps = 3;
    // foxtrotalgo sets its own defaults
  }

  StopMuFoxTrotConfig StopMuFoxTrotConfig::makeFromPSet( const larcv::PSet& ps ) {
    StopMuFoxTrotConfig cfg;
    cfg.min_num_steps = ps.get<int>("MinNumSteps");
    cfg.foxtrotalgo_cfg = FoxTrotTrackerAlgoConfig::makeFromPSet( ps );
    return cfg;
  }
  
}
