#include "Linear3DChargeTaggerConfig.h"

namespace larlitecv {

  // ================================================================
  // Configuration Class of Linear3DChargeTagger

  Linear3DChargeTaggerConfig::Linear3DChargeTaggerConfig() {
    trigger_tpc_tick = 3200.0;
    min_ADC_value = 10.0; // works for 6 tick downsampling, no wire downsampling
    step_size = 0.3; // cm
    neighborhood_square = 5;
    neighborhood_posttick = 10;
  }

  Linear3DChargeTaggerConfig Linear3DChargeTaggerConfig::makeFromPSet( const larcv::PSet& pset ) {
    Linear3DChargeTaggerConfig cfg;
    cfg.trigger_tpc_tick      = pset.get<float>( "TriggerTPCTick" );
    cfg.min_ADC_value         = pset.get<float>( "PixelThreshold");
    cfg.step_size             = pset.get<float>( "StepSize" );
    cfg.neighborhood_square   = pset.get<int>( "NeighborhoodSquareSize");
    cfg.neighborhood_posttick = pset.get<int>( "NeighborhoodPostTick" );
    return cfg;
  }

  // ================================================================
  
  
}
