#ifndef __LINEAR3D_CHARGE_TAGGER_CONFIG__
#define __LINEAR3D_CHARGE_TAGGER_CONFIG__

// larcv
#include "Base/PSet.h"

namespace larlitecv {

  class Linear3DChargeTaggerConfig {
  public:
    Linear3DChargeTaggerConfig();
    virtual ~Linear3DChargeTaggerConfig() {};

    float trigger_tpc_tick;
    float min_ADC_value;
    float step_size;
    int neighborhood_square;
    int neighborhood_posttick;// allows us to extend in later ticks to account for space-charge effect delay

    static Linear3DChargeTaggerConfig makeFromPSet( const larcv::PSet& pset );

  };

}

#endif
