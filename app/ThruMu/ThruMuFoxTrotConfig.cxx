#include "ThruMuFoxTrotConfig.h"

namespace larlitecv {

  ThruMuFoxTrotConfig::ThruMuFoxTrotConfig() {
    endpoint_radius = 5.0;
    maxiters = 10;
    hit_neighborhood = 1;
    use_thrumu_lead = true;
    pixel_thresholds.resize(3,10.0);
  }

}
