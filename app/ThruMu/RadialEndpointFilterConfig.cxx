#include "RadialEndpointFilterConfig.h"

namespace larlitecv {

  RadialEndpointFilterConfig::RadialEndpointFilterConfig() {
    // default constructor with some sensible defaults
    segment_radius = 10.0; // cm
    segment_min_width = 1; // pixels
    segment_frac_w_charge = 0.8; // fraction
    acceptance_angle = (5.0/180.0)*3.14159; // 5 degrees in radians
    pixel_thresholds.resize(3,10.0);
  }

}
