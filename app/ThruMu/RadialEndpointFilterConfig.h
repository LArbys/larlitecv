#ifndef __RADIAL_ENDPOINT_FILTER_CONFIG_H__
#define __RADIAL_ENDPOINT_FILTER_CONFIG_H__

/* ==================================================================
 *  This class holds parameters for running the radialendpointfilter 
 *  class. the algo doesn't require a config instance, but we provide
 *  one in order to provide a place where some sane default values
 *  are listed.
 *
 * author(s):
 *  taritree wongjirad taritree@mit.edu

 * rev:
 *  2017/04/24: initial draft
 * 
 * =================================================================== */

#include <vector>

namespace larlitecv {
  
  class RadialEndpointFilterConfig {
  public:
    RadialEndpointFilterConfig();
    virtual ~RadialEndpointFilterConfig() {};

    float segment_radius;    // half-edge of box around 3D point that we look for tracks to cross, defining a segment
    float segment_min_width; // min width in pixels that a clump of charge on the box boundary needs to be
    float segment_frac_w_charge; // if we draw a line from the box center point to a crossing charge cluster, this is the fraction of pixels along line that must have charge
    float acceptance_angle; // if two segments found, we check if the angle between is PI +/- acceptance_angle in radians
    std::vector<float> pixel_thresholds;

  };
  
}

#endif
