#ifndef __RADIAL_ENDPOINT_FILTER_H__
#define __RADIAL_ENDPOINT_FILTER_H__

/* ==================================================================
 *  This class filters out boundary space points
 *  based on if the end points are on the end of
 *  a track or not.  Uses ChargeSegmentAlgos/RadialSegmentSearchAlgo
 *  to check around a 3D point. The algo defines a 3D cube
 *  around the point, and if we can find two segments from the
 *  center point to the edges that are 180 deg. of one another, then
 *  the point is considered to be along a track and rejected.
 *
 * author(s):
 *  taritree wongjirad taritree@mit.edu

 * rev:
 *  2017/04/24: initial draft
 * 
 * =================================================================== */

#include <vector>

#include "DataFormat/Image2D.h"

#include "RadialEndpointFilterConfig.h"

namespace larlitecv {

  class RadialEndpointFilter {
  public:
    RadialEndpointFilter() {};
    virtual ~RadialEndpointFilter() {};

    bool isWithinStraightSegment( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
				  const float segment_radius, const int segment_min_width, const int hit_neighborhood, const float segment_frac_w_charge,
				  const std::vector<float>& pixel_thresholds, const float acceptance_angle, int& num_seg3d );

    bool isWithinStraightSegment( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
				  const RadialEndpointFilterConfig& config, int& num_seg3d );
    
    
  };

}

#endif
