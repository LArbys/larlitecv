#include "RadialEndpointFilter.h"

// larlitecv
#include "ChargeSegmentAlgos/Segment3DAlgoTypes.h"
#include "ChargeSegmentAlgos/RadialSegmentSearch.h"

namespace larlitecv {

  bool RadialEndpointFilter::isWithinStraightSegment( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
						      const float segment_radius, const float segment_min_width, const float segment_frac_w_charge,
						      const std::vector<float>& pixel_thresholds, const float acceptance_angle ) {

    larlitecv::RadialSegmentSearch segalgo;
    std::vector< Segment3D_t > seg3d_v = segalgo.find3Dsegments( img_v, badch_v, pos3d, segment_radius, pixel_thresholds, segment_min_width, segment_frac_w_charge );

    
    return true;
  }

  bool RadialEndpointFilter::isWithinStraightSegment( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
						      const RadialEndpointFilterConfig& config ) {
    return isWithinStraightSegment( pos3d, img_v, badch_v, config.segment_radius, config.segment_min_width, config.segment_frac_w_charge, config.pixel_thresholds, config.acceptance_angle );
  }
  

  
}
