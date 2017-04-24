#ifndef __RADIALSEGMENTSEARCH__
#define __RADIALSEGMENTSEARCH__

/*
 * basic idea here is given a 3D point, we project back into the planes, 
 * draw a 2D circle around each projection point, then find common intersections
 * in those circles that are common in time.
 */

#include <vector>

#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2D.h"

#include "Segment3DAlgo.h"
#include "RadialSegmentSearchTypes.h"

namespace larlitecv {

  class RadialSegmentSearch {
  public:
    RadialSegmentSearch() {};
    virtual ~RadialSegmentSearch() {};

    std::vector< larlitecv::RadialHit_t > findIntersectingChargeClusters( const larcv::Image2D& img, const larcv::Image2D& badch, const std::vector<float>& pos3d, float radius, float threshold );

    std::vector< Segment2D_t > make2Dsegments( const larcv::Image2D& img, const larcv::Image2D& badch, const std::vector<RadialHit_t>& hitlist,
					       const std::vector<float>& pos3d, const float threshold, const int min_hit_width, const float frac_w_charges );
    

    std::vector< Segment3D_t > find3Dsegments( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					       const std::vector<float>& pos3d, const float search_radius, const std::vector<float>& pixel_thresholds,
					       const int min_hit_width, const float segment_frac_w_charge );
    

  };

}

#endif
