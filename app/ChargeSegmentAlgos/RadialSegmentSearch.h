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

#include "RadialSegmentSearchTypes.h"

namespace larlitecv {

  class RadialSegmentSearch {
  public:
    RadialSegmentSearch() {};
    virtual ~RadialSegmentSearch() {};

    std::vector< larlitecv::RadialHit_t > findIntersectingChargeClusters( const larcv::Image2D& img, const larcv::Image2D& badch,
									  int row, int col, int radius, float threshold, int arc_divs );
    std::vector< larlitecv::RadialHit_t > findIntersectingChargeClustersX( const larcv::Image2D& img, const larcv::Image2D& badch,
									   const std::vector<float>& pos3d, float radius, float threshold );


  };

}

#endif
