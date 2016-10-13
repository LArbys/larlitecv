#ifndef __BOUNDARY_INTERSECTION_ALGO_H_
#define __BOUNDARY_INTERSECTION_ALGO_H_

#include <vector>
#include "BoundaryEndPt.h"

namespace larlitecv {

  class BoundaryIntersectionAlgo {

  public:
    BoundaryIntersectionAlgo() {};
    virtual ~BoundaryIntersectionAlgo() {};


    void determine3Dpoint( const std::vector<int>& plane_wires, std::vector<float>& vertex, BoundaryEndPt::BoundaryEnd_t endpt_type ) {};

  };

}

#endif
