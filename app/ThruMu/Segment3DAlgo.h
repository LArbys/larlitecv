#ifndef __SEGMENT3D_ALGO__
#define __SEGMENT3D_ALGO__

/*
  Segment3DAlgo: within a range of time, find 3d-consistent line segment across 3 planes
 */

// std lib
#include <vector>

// larcv
#include "DataFormat/Image2D.h"

// larlitecv
#include "Segment3DAlgoTypes.h"

namespace larlitecv {

  class Segment3DAlgo {

  public:
    
    Segment3DAlgo() {};
    virtual ~Segment3DAlgo() {};

    std::vector< Segment3D_t > find3DSegments( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					       const int row_a, const int row_b, const std::vector<float>& thresholds, const int min_hit_width );

    std::vector<int> findHits( const larcv::Image2D& img, const int row, const float& threshold, const int min_hit_width );
    
    std::vector< Segment2D_t > make2DSegments( const larcv::Image2D& img, const larcv::Image2D& badch, const int lowrow, const std::vector<int>& hits_low,
					       const int highrow, const std::vector<int>& hits_high, const float threshold,
					       const int hit_width, const float frac_good );
    
    void checkSegmentCharge( const larcv::Image2D& img, const larcv::Image2D& badch, const int low_row, const int low_col, const int high_row, const int high_col, const int hit_width,
			     const float threshold, int& nrows_w_charge, int& num_rows );

    void combine2Dinto3D( const std::vector< std::vector<Segment2D_t> >& plane_segments2d, const std::vector< larcv::Image2D >& img_v, const std::vector<larcv::Image2D>& badch_v,
			  const int hit_width, const std::vector<float>& threshold, float good_frac, std::vector<Segment3D_t>& segments );
    
  };

}

#endif
