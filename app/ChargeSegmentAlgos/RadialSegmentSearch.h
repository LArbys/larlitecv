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

    // Primary Method
    // --------------
    std::vector< Segment3D_t > find3Dsegments( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					       const std::vector<float>& pos3d, const float search_radius, const std::vector<float>& pixel_thresholds,
					       const int min_hit_width, const int hit_neighborhood, const float segment_frac_w_charge, int verbosity=0 );


    // Available as public in case other algorithms find useful
    // ---------------------------------------------------------
    std::vector< std::vector<larlitecv::RadialHit_t > > findIntersectingChargeClusters( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
											const std::vector<float>& pos3d, const float radius,
											const std::vector<float>& thresholds );

    std::vector< larlitecv::RadialHit_t > findIntersectingChargeClusters( const larcv::Image2D& img, const larcv::Image2D& badch,
									  const float threshold, const std::vector<larcv::Pixel2D>& pixelring );

    std::vector< Segment2D_t > make2Dsegments( const larcv::Image2D& img, const larcv::Image2D& badch, const std::vector<RadialHit_t>& hitlist,
					       const std::vector<float>& pos3d, const float threshold, const int min_hit_width, const float frac_w_charges, const int verbosity=0 );

    std::vector< std::vector< larcv::Pixel2D > > defineProjectedDiamondBoundary( const std::vector<float>& center3d,
										 const std::vector<larcv::Image2D>& img_v, const float radius );

    std::vector< std::vector< larcv::Pixel2D > > defineProjectedBoxBoundary( const std::vector<float>& center3d, const std::vector<larcv::Image2D>& img_v, const float radius );

    // intermediate algo data, provided for debugging: filled when find3Dsegments is called
    // -------------------------------------------------------------------------------------
    std::vector< std::vector<RadialHit_t> > m_planehits; //< hits found along the box/diamond boundary for each plane


  };

}

#endif
