#ifndef __LINE_REGION_TEST_H__
#define __LINE_REGION_TEST_H__

#include "TrackTests.h"

#include <vector>

#include "DataFormat/Image2D.h"
#include "BoundaryEndPt.h"
#include "BoundarySpacePoint.h"
#include "BMTrackCluster2D.h"

namespace larlitecv {

  class LineRegionTest : public TrackTestBase {
    
  public:

    LineRegionTest(int regionwidth=10, float frac_threshold=0.9, float pix_threshold=10.0 ) : TrackTestBase() 
      { 
	fRegionWidth = regionwidth; 
	fFractionThreshold = frac_threshold;  
	fPixelThreshold = pix_threshold;
	verbose_debug = false;
	for (int i=0; i<3; i++)
	  last_fractions[i] = 0.;
      };
    virtual  ~LineRegionTest() {};

    virtual bool test( const BoundarySpacePoint& start_v, const BoundarySpacePoint& end_v,
                       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                       std::vector< BMTrackCluster2D >* opt_track2d=NULL );
    int fRegionWidth;
    float fFractionThreshold;
    float fPixelThreshold;
    float last_fractions[3];
    bool verbose_debug;
    
  };

}

#endif
