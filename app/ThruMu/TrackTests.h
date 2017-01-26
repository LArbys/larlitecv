#ifndef __TRACK_TESTS_BASE__H__
#define __TRACK_TESTS_BASE__H__

#include <vector>
#include "DataFormat/Image2D.h"
#include "BoundaryEndPt.h"
#include "BoundarySpacePoint.h"
#include "BMTrackCluster2D.h"

namespace larlitecv {
  
  class TrackTestBase {

  public:
    TrackTestBase() {};
    virtual ~TrackTestBase() {};

    virtual bool test( const BoundarySpacePoint& start_v, const BoundarySpacePoint& end_v,
		       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
		       std::vector< BMTrackCluster2D >* opt_track2d=NULL ) = 0;
    
  };


}

#endif
