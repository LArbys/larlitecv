#ifndef __STOPMU_SKELETON__
#define __STOPMU_SKELETON__

#include "DataFormat/Image2D.h"

namespace larlitecv {

  class StopMuSkeleton {
  public:
    StopMuSkeleton() {};
    virtual ~StopMuSkeleton() {};
    
    larcv::Image2D skeletonize( const larcv::Image2D& img, const float thresh, const int kernel_size );
    
  };
  
}

#endif
