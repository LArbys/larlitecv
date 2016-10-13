#ifndef __EMPTY_CHANNEL_ALGO__
#define __EMPTY_CHANNEL_ALGO__

// larcv
#include "DataFormat/Image2D.h"

namespace larlitecv {

  class EmptyChannelAlgo {
  public:
    
    EmptyChannelAlgo() {};
    virtual ~EmptyChannelAlgo() {};

    larcv::Image2D labelEmptyChannels( float threshold, const larcv::Image2D& tpcimg );


  };

}

#endif

