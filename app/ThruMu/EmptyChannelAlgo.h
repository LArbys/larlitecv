#ifndef __EMPTY_CHANNEL_ALGO__
#define __EMPTY_CHANNEL_ALGO__

#include <vector>

// larlite
#include "Base/DataFormatConstants.h"
#include "DataFormat/chstatus.h"

// larcv
#include "DataFormat/Image2D.h"


namespace larlitecv {

  class EmptyChannelAlgo {
  public:
    
    EmptyChannelAlgo() {};
    virtual ~EmptyChannelAlgo() {};

    larcv::Image2D labelEmptyChannels( float threshold, const larcv::Image2D& tpcimg );
    std::vector<larcv::Image2D> makeBadChImage( int minstatus, int nplanes, int start_tick, int nticks, int nchannels, 
						int time_downsample_factor, int wire_downsample_factor,
						const larlite::event_chstatus& ev_status );
    
  };
  
}

#endif

