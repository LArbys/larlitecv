#ifndef __SHOWER_RECO_UTIL_H__
#define __SHOWER_RECO_UTIL_H__

#include <vector>

// larcv
#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/EventPGraph.h"


namespace larlitecv {
namespace ssnetshowerreco {

  class ShowerRecoUtil {

  public:

    ShowerRecoUtil() {};
    virtual ~ShowerRecoUtil() {};
    
    bool process_event( larcv::IOManager& iolcv,
                        std::vector< std::vector<larcv::Image2D> >& adccrop_vv,
                        std::vector< std::vector<larcv::Image2D> >& statuscrop_vv );

    void makeVertexImageCrop( const std::vector<larcv::Image2D>& adc_v,
                              const std::vector<larcv::Image2D>& shower_img_v,
                              const larcv::EventChStatus& ev_wirestatus,
                              const larcv::PGraph& vtxinfo,
                              const float adc_threshold,
                              const float showerscore_threshold,
                              std::vector<larcv::Image2D>& cropimg_v,
                              std::vector<larcv::Image2D>& statusimg_v );

    larcv::Image2D GetChStatusImage(int L, const larcv::ChStatus& ch, int leftCol);
    larcv::Image2D GetChStatusVtxImage(int L, const larcv::ChStatus& ch, int leftCol, int vtxRow, int vtxCol);    
    
  };

}
}

#endif
