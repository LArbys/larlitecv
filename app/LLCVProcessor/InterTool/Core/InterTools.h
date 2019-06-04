#ifndef __INTERTOOLS_H__
#define __INTERTOOLS_H__

#include "DataFormat/PGraph.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2DCluster.h"

#include <map>

namespace llcv {
  
  void MaskImage(const std::vector<larcv::PGraph>& pgraph_v,
  		 const std::map<larcv::PlaneID_t, std::vector<larcv::Pixel2DCluster> >& pix_m,
  		 const std::map<larcv::PlaneID_t, std::vector<larcv::ImageMeta> >&   pix_meta_m,
  		 std::vector<larcv::Image2D>& img_v);


}

#endif
