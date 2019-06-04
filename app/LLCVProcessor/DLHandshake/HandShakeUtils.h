#ifndef __HANDSHAKEUTILS_H__
#define __HANDSHAKEUTILS_H__

#include "LArOpenCV/ImageCluster/Base/ImageClusterTypes.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventPixel2D.h"

namespace llcv {

  larocv::GEO2D_Contour_t as_contour(const larcv::Pixel2DCluster& ctor);
  
  larocv::GEO2D_ContourArray_t as_contour_array(const std::vector<larcv::Pixel2DCluster>& ctor_v);

  larcv::Image2D contour_to_image(const larcv::Pixel2DCluster& contour, 
				  const larcv::ImageMeta& meta);

  std::vector<larcv::Image2D> as_image_array(const std::vector<larcv::Pixel2DCluster>& contour_v,
					     const std::vector<larcv::ImageMeta>& meta_v);

}
#endif
