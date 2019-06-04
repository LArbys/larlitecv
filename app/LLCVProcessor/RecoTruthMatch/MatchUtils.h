#ifndef MATCHUTILS_H
#define MATCHUTILS_H

// larcv
#include "DataFormat/Image2D.h"

namespace llcv {

  void Project3D(const larcv::ImageMeta& meta,
		 double parent_x,
		 double parent_y,
		 double parent_z,
		 double parent_t,
		 uint plane,
		 double& xpixel, double& ypixel);


  void mask_image(larcv::Image2D& target, const larcv::Image2D& ref);

  float TestPixelType(int row, int col, 
		      const larcv::Image2D& adc_img, const larcv::Image2D& pgraph_img,
		      bool ignore_zero=false);
  
}
#endif


