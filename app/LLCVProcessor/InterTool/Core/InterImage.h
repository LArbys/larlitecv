#ifndef __INTERIMAGE_H__
#define __INTERIMAGE_H__

#include "LArOpenCV/Core/ImageMeta.h"
#include "DataFormat/Image2D.h"

namespace llcv {

  class InterImage {
  public:
    
    cv::Mat mat;
    larocv::ImageMeta meta;
    larcv::Image2D img2d;
    
    template<class T> T* get() { return nullptr; }
    
  };

  template<> inline cv::Mat* InterImage::get<cv::Mat> ()  { return &mat; }
  template<> inline larocv::ImageMeta* InterImage::get<larocv::ImageMeta> () { return &meta; }
  template<> inline larcv::Image2D* InterImage::get<larcv::Image2D> () { return &img2d; }

}


#endif
