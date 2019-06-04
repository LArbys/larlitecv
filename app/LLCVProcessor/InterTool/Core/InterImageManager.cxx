#ifndef __INTERIMAGEMANAGER_CXX__
#define __INTERIMAGEMANAGER_CXX__

#include "InterImageManager.h"
#include "InterImageManager.imp.h"

#include "InterTool_Util/InterImageUtils.h"

#include <cassert>

namespace llcv {

  void InterImageManager::SetIIT(InterImageType iitype,const std::pair<int,int>& cpair) {
    switch(iitype) {
    case kImageADC:
      _iimg_v = &(_inter_adc_m[cpair]);
      break;
    case kImageTrack:
      _iimg_v = &(_inter_trk_m[cpair]);
      break;
    case kImageShower:
      _iimg_v = &(_inter_shr_m[cpair]);
      break;
    case kImageDead:
      _iimg_v = &(_inter_dead_m[cpair]);
      break;
    default: 
      throw llcv_err("invalid InterImageType requested");
    }
  }
  
  void InterImageManager::SetIIT(InterImageType iitype,int cropx, int cropy) {
    SetIIT(iitype,std::make_pair(cropx,cropy));
  }

  void InterImageManager::InitializeOCVImage(InterImageType iitype) {
    LLCV_DEBUG() << "start" << std::endl;     

    SetIIT(iitype,-1,-1);
    
    for(size_t plane=0; plane<3; ++plane) {
      auto& iimg  = (*_iimg_v).at(plane);

      auto mat_meta = _larmkr.ExtractImage(iimg.img2d,plane);
      iimg.mat  = std::move(std::get<0>(mat_meta));
      iimg.meta = std::get<1>(mat_meta);
    }

    LLCV_DEBUG() << "end" << std::endl;     
    return;
  }
  
  void InterImageManager::CropImage(const std::pair<int,int>& cpair, InterImageType iitype) {
    CropImage(cpair.first,cpair.second,iitype);
  }

  void InterImageManager::CropImage(int cropx, int cropy, InterImageType iitype) {
    LLCV_DEBUG() << "start" << std::endl;         
    
    SetIIT(iitype,-1,-1);
    const auto& orig_iimg_v = (*_iimg_v);
    if (orig_iimg_v.size()!=3) {
      LLCV_CRITICAL() << "@iitype=" << (int)iitype 
		      << " orig_iimg_v sz=" << orig_iimg_v.size() << std::endl;
      throw llcv_err("die");
    }

    SetIIT(iitype,cropx,cropy);
    _iimg_v->resize(3);

    // crop the image by defined number of pixels left and right in user config
    for(size_t plane=0; plane<3; ++plane) {

      auto row = _vtx_pixel_v[plane].first;
      auto col = _vtx_pixel_v[plane].second;

      const auto& img  = orig_iimg_v.at(plane).img2d;

      const auto& meta = img.meta();

      LLCV_DEBUG() << "@plane=" << plane << std::endl;
      LLCV_DEBUG() << "(row,col)="<< row << "," << col << ")" << std::endl;
      LLCV_DEBUG() << "pixel (width,height)=(" << meta.pixel_width() << "," << meta.pixel_height() << ")" << std::endl;

      double width  = cropx * meta.pixel_width();
      double height = cropy * meta.pixel_height();

      LLCV_DEBUG() << "(width,height)=("<< width << "," << height << ")" << std::endl;

      size_t row_count = (size_t) cropx;
      size_t col_count = (size_t) cropy;

      LLCV_DEBUG() << "(rows,cols)=("<< row_count << "," << col_count << ")" << std::endl;
      LLCV_DEBUG() << "origin (x,y)=( " << meta.tl().x << "," << meta.tl().y << ")" << std::endl;

      double origin_x = meta.tl().x + meta.pixel_width()  * ((double)col);
      double origin_y = meta.tl().y - meta.pixel_height() * ((double)row);

      LLCV_DEBUG() << "0 (" << origin_x << "," << origin_y << ")" << std::endl;

      origin_x -= width/2.0;
      origin_y += height/2.0;

      LLCV_DEBUG() << "tl: " << meta.tl().x << "," << meta.tl().y << std::endl;
      LLCV_DEBUG() << "tr: " << meta.tr().x << "," << meta.tr().y << std::endl;
      LLCV_DEBUG() << "bl: " << meta.bl().x << "," << meta.bl().y << std::endl;
      LLCV_DEBUG() << "br: " << meta.br().x << "," << meta.br().y << std::endl;
	
      // check if origin is on the edge, move it to the edge if needed
      if (origin_x < meta.tl().x) origin_x = meta.tl().x;
      if (origin_x > meta.tr().x) origin_x = meta.tr().x;

      if (origin_y > meta.tl().y) origin_y = meta.tl().y;
      if (origin_y < meta.br().y) origin_y = meta.br().y;

      // check if origin on the other edge, move if on the edge
      LLCV_DEBUG() << "2 (" << origin_x << "," << origin_y << ")" << std::endl;

      auto max_x = origin_x + width;
      auto min_y = origin_y - height;

      if (max_x > meta.max_x()) {
	auto dist = meta.max_x() - max_x;
	origin_x += dist;
      }

      if (min_y < meta.min_y()) {
	auto dist = meta.min_y() - min_y;
	origin_y += dist;
      }

      LLCV_DEBUG() << "3 (" << origin_x << "," << origin_y << ")" << std::endl;

      larcv::ImageMeta crop_meta(width,height,
				 row_count,col_count,
				 origin_x,origin_y,
				 plane);
	
      LLCV_DEBUG() << meta.dump();
      LLCV_DEBUG() << crop_meta.dump();

      auto& crop_iimg = (*_iimg_v).at(plane);
      
      crop_iimg.img2d = img.crop(crop_meta);

      auto mat_meta = _larmkr.ExtractImage(crop_iimg.img2d,plane);

      crop_iimg.mat  = std::move(std::get<0>(mat_meta));
      crop_iimg.meta = std::get<1>(mat_meta);
    }

    LLCV_DEBUG() << "end" << std::endl;     
    return;
  }

  void InterImageManager::SetImage(const std::vector<larcv::Image2D>& img_v, llcv::InterImageType iitype) {
    
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "InterImageType=" << (int)iitype << std::endl;

    std::pair<int,int> cpair = std::make_pair(-1,-1);

    SetIIT(iitype,cpair);

    if ((*_iimg_v).size()<3) 
      (*_iimg_v).resize(3);

    for(size_t plane=0; plane<3; ++plane) 
      (*_iimg_v)[plane].img2d = img_v[plane];
    
    InitializeOCVImage(iitype);

    LLCV_DEBUG() << "end" << std::endl; 
    return;
  }

  void InterImageManager::Reset() {
    LLCV_DEBUG() << "start" << std::endl; 
    _inter_adc_m.clear();
    _inter_shr_m.clear();
    _inter_trk_m.clear();
    _inter_dead_m.clear();

    _vtx_pixel_v.clear();
    _vtx_pixel_v.resize(3);
    LLCV_DEBUG() << "end" << std::endl; 
  }

  
  void InterImageManager::SetVertex(float x, float y, float z) {
    
    SetIIT(kImageADC,-1,-1);
    
    Erase();

    static int xx;
    static int yy;

    for(size_t plane=0; plane<3; ++plane) {
      xx = yy = kINVALID_INT;
      const auto& meta = (*_iimg_v)[plane].img2d.meta();
      ProjectImage2D(meta,x,y,z,xx,yy);
      SetPixel(yy,xx,plane);
      LLCV_DEBUG() << "@plane=" << plane << " (" << yy << "," << xx << ")" << std::endl;
    }
    
  }

  void InterImageManager::SetPixel(int row, int col, size_t plane) {
    if (plane>=_vtx_pixel_v.size()) throw llcv_err("Invalid plane");
    _vtx_pixel_v[plane] = std::make_pair(row,col);
  }

  void InterImageManager::Erase() {
    if (_inter_adc_m.size()>1) {
      auto adc_iter = std::begin(_inter_adc_m);
      std::advance(adc_iter, 1);
      _inter_adc_m.erase(adc_iter,std::end(_inter_adc_m));
    }

    if (_inter_shr_m.size()>1) {
      auto shr_iter = std::begin(_inter_shr_m);
      std::advance(shr_iter, 1);
      _inter_shr_m.erase(shr_iter,std::end(_inter_shr_m));
    }

    if (_inter_trk_m.size()>1) {
      auto trk_iter = std::begin(_inter_trk_m);
      std::advance(trk_iter, 1);
      _inter_trk_m.erase(trk_iter,std::end(_inter_trk_m));
    }

    if (_inter_dead_m.size()>1) {
      auto dead_iter = std::begin(_inter_dead_m);
      std::advance(dead_iter,1);
      _inter_dead_m.erase(dead_iter,std::end(_inter_dead_m));
    }
      
  }

  template std::vector<cv::Mat*> InterImageManager::Image<cv::Mat>(llcv::InterImageType iitype, int cropx, int cropy);
  template std::vector<larocv::ImageMeta*> InterImageManager::Image<larocv::ImageMeta>(llcv::InterImageType iitype, int cropx, int cropy);
  template std::vector<larcv::Image2D*> InterImageManager::Image<larcv::Image2D>(llcv::InterImageType iitype, int cropx, int cropy);
  
}

#endif
