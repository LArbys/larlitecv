#ifndef __INTERIMAGEUTILS_CXX__
#define __INTERIMAGEUTILS_CXX__

#include "InterImageUtils.h"

//ll
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"

#include <stdexcept>

namespace llcv {

  cv::Mat As8UC3(const cv::Mat& img,int channel) {
    cv::Mat matc(img.rows,img.cols,CV_8UC3);
    uchar val0,val1,val2;
    val0 = val1 = val2 = 0;
    for(int row=0;row<img.rows;++row) {
      for(int col=0;col<img.cols;++col) { 
	int val = (int)((uchar)img.at<uchar>(row,col));
	val0 = val1 = val2 = 0;
	switch (channel) {
	case -1: val0=val1=val2=(uchar)val; break;
	case  0: val0=(uchar)val; break;
	case  1: val1=(uchar)val; break;
	case  2: val2=(uchar)val; break;
	default: throw std::runtime_error("case must be -1,0,1,2");
	}
	matc.at<cv::Vec3b>(row,col) = {val0,val1,val2};
      }
    }
    return matc;
  }
  
  
  void ProjectImage2D(const larcv::ImageMeta& meta,
		      double x,double y,double z,
		      int& xx, int& yy) {

    auto geohelp = larutil::GeometryHelper::GetME();
    auto larpro  = larutil::LArProperties::GetME();

    auto plane = meta.plane();

    auto pt2D = geohelp->Point_3Dto2D(x, y, z, plane );

    double x_compression  = (double)meta.width()  / (double)meta.cols();
    double y_compression  = (double)meta.height() / (double)meta.rows();
    double xpixel = (pt2D.w/geohelp->WireToCm() - meta.tl().x) / x_compression;
    double ypixel = (((x/larpro->DriftVelocity())*2+3200)-meta.br().y)/y_compression;

    xx = (int)(xpixel+0.5);
    yy = (int)(ypixel+0.5);

    yy = meta.rows() - yy - 1;

    return ;
  }

  void ProjectMat(const larcv::ImageMeta& meta,
		  double x,double y,double z,
		  int& xx, int& yy) {

    auto geohelp = larutil::GeometryHelper::GetME();
    auto larpro  = larutil::LArProperties::GetME();

    auto plane = meta.plane();

    auto pt2D = geohelp->Point_3Dto2D(x, y, z, plane );

    double x_compression  = (double)meta.width()  / (double)meta.cols();
    double y_compression  = (double)meta.height() / (double)meta.rows();
    double xpixel = (pt2D.w/geohelp->WireToCm() - meta.tl().x) / x_compression;
    double ypixel = (((x/larpro->DriftVelocity())*2+3200)-meta.br().y)/y_compression;

    xx = (int)(xpixel+0.5);
    yy = (int)(ypixel+0.5);

    return;
  }

  larocv::GEO2D_Contour_t AsContour(const larcv::Pixel2DCluster& pcluster,
				    const larocv::ImageMeta& meta) {
    larocv::GEO2D_Contour_t ctor;
    ctor.resize(pcluster.size());
    for(size_t pt_idx=0; pt_idx < pcluster.size(); ++pt_idx) {
      auto const& pt = pcluster[pt_idx];
      
      float x = pt.X();
      float y = pt.Y();

      x -= meta.origin().y;
      y -= meta.origin().x;
	
      //x = x / meta.pixel_height();
      y = y / meta.pixel_height();

      ctor[pt_idx].x = y;
      ctor[pt_idx].y = x;
    }

    return ctor;
  }
  

  geo2d::Vector<int> MatToImage2D(const geo2d::Vector<int>& pt, const cv::Mat& img) {
    geo2d::Vector<int> ret(-1,-1);

    int x = img.cols - pt.x - 1;
    int y = pt.y;

    ret.x = x;
    ret.y = y;

    return ret;
  }

  float MatToImage2DPixel(const geo2d::Vector<int>& pt, const cv::Mat& img, const larcv::Image2D& img2d) {
    float ret = -1;
    
    auto pt2 = MatToImage2D(pt,img);
    
    auto row = pt2.x;
    auto col = pt2.y;    
    
    ret = img2d.pixel(row,col);

    return ret;
  }

}

#endif
