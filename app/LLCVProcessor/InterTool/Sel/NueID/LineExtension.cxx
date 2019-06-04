#ifndef __LINE_EXTENSION_CXX__
#define __LINE_EXTENSION_CXX__

#include "LineExtension.h"

#include "LLCVBase/llcv_err.h"

namespace llcv {

  void LineExtension::SetImageDimension(const cv::Mat& img, const cv::Mat& dead) {

    if (_img.rows != img.rows) {
      _img       = cv::Mat::zeros(img.rows,img.cols, CV_8UC1);
      _white_img = cv::Mat::ones(img.rows,img.cols, CV_8UC1);
    }
    else {
      _img.setTo(cv::Scalar(0));
      _white_img.setTo(cv::Scalar(1));
    }
    
    for(int row=0; row<img.rows; row++) {
      for(int col=0; col<img.cols; col++) {

	// set the adc value as 3
	if ( ((int)img.at<uchar>(row,col)) > 0 )
	  _img.at<uchar>(row,col) = (uchar)3;
	
	// set the dead value as 2
	if ( ((int)dead.at<uchar>(row,col)) == 0)
	  _img.at<uchar>(row,col) = (uchar)2;
      }
    }
    
    return;
  }

  void LineExtension::SetCosmicPixels(const std::vector<larocv::GEO2D_Contour_t>& cosmic_ctor_v) {
    for(const auto& cosmic_ctor : cosmic_ctor_v) {
      auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(_white_img,cosmic_ctor,-1,false));
      for(const auto& nz_pt : nz_pt_v) 
	_img.at<uchar>(nz_pt.y, nz_pt.x) = (uchar)1;
    }
    return;
  }

  bool LineExtension::AtBoundary(const geo2d::Vector<float>& pt) const {

    static std::vector<int> row_v(3);
    static std::vector<int> col_v(3);

    int maxrow = _img.rows;
    int maxcol = _img.cols;

    int row = (int)(pt.y + 0.5);
    int col = (int)(pt.x + 0.5);
    
    row_v[0] = row;
    row_v[1] = row+1;
    row_v[2] = row-1;

    col_v[0] = col;
    col_v[1] = col+1;
    col_v[2] = col-1;

    for(auto& v : row_v) { 
      if(v<0)       v=0;
      if(v>=maxrow) v=maxrow-1;
    }

    for(auto& v : col_v) { 
      if(v<0)       v=0;
      if(v>=maxcol) v=maxcol-1;
    }
    
    bool ret = false;
    
    for(size_t rid=0; rid<3; ++rid) {
      for(size_t cid=0; cid<3; ++cid) {

	int val = (int)_img.at<uchar>(row_v[rid],col_v[cid]);
	
	if (val == 1 or val == 2)
	  ret = true;
      }
    }

    return ret;
  }


  bool LineExtension::AtParticle(const geo2d::Vector<float>& pt) const {
    
    static std::vector<int> row_v(3);
    static std::vector<int> col_v(3);

    int maxrow = _img.rows;
    int maxcol = _img.cols;

    int row = (int)(pt.y + 0.5);
    int col = (int)(pt.x + 0.5);
    
    row_v[0] = row;
    row_v[1] = row+1;
    row_v[2] = row-1;

    col_v[0] = col;
    col_v[1] = col+1;
    col_v[2] = col-1;

    for(auto& v : row_v) { 
      if(v<0)       v=0;
      if(v>=maxrow) v=maxrow-1;
    }

    for(auto& v : col_v) { 
      if(v<0)       v=0;
      if(v>=maxcol) v=maxcol-1;
    }
    
    bool ret = false;
    
    for(size_t rid=0; rid<3; ++rid) {
      for(size_t cid=0; cid<3; ++cid) {

	int val = (int)_img.at<uchar>(row_v[rid],col_v[cid]);
	
	if (val == 3)
	  ret = true;
      }
    }

    return ret;
  }

  bool LineExtension::InBoundary(const geo2d::Vector<float>& pt) const {
    bool ret = false;

    geo2d::Vector<float> pt_i((int)(pt.x+0.5),(int)(pt.y+0.5));

    if (!larocv::Contained(_img,pt_i)) return ret;
    
    int val = (int)_img.at<uchar>(pt_i.y,pt_i.x);

    if (val == 0) {
      return ret;
    }
    else if (val == 3) {
      return ret;
    }
    else if (val == 1 or val == 2) {
      ret = true;
    }
    else {
      throw llcv_err("invalid pixel value observed");
    }
    
    return ret;
  }


  bool LineExtension::InParticle(const geo2d::Vector<float>& pt) const {
    bool ret = false;

    geo2d::Vector<float> pt_i((int)(pt.x+0.5),(int)(pt.y+0.5));

    if (!larocv::Contained(_img,pt_i)) return ret;
    
    int val = (int)_img.at<uchar>(pt_i.y,pt_i.x);

    if (val == 0) {
      return ret;
    }
    else if (val == 3) {
      ret = true;
    }
    else if (val == 1 or val == 2) {
      return ret;
    }
    else {
      throw llcv_err("invalid pixel value observed");
    }
    
    return ret;
  }
  
  
  geo2d::Vector<float> LineExtension::ExtendAcross(const geo2d::Vector<float> start, const geo2d::Vector<float> edge) const {

    // no need to extend if not at edge
    if (!AtBoundary(edge))
      return edge;

    geo2d::Line<float> line(start,edge - start);

    auto ret = edge;

    bool right = false;
    bool isinf = false;
    bool up    = false;

    if (edge.x > start.x) {
      right = true;
    }
    else if (edge.x == start.x) {
      isinf = true;
      if (edge.y > start.y)
	up = true;
    }
    
    // step inside boundary
    if (!isinf) {
      if (right) {
	ret.x += 1;
	ret.y  = line.y(ret.x);
      }
      else {
	ret.x -= 1;
	ret.y = line.y(ret.x);
      }
    } else {
      if (up) {
	ret.y += 1;
      }
      else {
	ret.y -= 1;
      }
    }
    
    bool cross_particle = false;
    
    geo2d::Vector<float> pre_pt(-1,-1);

    // while inside boundary or at particle, keep stepping
    while(InBoundary(ret) or AtParticle(ret)) {

      if (!isinf) {
	if (right) {
	  ret.x += 1;
	  ret.y  = line.y(ret.x);
	}
	else {
	  ret.x -= 1;
	  ret.y = line.y(ret.x);
	}
      } else {
	if (up) {
	  ret.y += 1;
	}
	else {
	  ret.y -= 1;
	}
      }

      if(!larocv::Contained(_img,ret)) break;

      cross_particle |= AtParticle(ret);
      
      if (AtParticle(ret))
	pre_pt = ret;
      
    }

    // cross a particle? no - return same edge
    if(!cross_particle)
      return edge;

    if (pre_pt.x > 0)
      ret = pre_pt;

    return ret;
  }
  
  
}
#endif
