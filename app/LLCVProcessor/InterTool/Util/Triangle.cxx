#ifndef __TRIANGLE_CXX__
#define __TRIANGLE_CXX__

#include "Triangle.h"

#include <cassert>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include "Geo2D/Core/Geo2D.h"


namespace llcv { 
  
  Triangle::Triangle(const larocv::GEO2D_Contour_t& ctor,
		     const geo2d::Vector<float>& start) {
    _ctor = ctor;
    _apex = start;

    Construct();
  }

  void Triangle::Construct() {

    auto min_rect  = larocv::MinAreaRect(_ctor);
    geo2d::Vector<float> pt_v[4];
    min_rect.points(pt_v);

    // get the two closest points to the vertex
    float dist = larocv::kINVALID_FLOAT;
    
    size_t id1 = larocv::kINVALID_SIZE;
    size_t id2 = larocv::kINVALID_SIZE;

    for(size_t pid=0; pid<4; ++pid) {
      auto dist_tmp = geo2d::dist(_apex,pt_v[pid]);
      if (dist_tmp < dist) {
	dist = dist_tmp;
	id1 = pid;
      }
    }

    dist = larocv::kINVALID_FLOAT;

    for(size_t pid=0; pid<4; ++pid) {
      if (pid == id1) continue;
      auto dist_tmp = geo2d::dist(_apex,pt_v[pid]);
      if (dist_tmp < dist) {
	dist = dist_tmp;
	id2 = pid;
      }
    }
    
    geo2d::Vector<float> far_pt1(larocv::kINVALID_FLOAT,
				 larocv::kINVALID_FLOAT);
    geo2d::Vector<float> far_pt2(larocv::kINVALID_FLOAT,
				 larocv::kINVALID_FLOAT);

    for(size_t pid=0; pid<4; ++pid) {
      if (pid == id1) continue;
      if (pid == id2) continue;
      if (far_pt1.x == larocv::kINVALID_FLOAT) {
	far_pt1 = pt_v[pid];
	continue;
      }
      if (far_pt2.x == larocv::kINVALID_FLOAT) {
	far_pt2 = pt_v[pid];
	continue;
      }
    }
    
    _base_pt1 = far_pt1;
    _base_pt2 = far_pt2;
  }

  void Triangle::Expand(const cv::Mat& img, const float fraction) { 

    assert(fraction>0);
    assert(fraction<1);

    //
    // get the mid point between base1 and base2
    //
    auto mid_pt = MidPoint(_base_pt1,_base_pt2);

    // 
    // find the point fraction * length
    //
    auto apex_to_mid1 = mid_pt - _apex;
    apex_to_mid1 *= fraction;
    apex_to_mid1 += _apex;

    //
    // make the ctor_img
    //
    auto ctor_img = larocv::BlankImage(img,255);
    ctor_img = larocv::MaskImage(ctor_img,_ctor,-1,false);

    auto ctor_img_noline = ctor_img.clone();

    //
    // Chop off the pixels along the apex line
    //
    auto dir = _base_pt1 - _base_pt2;
    geo2d::Vector<float> apex_to_mid2;
      
    if (dir.x != 0) {
      geo2d::Line<float> line(apex_to_mid1,dir);
      
      apex_to_mid1.x = img.rows;
      apex_to_mid1.y = line.y(apex_to_mid1.x);

      apex_to_mid2.x = 0;
      apex_to_mid2.y = line.y(apex_to_mid2.x);
    }
    // it's a vertical line
    else {
      apex_to_mid1.x = apex_to_mid1.x;
      apex_to_mid1.y = img.rows;

      apex_to_mid2.x = apex_to_mid1.x;
      apex_to_mid2.y = 0;
    }
    
    cv::line(ctor_img,apex_to_mid1,apex_to_mid2,cv::Scalar(0),2);
    
    //
    // Find contours, choose one closest to mid_pt
    //
    auto ctor_v = larocv::FindContours(ctor_img);
    size_t close_id = larocv::kINVALID_SIZE;
    float dist = larocv::kINVALID_FLOAT;
    for(size_t cid=0; cid < ctor_v.size(); ++cid) {
      auto dist_tmp = larocv::Pt2PtDistance(mid_pt,ctor_v[cid]);
      if (dist_tmp < dist) {
	dist = dist_tmp;
	close_id = cid;
      }
    }
    
    larocv::GEO2D_Contour_t far_ctor;
    if (ctor_v.empty() or close_id == larocv::kINVALID_SIZE) {
      ctor_img = ctor_img_noline.clone();
      far_ctor = _ctor;
    } else {
      far_ctor = ctor_v[close_id];
    }

    ctor_img = larocv::MaskImage(ctor_img,far_ctor,-1,false);

    MovePt(ctor_img,_base_pt1);
    MovePt(ctor_img,_base_pt2);

    return;
  }
  
  
  void Triangle::MovePt(const cv::Mat& ctor_img, geo2d::Vector<float>& base_pt) {

    geo2d::Vector<float> dir(0,0);
    bool right = false;
    bool up = false;

    // input is base 1
    if (base_pt == _base_pt1) {
      dir = base_pt - _base_pt2;

      if (base_pt.x > _base_pt2.x)
	right = true;
      if (base_pt.y > _base_pt2.y)
	up = true;
    }
    
    // input is base 2
    else { 
      dir = base_pt - _base_pt1;    

      if (base_pt.x > _base_pt1.x)
	right = true;
      if (base_pt.y > _base_pt1.y)
	up = true;
    }

    bool isinf = false;
    geo2d::Line<float> line(base_pt,dir);

    if (dir.x == 0)
      isinf = true;

    while(larocv::Broken(ctor_img,_apex,base_pt,1)) {

      if (isinf) {
	if (up)
	  base_pt.y += 1;
	else
	  base_pt.y -= 1;
      }

      else {
	if (right)
	  base_pt.x += 1;
	else
	  base_pt.x -= 1;

	base_pt.y  = line.y(base_pt.x);
      }

      if (base_pt.x >= ctor_img.rows) break;
      if (base_pt.y >= ctor_img.cols) break;

      if(base_pt.x < 0) break;
      if(base_pt.y < 0) break;
    }
    
    
    return;
  }

  void Triangle::Tighten(const cv::Mat& img, const float radius, const float fraction) { 

    //
    // Get nearby nonzero pixels around the apex
    //
    geo2d::Circle<float> circle(_apex,radius);
    auto near_v = larocv::FindNonZero(larocv::MaskImage(img,circle,-1,false));

    std::vector<float> area_v(near_v.size(),-1);

    std::vector<geo2d::Vector<float> > base1_v(near_v.size());
    std::vector<geo2d::Vector<float> > base2_v(near_v.size());

    for(size_t nid=0; nid<near_v.size(); ++nid) {

      const auto& near = near_v[nid];

      auto& area = area_v[nid];
      auto& base1 = base1_v[nid];
      auto& base2 = base2_v[nid];
      
      _apex.x = near.x;
      _apex.y = near.y;

      Construct();
      Expand(img,fraction);
      
      base1 = _base_pt1;
      base2 = _base_pt2;

      area = 0.5 * geo2d::dist(base1,base2) * geo2d::dist(_apex,MidPoint(base1,base2));
    }

    if (area_v.empty()) return;
      
    auto min_iter = std::min_element(area_v.begin(),area_v.end());
    auto min_idx = std::distance(area_v.begin(),min_iter);

    _base_pt1 = base1_v[min_idx];
    _base_pt2 = base2_v[min_idx];
    _apex.x   = near_v[min_idx].x;
    _apex.y   = near_v[min_idx].y;
  
    return;
  }

  geo2d::Vector<float> Triangle::MidPoint(const geo2d::Vector<float>& pt1,const geo2d::Vector<float>& pt2) const {
    geo2d::Vector<float> res;

    auto bmax_x = std::max(pt1.x,pt2.x);
    auto bmin_x = std::min(pt1.x,pt2.x);

    auto bmax_y = std::max(pt1.y,pt2.y);
    auto bmin_y = std::min(pt1.y,pt2.y);

    auto bdx = bmax_x - bmin_x;
    auto bdy = bmax_y - bmin_y;
    
    bdx /= 2.0;
    bdy /= 2.0;
    
    res.x = bmin_x + bdx;
    res.y = bmin_y + bdy;

    return res;
  }

  float Triangle::StraightLineTest(const cv::Mat& img) const { 

    float res = -1;

    auto mid_pt = MidPoint(_base_pt1,_base_pt2);
    auto white_img = larocv::MaskImage(img,_ctor,-1,false);
    white_img = larocv::Threshold(white_img,10,255);
    
    float nzero_before = (float) larocv::CountNonZero(white_img);
    
    cv::line(white_img,_apex,mid_pt,cv::Scalar(0),3);

    float nzero_after = (float) larocv::CountNonZero(white_img);
    
    if (nzero_before != 0) 
      res = 1.0 - (nzero_after / nzero_before);
    
    return res;
  }

  larocv::GEO2D_Contour_t Triangle::AsContour() const {
    larocv::GEO2D_Contour_t res(3);

    res[0].x = (int)_apex.x;
    res[0].y = (int)_apex.y;

    res[1].x = (int)_base_pt1.x;
    res[1].y = (int)_base_pt1.y;

    res[2].x = (int)_base_pt2.x;
    res[2].y = (int)_base_pt2.y;

    return res;
  }

  void Triangle::Extend(const float fraction) {

    assert (fraction>1);

    //
    // get the mid point length
    //
    auto mid_pt = MidPoint(_base_pt1,_base_pt2);
    
    auto mid_dir   = mid_pt - _apex;
    auto edge1_dir = _base_pt1 - _apex;
    auto edge2_dir = _base_pt2 - _apex;
    auto base_dir  = _base_pt1 - _base_pt2;

    auto new_mid = _apex + mid_dir*fraction;
    
    geo2d::Line<float> new_base_line(new_mid,base_dir);
    geo2d::Line<float> edge1_line(_base_pt1,edge1_dir);
    geo2d::Line<float> edge2_line(_base_pt2,edge2_dir);

    auto base1_ipoint = geo2d::IntersectionPoint(new_base_line,edge1_line);
    auto base2_ipoint = geo2d::IntersectionPoint(new_base_line,edge2_line);

    _base_pt1 = base1_ipoint;
    _base_pt2 = base2_ipoint;

  }

  float Triangle::Height() const {
    float res = -1;
    
    auto mid_pt = MidPoint(_base_pt1,_base_pt2);
    res = geo2d::dist(mid_pt,_apex);
    
    return res;
  }
  

  float Triangle::EmptyAreaRatio() const {
    float res = -1;

    float tri_area = larocv::ContourArea(this->AsContour());
    float ctor_area = larocv::ContourArea(_ctor);
    
    if (ctor_area > 0)
      res = tri_area / ctor_area;
    
    return res;
  }

  float Triangle::EmptyArea() const {
    float res = -1;
    
    float tri_area = larocv::ContourArea(this->AsContour());
    float ctor_area = larocv::ContourArea(_ctor);

    res = tri_area - ctor_area;
    
    return res;
  }

  float Triangle::Coverage(const cv::Mat& img) const {
    float res = -1;
    
    auto mask = larocv::MaskImage(img,_ctor,-1,false);
    float mask_area = larocv::CountNonZero(mask);
    float mask_area_triangle = larocv::CountNonZero(larocv::MaskImage(mask,this->AsContour(),-1,false));

    if (mask_area > 0)
      res = mask_area_triangle / mask_area;

    return res;
  }

  float Triangle::BaseLength() const {
    float res = -1;
    
    res = geo2d::dist(_base_pt1,_base_pt2);

    return res;
  }

  float Triangle::Area() const {
    float res = -1;

    res = 0.5 * BaseLength() * Height();
    
    return res;
  }

  Triangle Triangle::RotateToPoint(const geo2d::Vector<float>& pt,float scale) const {
    // make a new triangle
    
    auto new_dir   = pt - this->_apex;
    auto length = geo2d::length(new_dir);
    new_dir = new_dir / length;

    auto perp_dir  = geo2d::Vector<float>(new_dir.y,-1.0*new_dir.x);
    float half_base_length = this->BaseLength();
    half_base_length *= scale;
    half_base_length /= 2.0;
    
    geo2d::Vector<float> base1, base2;
    base1 = base2 = pt;

    float half_length = 0.0;
    float step = 0.1;

    // look to positive 
    half_length = 0.0;
    while(half_length < half_base_length) {
      base1 += step * perp_dir;
      half_length = geo2d::dist(base1,pt);
    }

    // look to negative
    half_length = 0.0;
    while(half_length < half_base_length) {
      base2 -= step * perp_dir;
      half_length = geo2d::dist(base2,pt);
    }

    return Triangle(this->_apex,base1,base2);
  }


}


#endif
