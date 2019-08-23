#ifndef __POLYGON_CXX__
#define __POLYGON_CXX__

#include "Polygon.h"

#ifndef __CLING__
#ifndef __CINT__
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#endif
#endif

#include <cassert>
#include <stdexcept>

#include "InterImageUtils.h"

namespace llcv {
  
  Polygon::Polygon(const larocv::GEO2D_Contour_t& ctor,const geo2d::Vector<float>& start) {

    _ctor = ctor;
    _start = start;

    Construct();

  }      
      
  void Polygon::Construct() {

    _hull.clear();
    _hull_v.clear();
    _defect_info_v.clear();
    _defect_dist_v.clear();

    cv::convexHull(_ctor, _hull_v);
    _hull_v.push_back(_hull_v.front());

    cv::convexityDefects(_ctor,_hull_v,_defect_info_v);

    _defect_dist_v.resize(_defect_info_v.size(),-1);
    
    for(size_t did=0; did<_defect_info_v.size(); ++did)
      _defect_dist_v[did]  = ((float)_defect_info_v[did][3])/256.0;

    _hull.reserve(_hull_v.size());
    for(auto hid : _hull_v)
      _hull.emplace_back(_ctor.at(hid));
    
    VetoStartPoint();
  }

  void Polygon::VetoStartPoint() {
    _veto_ctor_id = larocv::kINVALID_INT;
    _veto_hull_id = larocv::kINVALID_INT;

    _veto_ctor_id = larocv::Pt2PtDistance(_start,_ctor);
    _veto_hull_id = larocv::Pt2PtDistance(_start,_hull);

    return;
  }

  int Polygon::NumberDefects(float dist_thresh) const {
    int res = 0;
    
    for(auto defect_dist : _defect_dist_v) {
      if (defect_dist > dist_thresh) {
	res += 1;
      }
    }

    return res;
  }


  int Polygon::NumberDefectsNoStart(float dist_thresh) const {
    int res = 0;
    
    for(size_t did=0; did<_defect_dist_v.size(); ++did) {
      
      auto start_idx = _defect_info_v[did][0];
      auto end_idx   = _defect_info_v[did][1];

      if (_veto_ctor_id >= start_idx and _veto_ctor_id <= end_idx) continue;
      
      auto defect_dist = _defect_dist_v[did];
      if (defect_dist > dist_thresh)
	res += 1;

    }

    return res;
  }


  float Polygon::LargestDefect() const {
    float res = 0;

    if (_defect_dist_v.empty()) return res;

    auto max_iter = std::max_element(_defect_dist_v.begin(),_defect_dist_v.end());
    
    res = *max_iter;
    
    return res;
  }

  float Polygon::LargestDefectNoStart() const {
    float res = 0;

    float largest = -1*larocv::kINVALID_FLOAT;
    
    if (_defect_dist_v.empty()) return res;

    for(size_t did=0; did<_defect_dist_v.size(); ++did) {
      
      auto start_idx = _defect_info_v[did][0];
      auto end_idx   = _defect_info_v[did][1];

      if (_veto_ctor_id >= start_idx and _veto_ctor_id <= end_idx) continue;

      auto defect_dist = _defect_dist_v[did];
      largest = std::max(largest, defect_dist);
    }

    res = largest;
    
    return res;
  }


  float Polygon::SmallestDefect() const {
    float res = 0;
    
    if (_defect_dist_v.empty()) return res;

    auto min_iter = std::min_element(_defect_dist_v.begin(),_defect_dist_v.end());

    res = *min_iter;

    return res;
  }

  float Polygon::SmallestDefectNoStart() const {
    float res = 0;

    float smallest = larocv::kINVALID_FLOAT;

    if (_defect_dist_v.empty()) return res;

    for(size_t did=0; did<_defect_dist_v.size(); ++did) {
      
      auto start_idx = _defect_info_v[did][0];
      auto end_idx   = _defect_info_v[did][1];

      if (_veto_ctor_id >= start_idx and _veto_ctor_id <= end_idx) continue;
      
      auto defect_dist = _defect_dist_v[did];
      smallest = std::min(smallest,defect_dist);
    }

    res = smallest;
    
    return res;
  }

  float Polygon::EmptyAreaRatio() const {
    float res = -1;
    
    float hull_area = larocv::ContourArea(_hull);
    float ctor_area = larocv::ContourArea(_ctor);
    
    if (ctor_area > 0)
      res = hull_area / ctor_area;
      
    
    return res;
  }


  float Polygon::EmptyArea() const {
    float res = -1;
    
    float hull_area = larocv::ContourArea(_hull);
    float ctor_area = larocv::ContourArea(_ctor);
    
    res = hull_area - ctor_area;
    
    return res;
  }

  float Polygon::PocketArea() const {
    float res = 0;

    for(const auto& defect_info : _defect_info_v) {
      larocv::GEO2D_Contour_t pocket;
      //pocket.reserve(defect_info[1] - defect_info[0]);

      for(int did=defect_info[0]; did<defect_info[1]; ++did)
	pocket.emplace_back(_ctor.at(did));
      
      if (!pocket.empty()) {
	pocket.emplace_back(pocket.front());
	res += larocv::ContourArea(pocket);
      }

    }

    return res;
  }

  float Polygon::PocketAreaNoStart() const {
    float res = 0;

    for(size_t dinfo=0; dinfo<_defect_info_v.size(); ++dinfo) {

      auto start_idx = _defect_info_v[dinfo][0];
      auto end_idx   = _defect_info_v[dinfo][1];

      if (_veto_ctor_id >= start_idx and _veto_ctor_id <= end_idx) continue;

      larocv::GEO2D_Contour_t pocket;
      //pocket.reserve(end_idx - start_idx);

      for(int did=start_idx; did<end_idx; ++did)
	pocket.emplace_back(_ctor.at(did));
      
      if (!pocket.empty()) {
	pocket.emplace_back(pocket.front());
	res += larocv::ContourArea(pocket);
      }
    }

    return res;
  }

  float Polygon::Area() const {
    float ret = -1.0*larocv::kINVALID_FLOAT;
    ret = larocv::ContourArea(_ctor);
    return ret;
  }

  float Polygon::Perimeter() const {
    float ret = -1.0*larocv::kINVALID_FLOAT;
    ret = larocv::ArcLength(_ctor);
    return ret;
  }
  
  float Polygon::Charge(const larcv::Image2D& img2d, const cv::Mat& img) const {
    float ret = -1.0*larocv::kINVALID_FLOAT;
    auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(img,_ctor,-1,false));

    ret = 0;
    for(const auto& nz_pt : nz_pt_v) {
      float charge = MatToImage2DPixel(nz_pt,img,img2d); 
      if (charge<=0) continue;
      ret += charge;
    }

    return ret;
  }

  void Polygon::DetectBranching(const cv::Mat& img, 
				const float rad,
				const float thickness,
				const int edge_size,
				const int branch_size) {

    _edge_pt_v.clear();
    _branch_pt_v.clear();

    if (_ctor.empty()) return;
    
    std::vector<std::vector<const geo2d::Vector<int>*> > xs_vv;
    xs_vv.clear();

    auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(img,_ctor,-1,false));

    for(const auto& nz_pt : nz_pt_v) {
      geo2d::Circle<float> circle(nz_pt.x, nz_pt.y,rad);
      auto ret_v = larocv::OnCircleGroups(img,circle);
      int ctr=0;
      for(const auto& ret : ret_v) {
	if (!larocv::Connected(img,ret,circle.center,thickness))
	  continue;
	ctr++;
      }

      if (ctr >= (int)xs_vv.size())
	xs_vv.resize(ctr+1);

      auto& xs_v = xs_vv[ctr];
      xs_v.push_back(&nz_pt);
    }
    
    if (xs_vv.size()<=1) return; // ignore 0 point crossings (0)

    auto blank_img = larocv::BlankImage(img,0);
    
    //
    // detect edges
    //
    const auto& xs_v = xs_vv.at(1);
    for(const auto xs : xs_v) {
      if (!xs) throw std::runtime_error("invalid point");
      blank_img.at<uchar>(xs->y,xs->x) = (uchar)255;
    }

    auto edge_ctor_v = larocv::FindContours(blank_img);

    for(const auto& edge_ctor : edge_ctor_v) {
      auto nz_edge_pt_v = larocv::FindNonZero(larocv::MaskImage(blank_img,edge_ctor,-1,false));
      if ((int)nz_edge_pt_v.size() < edge_size) continue;
      float mx,my;
      mx=my=0.0;
      for(const auto& nz_edge_pt : nz_edge_pt_v) {
	mx += nz_edge_pt.x;
	my += nz_edge_pt.y;
      }
      mx /= (float)nz_edge_pt_v.size();
      my /= (float)nz_edge_pt_v.size();
      geo2d::Vector<float> mid_pt(mx,my);
      geo2d::Vector<float> edge_pt;

      float min_dist = larocv::kINVALID_FLOAT;
      for(const auto& nz_edge_pt : nz_edge_pt_v) {
	float distance = geo2d::dist(geo2d::Vector<float>(nz_edge_pt.x,nz_edge_pt.y),mid_pt);
	if (distance < min_dist) {
	  min_dist = distance;
	  edge_pt.x = nz_edge_pt.x;
	  edge_pt.y = nz_edge_pt.y;
	}
      } // end edge point

      _edge_pt_v.emplace_back(edge_pt.x,edge_pt.y);
    } // end this grouping
    
    if (xs_vv.size() <= 3) return; // ignore 2 point or less crossing (0,1,2)

    blank_img.setTo(0);

    //
    // detect branches
    //
    for(size_t bid=3; bid < xs_vv.size(); ++bid) {
      const auto& xs_v = xs_vv.at(bid);
      for(const auto xs : xs_v) {
	if (!xs) throw std::runtime_error("invalid point");
	blank_img.at<uchar>(xs->y,xs->x) = (uchar)255;
      }
    }

    auto branch_ctor_v = larocv::FindContours(blank_img);
    for(const auto& branch_ctor : branch_ctor_v) {
      auto nz_branch_pt_v = larocv::FindNonZero(larocv::MaskImage(blank_img,branch_ctor,-1,false));
      if ((int)nz_branch_pt_v.size() < branch_size) continue;
      float mx,my;
      mx=my=0.0;
      for(const auto& nz_branch_pt : nz_branch_pt_v) {
	mx += nz_branch_pt.x;
	my += nz_branch_pt.y;
      }
      mx /= (float)nz_branch_pt_v.size();
      my /= (float)nz_branch_pt_v.size();
      geo2d::Vector<float> mid_pt(mx,my);
      geo2d::Vector<float> branch_pt; 
      float min_dist = larocv::kINVALID_FLOAT;
      for(const auto& nz_branch_pt : nz_branch_pt_v) {
	float distance = geo2d::dist(geo2d::Vector<float>(nz_branch_pt.x,nz_branch_pt.y),mid_pt);
	if (distance < min_dist) {
	  min_dist = distance;
	  branch_pt.x = nz_branch_pt.x;
	  branch_pt.y = nz_branch_pt.y;
	}
      } // end branch
      
      _branch_pt_v.emplace_back(branch_pt.x,branch_pt.y);
    } // end this grouping
    
    return;
  }

  float Polygon::Fraction(const cv::Mat& img1, const cv::Mat& img2) const {
    float res = 0.0;
    
    float top = larocv::CountNonZero(larocv::MaskImage(img1,_ctor,-1,false));
    float bot = larocv::CountNonZero(larocv::MaskImage(img2,_ctor,-1,false));

    if (bot != 0)
      res = top / bot;

    return res;
  }

}

#endif
