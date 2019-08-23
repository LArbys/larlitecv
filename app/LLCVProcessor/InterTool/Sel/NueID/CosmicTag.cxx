#ifndef __COSMICTAG_CXX__
#define __COSMICTAG_CXX__

#include "CosmicTag.h"

#ifndef __CLING__
#ifndef __CINT__
#include <opencv2/core.hpp>
#include "opencv2/imgproc.hpp"
#endif
#endif

namespace llcv {

  CosmicTag::CosmicTag()  : llcv_base("CosmicTag") {
    this->set_verbosity((msg::Level_t)0);
    return;
  }
  
  void CosmicTag::Reset() {
    _cosmic_ctor_v.clear();
    _cosmic_line_ctor_vv.clear();
    _cosmic_end_v.clear();
    return;
  }

  void CosmicTag::TagCosmic(const cv::Mat& img, const cv::Mat& dimg) {
    _LineFollow.SetImageDimension(img,dimg);
	
    auto edge_v = _LineFollow.EdgePoints();

    LLCV_DEBUG() << "got edge sz=" << edge_v.size() << std::endl;
    _cosmic_line_ctor_vv.resize(edge_v.size());
    _cosmic_end_v.reserve(edge_v.size());
    _cosmic_ctor_v.reserve(edge_v.size());

    for(size_t eid=0; eid<edge_v.size(); ++eid) {
      LLCV_DEBUG() << "@eid=" << eid << std::endl;

      const auto& edge = edge_v[eid];
      auto& cosmic_line_ctor_v = _cosmic_line_ctor_vv[eid];
	
      cosmic_line_ctor_v = _LineFollow.FollowEdgeLine(geo2d::Vector<float>(edge.x,edge.y));

      if (cosmic_line_ctor_v.empty()) continue; 

      auto cosmic_ctor = LinesToContour(cosmic_line_ctor_v);

      geo2d::Vector<float> cosmic_end_pt(-1,-1);
      larocv::FindEdge(cosmic_ctor,edge,cosmic_end_pt);

      _cosmic_end_v.emplace_back(std::move(cosmic_end_pt));
      _cosmic_ctor_v.emplace_back(std::move(cosmic_ctor));
    }
      
  }
  
  larocv::GEO2D_Contour_t CosmicTag::LinesToContour(const larocv::GEO2D_ContourArray_t& ctor_v) {
    larocv::GEO2D_Contour_t ret;
    
    auto& white_img = _LineFollow.WhiteImage();
    auto& black_img = _LineFollow.BlackImage();

    white_img.setTo(255);
    black_img.setTo(0);

    for(const auto& ctor : ctor_v) 
      white_img = larocv::MaskImage(white_img,ctor,-1,true);

    cv::bitwise_not(white_img,black_img);

    auto line_ctor_v = larocv::FindContours(black_img);

    ret = line_ctor_v.front();
    
    return ret;
  }

  void CosmicTag::DrawContours(cv::Mat& mat3d) const {
    // draw the green lines
    for(size_t lid=0; lid < _cosmic_ctor_v.size(); ++lid) {
      const auto& cosmic_ctor = _cosmic_ctor_v[lid];
      cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,cosmic_ctor),-1,cv::Scalar(0,255,0));
    }
  }

  void CosmicTag::DrawLines(cv::Mat& mat3d) const {
    // draw the red lines
    for(size_t eid=0; eid < _cosmic_line_ctor_vv.size(); ++eid) {
      const auto& cosmic_line_ctor_v = _cosmic_line_ctor_vv[eid];
      for(size_t lid=0; lid < cosmic_line_ctor_v.size(); ++lid) {
	const auto cosmic_line_ctor = cosmic_line_ctor_v[lid];
	if (lid == 0)
	  cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,cosmic_line_ctor),-1,cv::Scalar(255,0,255));
	else
	  cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,cosmic_line_ctor),-1,cv::Scalar(0,0,255));
      }
    }
  }
  
  size_t CosmicTag::NearestCosmicToPoint(const cv::Point_<int>& pt, float& distance) const {
    size_t ret = larocv::kINVALID_SIZE;
    distance = larocv::kINVALID_FLOAT;
    
    geo2d::Vector<float> pt_f(pt.x,pt.y);

    for(size_t cid=0; cid<_cosmic_ctor_v.size(); ++cid) {
      const auto& cosmic_ctor = _cosmic_ctor_v[cid];
      auto dist = larocv::Pt2PtDistance(pt_f,cosmic_ctor);
      if (dist < distance) {
	distance = dist;
	ret = cid;
      }
    }
    
    return ret;
  }

  geo2d::Vector<float> CosmicTag::NearestCosmicEndToPoint(const cv::Point_<int>& pt, float& distance, size_t& cosmic_id) const {

    geo2d::Vector<float> ret_pt(larocv::kINVALID_FLOAT,larocv::kINVALID_FLOAT);

    size_t ret = larocv::kINVALID_SIZE;
    distance = larocv::kINVALID_FLOAT;
    
    geo2d::Vector<float> pt_f(pt.x,pt.y);

    for(size_t cid=0; cid<_cosmic_end_v.size(); ++cid) {
      const auto& cosmic_end = _cosmic_end_v[cid];
      auto dist = geo2d::dist(pt_f,cosmic_end);
      if (dist < distance) {
	distance = dist;
	ret = cid;
      }
    }
    
    cosmic_id = ret;

    if (ret == larocv::kINVALID_SIZE)
      return ret_pt;

    ret_pt = _cosmic_end_v[ret];
    
    return ret_pt;
  }
    
  size_t CosmicTag::NearestCosmicToContour(const larocv::GEO2D_Contour_t& ctor, float& distance) const {
    size_t ret = larocv::kINVALID_SIZE;
    distance = larocv::kINVALID_FLOAT;
      
    for(size_t cid=0; cid<_cosmic_ctor_v.size(); ++cid) {
      const auto& cosmic_ctor = _cosmic_ctor_v[cid];
      auto dist = larocv::Pt2PtDistance(cosmic_ctor,ctor);
      if (dist < distance) {
	distance = dist;
	ret = cid;
      }
    }
      
    return ret;
  }


  geo2d::Vector<float> CosmicTag::NearestCosmicEndToContour(const larocv::GEO2D_Contour_t& ctor, float& distance, size_t& cosmic_id) const {
    
    geo2d::Vector<float> ret_pt(larocv::kINVALID_FLOAT,larocv::kINVALID_FLOAT);
    size_t ret = larocv::kINVALID_SIZE;
    distance = larocv::kINVALID_FLOAT;
      
    for(size_t cid=0; cid<_cosmic_end_v.size(); ++cid) {
      const auto& cosmic_end = _cosmic_end_v[cid];
      auto dist = larocv::Pt2PtDistance(cosmic_end,ctor);
      if (dist < distance) {
	distance = dist;
	ret = cid;
      }
    }
    
    cosmic_id = ret;

    if (ret == larocv::kINVALID_SIZE)
      return ret_pt;
    
    ret_pt = _cosmic_end_v[ret];

    return ret_pt;
  }

}

#endif 
