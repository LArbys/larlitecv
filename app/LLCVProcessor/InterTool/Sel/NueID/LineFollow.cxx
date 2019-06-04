#ifndef __LINEFOLLOW_CXX__
#define __LINEFOLLOW_CXX__

#include "LineFollow.h"

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <sstream>

#include <stdexcept>

namespace llcv {
  
  LineFollow::LineFollow() : llcv_base("LineFollow") {
    this->set_verbosity((msg::Level_t)2);
    _thickness = 2;
    _radius = 15;
    _radius_v.clear();
    _radius_v.resize(7);
    _radius_v[0] = 8;
    _radius_v[1] = 10;
    _radius_v[2] = 12;
    _radius_v[3] = 14;
    _radius_v[4] = 16;
    _radius_v[5] = 18;
    _radius_v[6] = 20;
    return;
  }
  
  void LineFollow::SetImageDimension(const cv::Mat& img, const cv::Mat& dead) {
    _img       = img.clone();
    _dead_img  = dead.clone();
    _black_img = larocv::BlankImage(img,0);
    _white_img = larocv::BlankImage(img,255);
    _bond_img  = _DeadWirePatch.WireBondage(_img,_dead_img);
  }

  larocv::GEO2D_ContourArray_t LineFollow::FollowEdgeLine(const geo2d::Vector<float>& start) {
    LLCV_DEBUG() << "start" << std::endl;
    larocv::GEO2D_ContourArray_t ret_v;

    float min_radius = 10;
    float max_radius = 10;

    //
    // find the first point off the boundary
    //

    geo2d::Vector<float> init_pt;
    auto first_exists = InitializeFirstPoint(start,init_pt);
    if (!first_exists) return ret_v;

    ret_v.emplace_back(AsLine(start,init_pt));
    
    //
    // Move along the line
    //

    float radius     = _radius;
    float min_angle  = 0;
    float _angle_tol = 20;
    float _dist_tol  = 10;

    geo2d::Vector<float> center_pt = init_pt;
    geo2d::Vector<float> prev_pt;
    geo2d::Circle<float> circle;
    int ctr = 0;

    bool first = true;
    std::vector<float> radius_v;

    while(1) {
      ctr+=1;

      //
      // get distance to edge
      //
      float edge_dist = DistanceToEdge(center_pt);

      if (edge_dist<=1) {
	LLCV_DEBUG() << "@edge break" << std::endl;
	break;
      }

      radius_v.clear();
      radius_v.reserve(_radius_v.size());

      for(auto rad : _radius_v) {
	if (edge_dist < rad)
	  radius_v.push_back(edge_dist - 1);
	else
	  radius_v.push_back(rad);
	LLCV_DEBUG() << "@rad=" << radius_v.back() << std::endl;
      }

      
      circle.center = center_pt;

      auto pt_vv = larocv::OnCircleGroupsOnCircleArray(_bond_img,circle.center,radius_v);

      //
      // choose the comparison point
      //
      
      if (first) {
	prev_pt = geo2d::mean(start,circle.center);
	first = false;
      } 

      float mmin_angle = larocv::kINVALID_FLOAT;
      size_t mmpid1 = larocv::kINVALID_SIZE;
      geo2d::Vector<float> far_pt;
      std::vector<float> ignore_v;

      geo2d::Line<float> line2(circle.center, prev_pt - circle.center);

      _PiRange.SetAngle(Angle(circle.center,prev_pt),45,45);
      
      LLCV_DEBUG() << "start=(" << start.x << "," << start.y << ")" << std::endl;
      LLCV_DEBUG() << "center=(" << circle.center.x << "," << circle.center.y << ")" << std::endl;
      LLCV_DEBUG() << "prev_pt=(" << prev_pt.x << "," << prev_pt.y << ")" << std::endl;
      LLCV_DEBUG() << "line2 angle=" << Angle(circle.center,prev_pt) << std::endl;

      for(size_t rid=0; rid < radius_v.size(); ++rid) {
	const auto& pt_v = pt_vv[rid];

	LLCV_DEBUG() << "@rid=" << rid << " rad=" << radius_v[rid] << " pt_v sz=" << pt_v.size() << std::endl;
	
	if (pt_v.size()<2) continue;
	
	ignore_v.clear();
	ignore_v.resize(pt_v.size(),false);
	
	for(size_t pid=0; pid<pt_v.size(); ++pid) {
	  const auto& pt = pt_v[pid];
	  geo2d::Line<float> line(circle.center, pt - circle.center);
	  auto angle = Angle(circle.center,pt);
	  LLCV_DEBUG() << "@pid= " << pid << " (" << pt.x << "," << pt.y << ") a=" << angle << std::endl;
	  if (_PiRange.Inside(angle)) {
	    LLCV_DEBUG() << "...inside" << std::endl;
	    ignore_v.at(pid) = true;
	  }
	}
	
	//
	// minimize angle between crossing points (closest to 180)
	//
	min_angle = larocv::kINVALID_FLOAT;
	
	size_t mpid1 = larocv::kINVALID_SIZE;
	
	for(size_t pid=0; pid < pt_v.size(); ++pid) {
	  if (ignore_v[pid]) continue;
	  const auto& pt1 = pt_v[pid];
	  auto line1 = geo2d::Line<float>(pt1, pt1 - circle.center);
	  float angle = std::fabs(geo2d::angle(line1) - geo2d::angle(line2));
	  if (angle > 90) angle = std::fabs(180 - angle);
	  LLCV_DEBUG() << "(" << pid << ") (" << pt1.x << "," << pt1.y << ")&(" << prev_pt.x << "," << prev_pt.y << ") a=" << angle << std::endl;
	  if (angle < min_angle) {
	    min_angle = angle;
	    mpid1 = pid;
	    LLCV_DEBUG() << "...accepted (" << min_angle << "," << mpid1 << ")" << std::endl;
	  }
	}

	if (min_angle < mmin_angle) {
	  mmin_angle = min_angle;
	  far_pt     = pt_v[mpid1];
	}

      }

      if (mmin_angle > _angle_tol) {
	LLCV_DEBUG() << "exit angle" << std::endl;
	break;
      }

      ret_v.emplace_back(AsLine(circle.center,far_pt));

      prev_pt   = circle.center;
      center_pt = far_pt;
      LLCV_DEBUG() << "next ctr=" << ctr << std::endl;

      if (ctr>100) break;
    }
    
    LLCV_DEBUG() << "end ctr=" << ctr << std::endl;
    return ret_v;
  }

  bool LineFollow::InitializeFirstPoint(geo2d::Vector<float> start, geo2d::Vector<float>& init_pt) {
    bool exists = false;

    // define the small image
    float pad_x = 20;
    float pad_y = 20;
    
    float min_x = start.x;
    float min_y = start.y;
    
    min_x -= pad_x;
    min_y -= pad_y;
    
    float dx = 2*pad_x;
    float dy = 2*pad_y;

    LLCV_DEBUG() << "start=(" << start.x << "," << start.y << ") : m=(" << min_x << "," << min_y << ") : d=(" << dx << "," << dy << ")" << std::endl;
    auto small_img = larocv::SmallImg(_img,geo2d::Vector<float>(min_x,min_y),dx,dy);

    float ddx = min_x-dx;
    float ddy = min_y-dy;

    int ecase = -1;
    bool repo1 = false;
    bool repo2 = false;
    bool repo3 = false;

    if (start.x  < (_img.cols-1) and start.y == 0) {
      
      if (start.x > (_img.rows - 1 - dx/2.0)) {
      	start.x *= -1;
      	start.x += (_img.rows - 1);
	start.x *= -1;
	start.x += 2*dx;
      	repo2 = true;
      }

      if (!repo2 and start.x > 1.5*dx){
	start.x -= ddx;
	repo1 = true;
      }
      
      ecase = 0;
    }
    else if (start.x < (_img.cols-1) and start.y >= (_img.rows-1)) {

      if (start.x > (_img.rows - 1 - dx/2.0)) {
      	start.x *= -1;
      	start.x += (_img.rows - 1);
	start.x *= -1;
	start.x += 2*dx;
      	repo2 = true;
      }

      if (!repo2 and start.x > 1.5*dx) {
	start.x -= ddx;
	repo1 = true;
      }

      start.y  = small_img.rows-1;

      ecase = 1;
    }
    else if (start.x == 0 and start.y < (_img.rows-1)) {

      if (start.y > (_img.rows - 1 - dy/2.0)) {
	start.y *= -1;
	start.y += (_img.rows - 1);
	start.y *= -1;
	start.y += 2*dy;
	repo2 = true;
      }

      if (!repo2 and start.y > 1.5*dy) {
	start.y -= ddy;
	repo1 = true;
      }

      ecase = 2;
    }
    else if (start.x >= (_img.cols-1) and start.y < (_img.rows-1)) {

      start.x  = small_img.cols-1;

      if (start.y > (_img.rows - 1 - dy/2.0)) {
	start.y *= -1;
	start.y += (_img.rows - 1);
	start.y *= -1;
	start.y += 2*dy;
	repo2 = true;
      }
      
      if (!repo2 and start.y > 1.5*dy) {
	start.y -= ddy;
	repo1 = true;
      }

      ecase = 3;
    }

    
    else {
      return false;
      throw std::runtime_error("unhandled case");
    }

    LLCV_DEBUG() << "ecase=" << ecase << " repo1=" << repo1 << " repo2=" << repo2 << std::endl;

    // define the circle
    LLCV_DEBUG() << "ddx: " << ddx << " ddy: " << ddy << std::endl;
    LLCV_DEBUG() << "small start=(" << start.x << "," << start.y << ")" << std::endl;
    auto small_cimg = larocv::BlankImage(small_img,0);
    cv::circle(small_cimg,cv::Point((int)(start.x+0.5),(int)(start.y+0.5)),(int)20,cv::Scalar(255),1);

    auto nz_pt_v = larocv::FindNonZero(small_cimg);
    
    auto small_bimg = larocv::BlankImage(small_cimg,0);
    auto small_wimg = larocv::BlankImage(small_cimg,255);

    float thickness = 3;

    float n_line = -1.0*larocv::kINVALID_FLOAT;

    geo2d::Vector<float> line_pt;

    for(auto& nz_pt : nz_pt_v) {
      
      small_bimg.setTo(cv::Scalar(0));
      small_wimg.setTo(cv::Scalar(255));

      cv::line(small_bimg,
	       cv::Point((int)(start.x + 0.5),(int)(start.y+0.5)),
	       cv::Point((int)(nz_pt.x + 0.5),(int)(nz_pt.y+0.5)),
	       cv::Scalar(255),
	       thickness);
      
      auto line_ctor_v = larocv::FindContours(small_bimg);
      const auto& line_ctor = line_ctor_v.front();
      
      auto line_mask = larocv::MaskImage(small_img,line_ctor,-1,false);

      float n_line_mask = larocv::CountNonZero(line_mask);
      
      if (n_line_mask > n_line) {
	n_line = n_line_mask;
      	line_pt = nz_pt;
      }
      
    }

    LLCV_DEBUG() << "line_pt=(" << line_pt.x << "," << line_pt.y << ") n_line=" << n_line << std::endl;
    
    init_pt = line_pt;

    switch(ecase) {

    case 0 : {
      if (repo1)
	init_pt.x += (ddx);

      if (repo2) {
	init_pt.x -= 2*dx;
	init_pt.x *= -1;
      	init_pt.x -= (_img.rows - 1);
      	init_pt.x *= -1;
      }

      break;
    }

    case 1: {
      if (repo1)
	init_pt.x += (ddx);

      if (repo2) {
	init_pt.x -= 2*dx;
	init_pt.x *= -1;
      	init_pt.x -= (_img.rows - 1);
      	init_pt.x *= -1;
      }

      init_pt.y  = _img.rows - (small_img.rows - init_pt.y) - 1;
      break;
    }

    case 2: {
      if (repo1)
	init_pt.y += (ddy);
      if (repo2) {
	init_pt.y -= 2*dy;
	init_pt.y *= -1;
	init_pt.y -= (_img.rows - 1);
	init_pt.y *= -1;
      }
      break;
    }

    case 3: {
      init_pt.x  = _img.cols - (small_img.cols - init_pt.x) - 1;
      if (repo1)
	init_pt.y += (ddy);
      if (repo2) {
	init_pt.y -= 2*dy;
	init_pt.y *= -1;
	init_pt.y -= (_img.rows - 1);
	init_pt.y *= -1;
      }
      break;
    }
      
    default:  {
      break;
    }
      
    }
    
    LLCV_DEBUG() << "init_pt=(" << init_pt.x << "," << init_pt.y << ")" << std::endl;
    
    exists = true;
    return exists;
  }

  larocv::GEO2D_Contour_t LineFollow::AsLine(geo2d::Vector<float> pt1, geo2d::Vector<float> pt2) {
    larocv::GEO2D_Contour_t ret;
    
    float pad_x = 30;
    float pad_y = 30;

    float min_x = std::min(pt1.x,pt2.x);
    float min_y = std::min(pt1.y,pt2.y);

    float dx = std::abs(pt1.x - pt2.x)/2.0;
    float dy = std::abs(pt1.y - pt2.y)/2.0;

    dx += pad_x;
    dy += pad_y;

    auto small_mat = cv::Mat((int)(dx+0.5),(int)(dy+0.5),_black_img.type(),cv::Scalar(0));

    float ox = min_x - pad_x/5;
    float oy = min_y - pad_y/5;

    pt1.x -= ox;
    pt1.y -= oy;
    
    pt2.x -= ox;
    pt2.y -= oy;

    float thickness = 3;
    cv::line(small_mat,
    	     cv::Point((int)(pt1.x+0.5),(int)(pt1.y+0.5)),
	     cv::Point((int)(pt2.x+0.5),(int)(pt2.y+0.5)),
	     cv::Scalar(255),
	     thickness);

    auto line_ctor_v = larocv::FindContours(small_mat);

    if (!line_ctor_v.empty()) {
      const auto& line_ctor = line_ctor_v.front();
      for(size_t lid=0; lid<line_ctor.size(); ++lid) {
	geo2d::Vector<int> cpt;
	cpt.x = line_ctor[lid].x;
	cpt.y = line_ctor[lid].y;
	
	cpt.x += ox;
	cpt.y += oy;

	ret.emplace_back(std::move(cpt));
      }
    }

    return ret;
  }
  
  larocv::GEO2D_Contour_t LineFollow::EdgePoints() {
    larocv::GEO2D_Contour_t ret;

    _black_img.setTo(0);
    
    // bottom
    cv::line(_black_img,cv::Point(0,0),cv::Point(_black_img.cols-1,0),cv::Scalar(255),1);

    // left
    cv::line(_black_img,cv::Point(0,0),cv::Point(0,_black_img.rows-1),cv::Scalar(255),1);

    // right
    cv::line(_black_img,cv::Point(_black_img.cols-1,0),cv::Point(_black_img.cols-1,_black_img.rows-1),cv::Scalar(255),1);
    
    // top
    cv::line(_black_img,cv::Point(0,_black_img.rows-1),cv::Point(_black_img.cols-1,_black_img.rows-1),cv::Scalar(255),1);

    _white_img.setTo(0);
    _img.copyTo(_white_img,_black_img);
    
    //
    // find edge points
    //
    
    auto ctor_v = larocv::FindContours(_white_img);
    ret.reserve(ctor_v.size());
    for(const auto& ctor : ctor_v) {
      auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(_white_img,ctor,-1,false));
      // calc the average
      float mean_x = 0;
      float mean_y = 0;
      for(const auto& nz_pt : nz_pt_v) {
	mean_x += nz_pt.x;
	mean_y += nz_pt.y;
      }
      mean_x /= (float)nz_pt_v.size();
      mean_y /= (float)nz_pt_v.size();

      ret.emplace_back((int)mean_x+0.5,(int)mean_y+0.5);
    }

    //
    // determine if dead wires exists on image boundary
    //

    //
    // check the bottom
    //
    int bgood = (int)_dead_img.at<uchar>(0,0);
    int bot_max = 10;
    if (!bgood) {
      _black_img.setTo(0);
      _white_img.setTo(0);
      int rowid = 0;
      for(int rid=0; rid<bot_max; ++rid) {
	uint d = (uint)_dead_img.at<uchar>(rid,0);
	rowid = rid;
	if (d) break;
      }

      LLCV_DEBUG() << "bottom in dead region rowid=" << rowid << std::endl;
      
      cv::line(_black_img,cv::Point(0,rowid),cv::Point(_black_img.cols-1,rowid),cv::Scalar(255),1);
      _img.copyTo(_white_img,_black_img);

      auto ctor_v = larocv::FindContours(_white_img);
      ret.reserve(ctor_v.size());
      for(const auto& ctor : ctor_v) {
	auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(_white_img,ctor,-1,false));
	float mean_x = 0;
	float mean_y = 0;
	for(const auto& nz_pt : nz_pt_v)
	  mean_x += nz_pt.x;

	mean_x /= (float)nz_pt_v.size();
	ret.emplace_back((int)mean_x+0.5,mean_y);

      }
    } // end bottom search

    //
    // check the top
    //
    int tgood = (int)_dead_img.at<uchar>(_black_img.rows-1,0);
    int top_max = 10;
    if (!tgood) {
      _black_img.setTo(0);
      _white_img.setTo(0);
      int rowid = _black_img.rows-1;
      for(int rid=(_black_img.rows - 1); rid>((_black_img.rows - 1) - top_max); --rid) {
	uint d = (uint)_dead_img.at<uchar>(rid,0);
	rowid = rid;
	if (d) break;
      }

      LLCV_DEBUG() << "top in dead region rowid=" << rowid << std::endl;
      
      cv::line(_black_img,cv::Point(0,rowid),cv::Point(_black_img.cols-1,rowid),cv::Scalar(255),1);
      _img.copyTo(_white_img,_black_img);

      auto ctor_v = larocv::FindContours(_white_img);
      ret.reserve(ctor_v.size());
      for(const auto& ctor : ctor_v) {
	auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(_white_img,ctor,-1,false));
	float mean_x = 0;
	float mean_y = _black_img.rows-1;
	for(const auto& nz_pt : nz_pt_v)
	  mean_x += nz_pt.x;

	mean_x /= (float)nz_pt_v.size();
	ret.emplace_back((int)mean_x+0.5,(int)mean_y);
      }
    } // end top search


    return ret;
  }

  float LineFollow::DistanceToEdge(const geo2d::Vector<float>& pt) const {
    float ret = larocv::kINVALID_FLOAT;
    
    ret = std::min(ret,(float)pt.x);
    ret = std::min(ret,(float)pt.y);
    ret = std::min(ret,(float)_img.rows - (float)pt.x - 1);
    ret = std::min(ret,(float)_img.cols - (float)pt.y - 1);
    
    return ret;
  }

  int LineFollow::Quadrant(const geo2d::Vector<float>& pt) const{
    int ret = larocv::kINVALID_INT;
    
    if (pt.x >= 0 and pt.y >= 0)
      ret = 0;
    else if (pt.x <= 0 and pt.y >= 0)
      ret = 1;
    else if (pt.x <= 0 and pt.y <= 0)
      ret = 2;
    else if (pt.x >=0 and pt.y <= 0)
      ret = 3;

    if (ret == larocv::kINVALID_INT)
      throw std::runtime_error("unhandled quadrant");

    return ret;
  }

  float LineFollow::Angle(const geo2d::Vector<float>& origin, const geo2d::Vector<float>& pt1) const {
    float ret = larocv::kINVALID_FLOAT;

    auto pt = pt1 - origin;

    float cos = pt.x;
    cos /= std::sqrt(pt.x*pt.x + pt.y*pt.y);
    cos = std::acos(cos);
    cos *= 180.0/3.14159;
    
    switch (Quadrant(pt)) {
    case 0 : { break; }
    case 1 : { break; }
    case 2 : { cos = 360 - cos; break; }
    case 3 : { cos = 360 - cos; break; }
    default : {
      throw std::runtime_error("unhandled quadrant");
      break;
    }
    }
    ret = cos;
    
    return ret;
  }


}

#endif


