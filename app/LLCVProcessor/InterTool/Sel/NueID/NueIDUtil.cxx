#ifndef __NUEIDUTIL_CXX__
#define __NUEIDUTIL_CXX__

#include "NueIDUtil.h"

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

namespace llcv {

  size_t FindClosestContour(const larocv::GEO2D_ContourArray_t& ctor_v,
			    const geo2d::Vector<float>& pt,
			    float& distance) {

    size_t res;
    
    // get the contour closest to the vertex
    float dist = kINVALID_FLOAT;
    size_t close_id = kINVALID_SIZE;
    for(size_t cid=0; cid<ctor_v.size(); ++cid) {
      const auto& ctor = ctor_v[cid];
      auto dist_tmp = larocv::Pt2PtDistance(pt,ctor);
      if (dist_tmp < dist) {
	dist = dist_tmp;
	close_id = cid;
      }
    }
    
    distance = dist;
    res = close_id;
    return res;
  }
  

  larocv::GEO2D_Contour_t MaximizeLine(const cv::Mat& timg3d_par,
				       const Triangle& triangle,
				       float& nline_pixels,
				       geo2d::Vector<float>& edge,
				       cv::Mat& white_img) {
    larocv::GEO2D_Contour_t res;

    nline_pixels = -1*larocv::kINVALID_FLOAT;

    const geo2d::Vector<float>* base_left  = nullptr;
    const geo2d::Vector<float>* base_right = nullptr;

    const geo2d::Vector<float>* base_up   = nullptr;
    const geo2d::Vector<float>* base_down = nullptr;

    bool isinf = false;
    
    if (triangle.Base1().x < triangle.Base2().x) {
      base_left  = &triangle.Base1();
      base_right = &triangle.Base2();
    }
    else if(triangle.Base1().x > triangle.Base2().x) {
      base_left  = &triangle.Base2();
      base_right = &triangle.Base1();
    }
    // it's infinity
    else {
      isinf = true;
      if (triangle.Base1().y < triangle.Base2().y) {
	base_down = &triangle.Base1();
	base_up   = &triangle.Base2();
      }
      else {
	base_down = &triangle.Base2();
	base_up   = &triangle.Base1();
      }
    }

    geo2d::LineSegment<float> base_seg;
    
    float lo;
    float hi; 

    if (!isinf)  {
      base_seg = geo2d::LineSegment<float>(*base_left,*base_right);
      lo = base_left->x;
      hi = base_right->x;
    }
    else {
      base_seg = geo2d::LineSegment<float>(*base_down,*base_up);
      lo = base_down->y;
      hi = base_up->y;
    }

    geo2d::Vector<float> pt;
    
    float step = 0.1;
    for(float rid=lo; rid<hi; rid+=step) {
      
      white_img.setTo(cv::Scalar(0));
      
      if (!isinf) {
	pt.x = rid;
	pt.y = -1;
	try {
	  pt.y = base_seg.y(pt.x);
	}
	catch (...) {
	  continue;
	}
      } 
      
      else {
	pt.x = base_down->x;
	pt.y = rid;
      }

      // draw a line
      cv::line(white_img,triangle.Apex(),pt,cv::Scalar(255),3);
      
      // find the contour
      auto line_ctor_v = larocv::FindContours(white_img);
      if (line_ctor_v.empty()) continue;
      const auto& line_ctor = line_ctor_v.front();
      
      if (line_ctor.empty()) continue;
      
      auto timg3d_par_line = larocv::MaskImage(timg3d_par,line_ctor,-1,false);
      float nline_pixels_tmp = (float)larocv::CountNonZero(timg3d_par_line);
      
      if (nline_pixels_tmp > nline_pixels) {
	nline_pixels = nline_pixels_tmp;
	edge = pt;
	res = line_ctor;
      }
      
    }

    return res;
    
  }


  larocv::GEO2D_Contour_t MaximizePolygonLine(const cv::Mat& timg3d_mask,
					      const std::vector<Polygon>& polygon_v,
					      Triangle& triangle,
					      float& nline_pixels,
					      float& npar_pixels,
					      geo2d::Vector<float>& edge,
					      cv::Mat& white_img) {
    auto timg3d_par = timg3d_mask.clone();
    timg3d_par.setTo(cv::Scalar(0));

    for(const auto& polygon : polygon_v)
      timg3d_par += larocv::MaskImage(timg3d_mask,polygon.Contour(),-1,false);
    
    npar_pixels = (float)larocv::CountNonZero(timg3d_par);

    size_t size_estimate = 0;
    for(const auto& polygon : polygon_v)
      size_estimate += polygon.Contour().size();

    larocv::GEO2D_Contour_t comb_ctor_v;
    comb_ctor_v.reserve(size_estimate);

    for(const auto& polygon : polygon_v)
      for(const auto& pt : polygon.Contour())
	comb_ctor_v.emplace_back(pt);

    triangle = Triangle(comb_ctor_v,triangle.Apex());

    return MaximizeLine(timg3d_par, triangle, nline_pixels, edge, white_img);
  }
  

  larocv::GEO2D_Contour_t MaximizeTriangleLine(const cv::Mat& timg3d_mask,
					       const Triangle& triangle,
					       float& nline_pixels,
					       float& npar_pixels,
					       geo2d::Vector<float>& edge,
					       cv::Mat& white_img) {
    

    auto timg3d_par = larocv::MaskImage(timg3d_mask,triangle.Contour(),-1,false);
    
    npar_pixels = (float)larocv::CountNonZero(timg3d_par);
    
    return MaximizeLine(timg3d_par,triangle,nline_pixels,edge,white_img);
  }


  float MinimizeToEdge(const larocv::GEO2D_Contour_t& ctor,
		       size_t cropx,
		       size_t cropy) {

    float res = larocv::kINVALID_FLOAT;
    for(const auto& pt : ctor) {
      res = std::min(res,(float)pt.x);
      res = std::min(res,(float)pt.y);
      res = std::min(res,(float)cropx - (float)pt.x);
      res = std::min(res,(float)cropy - (float)pt.y);
    }
    return res;
  }

  float NearestPolygonToCosmic(const std::vector<Polygon>& polygon_v,
			       const std::vector<CosmicTag>& CosmicTag_v,
			       size_t plane) {
		
    float ret = larocv::kINVALID_FLOAT;
    const auto& CosmicTag_ = CosmicTag_v[plane];
    for(const auto& polygon : polygon_v) {
      float distance = larocv::kINVALID_FLOAT;
      auto id = CosmicTag_.NearestCosmicToContour(polygon.Contour(),distance);
      ret = std::min(ret,distance);
    }
    return ret;
  }

  float NearestPolygonToCosmicEnd(const std::vector<Polygon>& polygon_v,
				  const std::vector<CosmicTag>& CosmicTag_v,
				  size_t plane) {
    
    float ret = larocv::kINVALID_FLOAT;
    const auto& CosmicTag_ = CosmicTag_v[plane];
    for(const auto& polygon : polygon_v) {
      float distance = larocv::kINVALID_FLOAT;
      size_t cosmic_id = larocv::kINVALID_SIZE;
      auto id = CosmicTag_.NearestCosmicEndToContour(polygon.Contour(),distance,cosmic_id);
      ret = std::min(ret,distance);
    }
    return ret;
  }

  float PointCosmicDistance(const geo2d::Vector<float>& pt,
			    const std::vector<CosmicTag>& CosmicTag_v,
			    const size_t plane) {

    float ret = larocv::kINVALID_FLOAT;
    const auto& CosmicTag_ = CosmicTag_v[plane];
    geo2d::Vector<int> pt_i((int)(pt.x+0.5),(int)(pt.y+0.5));
    auto id = CosmicTag_.NearestCosmicToPoint(pt_i,ret);
    return ret;
  }

  float PointCosmicEndDistance(const geo2d::Vector<float>& pt,
			       const std::vector<CosmicTag>& CosmicTag_v,
			       const size_t plane) {

    float ret = larocv::kINVALID_FLOAT;
    const auto& CosmicTag_ = CosmicTag_v[plane];
    geo2d::Vector<int> pt_i((int)(pt.x+0.5),(int)(pt.y+0.5));
    size_t cosmic_id = larocv::kINVALID_SIZE;
    auto id = CosmicTag_.NearestCosmicEndToPoint(pt_i,ret,cosmic_id);
    return ret;
  }

  int DetectBrem(Triangle& triangle, 
		 const larocv::GEO2D_ContourArray_t& other_ctor_v,
		 std::vector<size_t>& id_v,
		 std::vector<size_t>& inside_v,
		 float brem_dist,
		 float brem_size) {
    int res = 0;

    //
    // count the number of other contours inside or touching the expanded triangle
    //
    
    float fraction = (triangle.Height() + brem_dist) / triangle.Height();

    triangle.Extend(fraction);
    
    auto tri_ctor = triangle.AsContour();
    
    id_v.reserve(other_ctor_v.size());
    inside_v.reserve(other_ctor_v.size());

    for(size_t oid=0; oid < other_ctor_v.size(); ++oid) {
      const auto& other_ctor = other_ctor_v[oid];
      
      auto common_area = larocv::AreaOverlap(other_ctor,tri_ctor);
      if (common_area == 0) continue;
      
      int pixel_area = larocv::ContourPixelArea(other_ctor);
      if (pixel_area > brem_size) {
	res += 1;
	id_v.push_back(oid);
      }
      inside_v.push_back(oid);
    }
    
    return res;
  }
  
  std::array<float,3> EstimateDirection(const Object2DCollection& obj_col,
					cv::Mat& white_img,
					ContourScan& ContourScan_,
					ShowerTools& ShowerTools_) {

    std::array<float,3> ret_v;
    ret_v[0] = -1*larocv::kINVALID_FLOAT;
    ret_v[1] = -1*larocv::kINVALID_FLOAT;
    ret_v[2] = -1*larocv::kINVALID_FLOAT;

    ContourScan_.Reset();

    // register the lines for this object
    for(const auto& obj2d : obj_col) {
      white_img.setTo(cv::Scalar(255));
      ContourScan_.RegisterContour(white_img,obj2d.Line(),obj2d.Plane(),-1);
    }

    // scan 
    ContourScan_.Scan();
    
    // return voxels
    auto vox_v = ContourScan_.Voxelize(0.3,0.3,0.3);

    LLCV_DEBUG() << "vox_v sz=" << vox_v.size() << std::endl;

    // calculate PCA
    ret_v = ShowerTools_.ComputePCA(vox_v,obj_col);
 
    return ret_v;
  }

  void ComputeLineParameters(Object2DCollection& obj_col,
			     const std::array<cv::Mat,3>& cimg_v,
			     cv::Mat& white_img) {
      
    auto black_img = white_img.clone();
    
    for(auto& obj2d : obj_col) {

      white_img.setTo(cv::Scalar(255));
      black_img.setTo(cv::Scalar(0));      
      
      auto plane = obj2d.Plane();

      const auto& cimg = cimg_v[plane];

      // mask the polygon
      for(const auto& poly : obj2d._polygon_v)
	black_img += larocv::MaskImage(cimg,poly.Contour(),-1,false);

      // mask out the line
      black_img = larocv::MaskImage(black_img,obj2d.Line(),-1,true);

      // draw a line of thickness 1 from start to end
      white_img.setTo(cv::Scalar(0));
      cv::line(white_img,
	       geo2d::Vector<int>((int)(obj2d.Start().x+0.5),
				  (int)(obj2d.Start().y+0.5)),
	       geo2d::Vector<int>((int)(obj2d.Edge().x+0.5),
				  (int)(obj2d.Edge().y+0.5)),
	       cv::Scalar(255),
	       1);

      auto nz_line_pt_v = larocv::FindNonZero(white_img);
      auto nz_ctor_pt_v = larocv::FindNonZero(black_img);
      
      float mean_dist = 0;
      float max_dist = -1.0*larocv::kINVALID_FLOAT;
      
      geo2d::Vector<float> pt1_f, pt2_f;

      for(const auto& nz_ctor_pt : nz_ctor_pt_v) {
	pt1_f.x = nz_ctor_pt.x;
	pt1_f.y = nz_ctor_pt.y;
	
	float ctor_min_dist = larocv::kINVALID_FLOAT;

	for(const auto& nz_line_pt : nz_line_pt_v) {
	  pt2_f.x = nz_line_pt.x;
	  pt2_f.y = nz_line_pt.y;
	  
	  auto dist = geo2d::dist(pt1_f,pt2_f);
	  ctor_min_dist = std::min(dist,ctor_min_dist);
	}

	max_dist   = std::max(ctor_min_dist,max_dist);
	mean_dist += ctor_min_dist;
      }

      if (!nz_ctor_pt_v.empty())
	mean_dist /= ((float) nz_ctor_pt_v.size());

      obj2d._line_max_dist = max_dist;
      obj2d._line_mean_dist = mean_dist;
      
    } // end obj2d
    
    return;
  }
  
  
  void SplitLineParameters(Object2DCollection& obj_col,
			   const std::array<cv::Mat,3>& cimg_v,
			   cv::Mat& white_img) {
    
    for(auto& obj2d : obj_col) {
      white_img.setTo(cv::Scalar(255));
      auto plane = obj2d.Plane();
      const auto& cimg = cimg_v[plane];
      
      // get the perp direction
      auto dx = obj2d.LinedX();
      auto dy = obj2d.LinedY();
      geo2d::Vector<float> dir_perp(dy,-1.0*dx);
      
      // get the middle of the line
      auto mid_pt = obj2d.Edge() + obj2d.Start();
      mid_pt /= (float)2.0;

      // make the line
      geo2d::Line<float> mid_line(mid_pt,dir_perp);
      
      // gather the points
      larocv::GEO2D_Contour_t nz_pt_v;
      for(const auto& poly : obj2d._polygon_v) {
	auto local_nz_pt_v = larocv::FindNonZero(larocv::MaskImage(cimg,poly.Contour(),-1,false));
	for(auto& local_nz_pt : local_nz_pt_v) 
	  nz_pt_v.emplace_back(std::move(local_nz_pt));
      }

      std::vector<const cv::Point_<int>* > side1_v;
      std::vector<const cv::Point_<int>* > side2_v;

      size_t nz_pt_half = (size_t)(((float)nz_pt_v.size())/2.0);
      side1_v.reserve(nz_pt_half);
      side2_v.reserve(nz_pt_half);
      
      // check for infinite slope
      bool isinf = false;
      if (mid_line.dir.x == 0)
	isinf = true;
      
      for(const auto& nz_pt : nz_pt_v) {
	
	if (!isinf) {
	  float y = mid_line.y(nz_pt.x);	  
	
	  if ((float)nz_pt.y >= y)
	    side1_v.push_back(&nz_pt);
	  else
	    side2_v.push_back(&nz_pt);
	  
	} else {
	  
	  if (nz_pt.x >= mid_pt.x) 
	    side1_v.push_back(&nz_pt);
	  else
	    side2_v.push_back(&nz_pt);
	}
      } // end sorting points  

      assert((side1_v.size() + side2_v.size()) == nz_pt_v.size());

      // determine which side the start is on
      std::vector<const cv::Point_<int>* >* near_side_v = nullptr;
      std::vector<const cv::Point_<int>* >* far_side_v = nullptr;
      
      if (!isinf) {
	float y = mid_line.y(obj2d.Start().x);
	if (obj2d.Start().y >= y) {
	  near_side_v = &side1_v;
	  far_side_v  = &side2_v;
	} 
	else {
	  near_side_v = &side2_v;
	  far_side_v  = &side1_v;
	}

      } else {

	if (obj2d.Start().x >= mid_pt.x) {
	  near_side_v = &side1_v;
	  far_side_v  = &side2_v;
	}
	else {
	  near_side_v = &side2_v;
	  far_side_v  = &side1_v;
	}

      }

      assert(near_side_v);
      assert(far_side_v);

      // paint the first half
      white_img.setTo(cv::Scalar(0));
      for(const auto near_side : *near_side_v)
	white_img.at<uchar>(near_side->y,near_side->x) = (uchar)255;
      
      // count the fraction of pixels contained in the line
      float near_cnt = (float)near_side_v->size();
      float near_line_cnt = larocv::CountNonZero(larocv::MaskImage(white_img,obj2d.Line(),-1,false));
      
      float first_half_linefrac = 0;
      if (near_cnt > 0)
	first_half_linefrac = near_line_cnt / near_cnt;
      else
	first_half_linefrac = -1;

      // paint the second hald half
      white_img.setTo(cv::Scalar(0));
      for(const auto far_side : *far_side_v)
	white_img.at<uchar>(far_side->y,far_side->x) = (uchar)255;
      
      // count the fraction of pixels contained in the line
      float far_cnt = (float)far_side_v->size();
      float far_line_cnt = larocv::CountNonZero(larocv::MaskImage(white_img,obj2d.Line(),-1,false));
      
      float second_half_linefrac = 0;
      if (far_cnt > 0)
	second_half_linefrac = far_line_cnt / far_cnt;
      else
	second_half_linefrac = -1;

      obj2d._line_first_half_linefrac = first_half_linefrac;
      obj2d._line_second_half_linefrac = second_half_linefrac;

    }
    
    
    return;
  }

  void FillTrack(const Object2DCollection& obj_col,
		 larlite::track& out_track) {
    
    const auto& start3d = obj_col.Start();

    TVector3 dir3d(obj_col.dX(),
		   obj_col.dY(),
		   obj_col.dZ());

    for(size_t plane=0; plane<3; ++plane) {
      
      std::vector<double> dqdx_v;
      
      if (plane==0) { // skip U plane
	out_track.add_dqdx(dqdx_v); 
	continue;
      }

      if (!obj_col.HasObject(plane)) {
	out_track.add_dqdx(dqdx_v);
	continue;
      }
      
      const auto& obj2d = obj_col.PlaneObject(plane);

      dqdx_v.reserve(obj2d.dQdxProfile().size());

      for(size_t xid=0; xid < obj2d.dQdxProfile().size(); ++xid) {

	auto dx = obj2d.dxProfile()[xid];

	auto step3d = dir3d*dx;
	
	auto xyz = start3d + step3d;
	
	// retrieve plane 1 from "DirectionAtPoint(p)"
	if (plane==1) 
	  out_track.add_direction(xyz);
	
	// retrieve plane 2 from "TrajectoryAtPoint(p)"
	if (plane==2) 
	  out_track.add_vertex(xyz);

	dqdx_v.push_back(obj2d.dQdxProfile()[xid]);
      }
      
      out_track.add_dqdx(dqdx_v);
    }

    return;
  }
  

}

#endif
