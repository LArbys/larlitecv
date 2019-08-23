#ifndef __SELTRIANGLESTUDY_CXX__
#define __SELTRIANGLESTUDY_CXX__

#include "SelTriangleStudy.h"

#include "InterTool_Util/InterImageUtils.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"

#include <iomanip>
#include <sstream>
#include <cstdlib>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

namespace llcv {

  void SelTriangleStudy::Configure(const larcv::PSet &pset) {

    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    
    _cropx = 400;
    _cropy = 400;

    _twatch.Stop();

    _n_neighbors  = pset.get<size_t>("NNeighbors",3);
    _brem_dist    = pset.get<float>("BremDistance",10);
    _brem_size    = pset.get<int>("BremSize",4);
    _min_defect_sz= pset.get<float>("MinDefectSize",5);

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelTriangleStudy::Initialize() {
    LLCV_DEBUG() << "start" << std::endl;

    _outtree = new TTree("SelTriangleStudy","");
    AttachRSEV(_outtree);

    _outtree->Branch("triangle_height_U_v",&_triangle_height_U_v);
    _outtree->Branch("triangle_height_V_v",&_triangle_height_V_v);
    _outtree->Branch("triangle_height_Y_v",&_triangle_height_Y_v);

    _outtree->Branch("triangle_base_U_v",&_triangle_base_U_v);
    _outtree->Branch("triangle_base_V_v",&_triangle_base_V_v);
    _outtree->Branch("triangle_base_Y_v",&_triangle_base_Y_v);

    _outtree->Branch("triangle_area_U_v",&_triangle_area_U_v);
    _outtree->Branch("triangle_area_V_v",&_triangle_area_V_v);
    _outtree->Branch("triangle_area_Y_v",&_triangle_area_Y_v);

    _outtree->Branch("triangle_straight_U_v",&_triangle_straight_U_v);
    _outtree->Branch("triangle_straight_V_v",&_triangle_straight_V_v);
    _outtree->Branch("triangle_straight_Y_v",&_triangle_straight_Y_v);

    _outtree->Branch("triangle_empty_area_ratio_U_v",&_triangle_empty_area_ratio_U_v);
    _outtree->Branch("triangle_empty_area_ratio_V_v",&_triangle_empty_area_ratio_V_v);
    _outtree->Branch("triangle_empty_area_ratio_Y_v",&_triangle_empty_area_ratio_Y_v);

    _outtree->Branch("triangle_empty_area_U_v",&_triangle_empty_area_U_v);
    _outtree->Branch("triangle_empty_area_V_v",&_triangle_empty_area_V_v);
    _outtree->Branch("triangle_empty_area_Y_v",&_triangle_empty_area_Y_v);

    _outtree->Branch("triangle_brem_U_v",&_triangle_brem_U_v);
    _outtree->Branch("triangle_brem_V_v",&_triangle_brem_V_v);
    _outtree->Branch("triangle_brem_Y_v",&_triangle_brem_Y_v);

    _outtree->Branch("polygon_number_defects_U_v",&_polygon_number_defects_U_v);
    _outtree->Branch("polygon_number_defects_V_v",&_polygon_number_defects_V_v);
    _outtree->Branch("polygon_number_defects_Y_v",&_polygon_number_defects_Y_v);

    _outtree->Branch("polygon_number_defects_no_start_U_v",&_polygon_number_defects_no_start_U_v);
    _outtree->Branch("polygon_number_defects_no_start_V_v",&_polygon_number_defects_no_start_V_v);
    _outtree->Branch("polygon_number_defects_no_start_Y_v",&_polygon_number_defects_no_start_Y_v);

    _outtree->Branch("polygon_largest_defect_U_v",&_polygon_largest_defect_U_v);
    _outtree->Branch("polygon_largest_defect_V_v",&_polygon_largest_defect_V_v);
    _outtree->Branch("polygon_largest_defect_Y_v",&_polygon_largest_defect_Y_v);

    _outtree->Branch("polygon_smallest_defect_U_v",&_polygon_smallest_defect_U_v);
    _outtree->Branch("polygon_smallest_defect_V_v",&_polygon_smallest_defect_V_v);
    _outtree->Branch("polygon_smallest_defect_Y_v",&_polygon_smallest_defect_Y_v);
    
    _outtree->Branch("polygon_largest_defect_no_start_U_v",&_polygon_largest_defect_no_start_U_v);
    _outtree->Branch("polygon_largest_defect_no_start_V_v",&_polygon_largest_defect_no_start_V_v);
    _outtree->Branch("polygon_largest_defect_no_start_Y_v",&_polygon_largest_defect_no_start_Y_v);

    _outtree->Branch("polygon_smallest_defect_no_start_U_v",&_polygon_smallest_defect_no_start_U_v);
    _outtree->Branch("polygon_smallest_defect_no_start_V_v",&_polygon_smallest_defect_no_start_V_v);
    _outtree->Branch("polygon_smallest_defect_no_start_Y_v",&_polygon_smallest_defect_no_start_Y_v);

    _outtree->Branch("polygon_empty_area_ratio_U_v",&_polygon_empty_area_ratio_U_v);
    _outtree->Branch("polygon_empty_area_ratio_V_v",&_polygon_empty_area_ratio_V_v);
    _outtree->Branch("polygon_empty_area_ratio_Y_v",&_polygon_empty_area_ratio_Y_v);

    _outtree->Branch("polygon_empty_area_U_v",&_polygon_empty_area_U_v);
    _outtree->Branch("polygon_empty_area_V_v",&_polygon_empty_area_V_v);
    _outtree->Branch("polygon_empty_area_Y_v",&_polygon_empty_area_Y_v);

    _outtree->Branch("polygon_pocket_area_U_v",&_polygon_pocket_area_U_v);
    _outtree->Branch("polygon_pocket_area_V_v",&_polygon_pocket_area_V_v);
    _outtree->Branch("polygon_pocket_area_Y_v",&_polygon_pocket_area_Y_v);

    _outtree->Branch("polygon_pocket_area_no_start_U_v",&_polygon_pocket_area_no_start_U_v);
    _outtree->Branch("polygon_pocket_area_no_start_V_v",&_polygon_pocket_area_no_start_V_v);
    _outtree->Branch("polygon_pocket_area_no_start_Y_v",&_polygon_pocket_area_no_start_Y_v);
    

    LLCV_DEBUG() << "end" << std::endl;

  }

  double SelTriangleStudy::Select() {
    LLCV_DEBUG() << "start" << std::endl;

    LLCV_DEBUG() << "(RSEV)=("<< Run() << "," << SubRun() << "," << Event() << "," << VertexID() << ")" << std::endl;    
    LLCV_DEBUG() << "VTX=(" 
		 << Data().Vertex()->X() << "," 
		 << Data().Vertex()->Y() << "," 
		 << Data().Vertex()->Z() << ")" << std::endl;

    int isnue = Tree().Scalar<int>("isnue");
    LLCV_DEBUG() << "Is this nue?: " <<  isnue << std::endl;
    if (isnue == 0)  {
      LLCV_DEBUG() << "skip" << std::endl;
      return 0.0;
    }
    
    //
    // get the image centered around the vertex
    //
    
    auto mat_v  = Image().Image<cv::Mat>(kImageADC,_cropx,_cropy);
    auto meta_v = Image().Image<larocv::ImageMeta>(kImageADC,_cropx,_cropy);
    auto img_v  = Image().Image<larcv::Image2D>(kImageADC,_cropx,_cropy);
    auto dead_v = Image().Image<cv::Mat>(kImageDead,_cropx,_cropy);

    for(const auto& meta : meta_v)
      _PixelScan3D.SetPlaneInfo(*meta);

    std::vector<cv::Mat> thresh_mat_v;
    std::vector<cv::Mat> mat3d_v;
    mat3d_v.reserve(3);

    thresh_mat_v.reserve(3);

    for(size_t plane=0; plane<3; ++plane) {
      const auto mat = mat_v[plane];
      auto img = larocv::Threshold(*mat,10,255);
      auto img3d = As8UC3(img);
      thresh_mat_v.emplace_back(std::move(img));      
      mat3d_v.emplace_back(std::move(img3d));
    }
    
    larocv::data::Vertex3D vtx3d;
    vtx3d.x = Data().Vertex()->X();
    vtx3d.y = Data().Vertex()->Y();
    vtx3d.z = Data().Vertex()->Z();
    
    //
    // Project the vertex onto the plane
    //
    larocv::GEO2D_Contour_t vertex_ctor;
    vertex_ctor.resize(3);

    for(size_t plane=0; plane<3; ++plane) {
      int px_x = kINVALID_INT;
      int px_y = kINVALID_INT;

      ProjectMat(img_v[plane]->meta(),
    		 Data().Vertex()->X(),
		 Data().Vertex()->Y(),
		 Data().Vertex()->Z(),
    		 px_x, px_y);

      vertex_ctor[plane] = cv::Point_<int>(px_y,px_x);
    }


    std::array<cv::Mat,3> timg_v;
    for(size_t plane=0; plane<3; ++plane) {
      const auto& thresh_mat = thresh_mat_v[plane];
      timg_v[plane] = thresh_mat;
      //timg_v[plane] = larocv::BlankImage(thresh_mat,255);
    }

    const auto& shr_v = Data().Showers();

    ResizeOutput(shr_v.size());

    for(size_t shrid=0; shrid < shr_v.size(); ++shrid) {

      const auto& shr = (*shr_v[shrid]);

      auto theta_phi = ShowerAngle(shr);
      
      float pi8 = 3.14159 / 8;
      float pi4 = 3.14159 / 4;
      
      float rad_min  = 0.5;
      float rad_max  = 50;
      float rad_step = 0.5;
      
      float theta_min = theta_phi.first - pi8;
      float theta_max = theta_phi.first + pi8;
      
      float phi_min = theta_phi.second - pi4;
      float phi_max = theta_phi.second + pi4;
      
      LLCV_DEBUG() << "Scanning Spheres" << std::endl;
      LLCV_DEBUG() << "rad_min=" << rad_min << " rad_max=" << rad_max << " rad_step=" << rad_step << std::endl;
      LLCV_DEBUG() << "theta_min=" << theta_min << " theta_max=" << theta_max << std::endl;
      LLCV_DEBUG() << "phi_min=" << phi_min << " phi_max=" << phi_max << std::endl;

      //
      // Build 3D object
      //
      _twatch.Start();
      
      _PixelScan3D.Reconfigure(rad_min, rad_max, rad_step,
			       theta_min, theta_max,
			       phi_min, phi_max);
            
      std::vector<larocv::data::Vertex3D> reg_v;
      
      reg_v = _PixelScan3D.SphereScan3D(timg_v,dead_v,vtx3d);
      
      _twatch.Stop();
      
      LLCV_DEBUG() << reg_v.size() << " spheres scanned in " 
		   << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;
    
      if (reg_v.empty()) continue;

      std::array<cv::Mat,3> timg3d_v;
      for(size_t plane=0; plane<3; ++plane)
	timg3d_v[plane] = larocv::BlankImage(timg_v[plane],0);

      for(size_t pid=0; pid < reg_v.size(); ++pid) {
	const auto& pt = reg_v[pid];
	for(size_t plane=0; plane<3; ++plane) {
	  int px_x = pt.vtx2d_v[plane].pt.x;
	  int px_y = pt.vtx2d_v[plane].pt.y;

	  auto& timg3d = timg3d_v[plane];

	  if (px_x < 0) continue;
	  if (px_y < 0) continue;
	  if (px_x >= timg3d.rows) continue;
	  if (px_y >= timg3d.cols) continue;

	  timg3d.at<uchar>(px_y,px_x) = (uchar)255;
	} // end plane
      } // end sphere point


      // Fill
      for(size_t fnid=0; fnid<_n_neighbors; ++fnid)
	FillNeighbors(timg_v,timg3d_v,vertex_ctor);
      
      for(size_t plane=0; plane<3; ++plane) {
	auto& mat3d = mat3d_v[plane];
	auto pts_v = larocv::FindNonZero(timg3d_v[plane]);
	
	for(const auto& pt : pts_v) 
	  mat3d.at<cv::Vec3b>(pt.y,pt.x) = {255,255,0};

      }

      // Cluster
      for(size_t plane=0; plane<3; ++plane) {
	const auto& timg3d = timg3d_v[plane];
	auto ctor_v = larocv::FindContours(timg3d);

	LLCV_DEBUG() << "@plane=" << plane << " found " << ctor_v.size() << " ctors" << std::endl;

	const auto vertex_pt = vertex_ctor[plane];

	//
	// get the contour closest to the vertex
	//
	float dist = kINVALID_FLOAT;
	size_t close_id = kINVALID_SIZE;
	for(size_t cid=0; cid< ctor_v.size(); ++cid) {
	  const auto& ctor = ctor_v[cid];
	  auto dist_tmp = larocv::Pt2PtDistance(vertex_pt,ctor);
	  if (dist_tmp < dist) {
	    dist = dist_tmp;
	    close_id = cid;
	  }
	}
	
	if (close_id == larocv::kINVALID_SIZE) {
	  LLCV_CRITICAL() << "No contour? Skip..." << std::endl;
	  continue;
	}
      
	larocv::GEO2D_ContourArray_t other_ctor_v;
	for(size_t cid=0; cid<ctor_v.size(); ++cid) {
	  if (close_id == cid) continue;
	  other_ctor_v.emplace_back(ctor_v[cid]);
	}

	Triangle tri(ctor_v[close_id],vertex_pt);
	tri.Tighten(timg3d,3,0.25);
	
	int brem_exists = DetectBrem(tri,other_ctor_v);
	int* brem_ptr = nullptr;
	if (plane==0) brem_ptr = &(_triangle_brem_U_v[shrid]);
	if (plane==1) brem_ptr = &(_triangle_brem_V_v[shrid]);
	if (plane==2) brem_ptr = &(_triangle_brem_Y_v[shrid]);
	(*brem_ptr) = brem_exists;

	Polygon poly(ctor_v[close_id],vertex_pt);
      
	auto& img3d = mat3d_v[plane];	
      
	cv::line(img3d,tri.Base1(),tri.Base2(),cv::Scalar(238,130,238),1);
	cv::line(img3d,tri.Apex() ,tri.Base1(),cv::Scalar(0,255,0),1);
	cv::line(img3d,tri.Apex() ,tri.Base2(),cv::Scalar(0,255,0),1);
	
	FillTriangle(timg3d,tri,shrid,plane);
	FillPolygon(poly,shrid,plane);

      } // end plane

    } // end shower
    
    LLCV_DEBUG() << "done!" << std::endl;
    
    //
    // Draw the pixel clusters
    //
    std::vector<larocv::GEO2D_ContourArray_t> ctor_vv;
    ctor_vv.resize(3);

    const auto& par_vv = Data().Particles();

    LLCV_DEBUG() << "n particles=" << par_vv.size() << std::endl;
    for(size_t parid=0; parid < par_vv.size(); ++parid) {
      LLCV_DEBUG() << "@parid=" << parid << std::endl;
      const auto& par_v = par_vv[parid];
      for(size_t plane=0; plane<3; ++plane) {
	LLCV_DEBUG() << "@plane=" << plane << std::endl;
	const auto& par = par_v[plane];
	if (par.empty()) continue;
	const auto meta   = meta_v[plane];
	auto ctor = AsContour(par,*meta);
	auto& ctor_v = ctor_vv[plane];
	LLCV_DEBUG() << "saved ctor sz=" << ctor.size() << std::endl;
	ctor_v.emplace_back(std::move(ctor));
      }
    }
    
    for(size_t plane=0; plane<3; ++plane) {
      auto& img3d = mat3d_v[plane];
      const auto& ctor_v = ctor_vv[plane];
      for(size_t cid=0; cid<ctor_v.size(); ++cid) {
	LLCV_DEBUG() << "draw cid=" << cid << " sz=" << ctor_v[cid].size() << std::endl;
	cv::Scalar color(std::rand() % 255,std::rand() % 255,std::rand() % 255);
	cv::drawContours(img3d,ctor_v,cid,color);
      }
    }

    //
    // Write the image
    //
    for(size_t plane=0; plane<3; ++plane) {
      auto& img3d = mat3d_v[plane];
      std::stringstream ss; 
      ss << "cpng/plane_img_" << Run() << "_" << SubRun() << "_" << Event() << "_" << VertexID() << "_" << plane << ".png";
      cv::imwrite(ss.str(),img3d);
    }
    
    _outtree->Fill();
    
    LLCV_DEBUG() << "end" << std::endl;
    return 0.0;
  }
  
  std::pair<float,float> SelTriangleStudy::ShowerAngle(const larlite::shower& shower) {
    
    std::pair<float,float> res;
    
    const auto& dir = shower.Direction();

    res.first  = std::acos( dir.Z() / dir.Mag());
    res.second = std::atan2(dir.X(),dir.Y()); // wtf
    
    return res;
  }

  
  void SelTriangleStudy::FillNeighbors(const std::array<cv::Mat,3>& parent_v,
				       std::array<cv::Mat,3>& child_v,
				       const larocv::GEO2D_Contour_t& vertex_ctor) const {
    
    for(size_t plane=0; plane<3; ++plane) {
      const auto& parent = parent_v[plane];
      auto& child = child_v[plane];
      const auto& vertex_pt = vertex_ctor[plane];

      child = larocv::MaskImage(child,geo2d::Circle<float>(vertex_pt,2),-1,true);

      auto pts_v = larocv::FindNonZero(child);
      
      for(const auto& pt : pts_v) {
	auto spt_v = Neighbors(pt,child);
	for(const auto& spt : spt_v) {
	  uint ret = (uint)parent.at<uchar>(spt.x,spt.y);
	  if (ret==0) continue;
	  child.at<uchar>(spt.x,spt.y) = (uchar)255;
	} // end neighbors
      } // end 3D point
      
    } // end plane
    
  }
  
  larocv::GEO2D_Contour_t SelTriangleStudy::Neighbors(const geo2d::Vector<int>& pt, const cv::Mat& mat) const {
    larocv::GEO2D_Contour_t res;
    res.reserve(9);
    
    int row = pt.y;
    int col = pt.x;

    static std::array<int, 3> row_v;
    static std::array<int, 3> col_v;

    row_v[0] = row;
    row_v[1] = row+1;
    row_v[2] = row-1;
    col_v[0] = col;
    col_v[1] = col+1;
    col_v[2] = col-1;
    
    for(size_t rid=0; rid<3; ++rid) {
      for(size_t cid=0; cid<3; ++cid) {
	if (row_v[rid] < 0 or row_v[cid]>mat.rows) continue;
	if (col_v[cid] < 0 or col_v[cid]>mat.cols) continue;
	res.emplace_back(row_v[rid],col_v[cid]);
      }
    }
    
    return res;
  }

  int SelTriangleStudy::DetectBrem(Triangle triangle, const larocv::GEO2D_ContourArray_t& other_ctor_v) {
    int res = 0;

    //
    // count the number of other contours inside or touching the expanded triangle
    //
    
    float fraction = (triangle.Height() + _brem_dist) / triangle.Height();

    triangle.Extend(fraction);
    
    auto tri_ctor = triangle.AsContour();

    for(size_t oid=0; oid < other_ctor_v.size(); ++oid) {
      const auto& other_ctor = other_ctor_v[oid];

      auto common_area = larocv::AreaOverlap(other_ctor,tri_ctor);
      if (common_area == 0) continue;

      int pixel_area = larocv::ContourPixelArea(other_ctor);
      if (pixel_area > _brem_size)
	res += 1;
      
    }
    
    return res;
  }


  void SelTriangleStudy::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _outtree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }
  
  void SelTriangleStudy::FillTriangle(const cv::Mat& img, const Triangle& tri, const size_t idx, const size_t plane) {
    
    switch(plane) {
      
    case 0: {

      auto& triangle_height_U           = _triangle_height_U_v[idx];
      auto& triangle_base_U             = _triangle_base_U_v[idx];
      auto& triangle_area_U             = _triangle_area_U_v[idx];
      auto& triangle_straight_U         = _triangle_straight_U_v[idx];
      auto& triangle_empty_area_ratio_U = _triangle_empty_area_ratio_U_v[idx];
      auto& triangle_empty_area_U       = _triangle_empty_area_U_v[idx];

      triangle_height_U           = tri.Height();
      triangle_base_U             = tri.BaseLength();
      triangle_area_U             = tri.Area();
      triangle_straight_U         = tri.StraightLineTest(img);
      triangle_empty_area_ratio_U = tri.EmptyAreaRatio();
      triangle_empty_area_U       = tri.EmptyArea();

      break;
    }

    case 1: {

      auto& triangle_height_V           = _triangle_height_V_v[idx];
      auto& triangle_base_V             = _triangle_base_V_v[idx];
      auto& triangle_area_V             = _triangle_area_V_v[idx];
      auto& triangle_straight_V         = _triangle_straight_V_v[idx];
      auto& triangle_empty_area_ratio_V = _triangle_empty_area_ratio_V_v[idx];
      auto& triangle_empty_area_V       = _triangle_empty_area_V_v[idx];

      triangle_height_V           = tri.Height();
      triangle_base_V             = tri.BaseLength();
      triangle_area_V             = tri.Area();
      triangle_straight_V         = tri.StraightLineTest(img);
      triangle_empty_area_ratio_V = tri.EmptyAreaRatio();
      triangle_empty_area_V       = tri.EmptyArea();

      break;
    }

    case 2: {
      
      auto& triangle_height_Y           = _triangle_height_Y_v[idx];
      auto& triangle_base_Y             = _triangle_base_Y_v[idx];
      auto& triangle_area_Y             = _triangle_area_Y_v[idx];
      auto& triangle_straight_Y         = _triangle_straight_Y_v[idx];
      auto& triangle_empty_area_ratio_Y = _triangle_empty_area_ratio_Y_v[idx];
      auto& triangle_empty_area_Y       = _triangle_empty_area_Y_v[idx];

      triangle_height_Y           = tri.Height();
      triangle_base_Y             = tri.BaseLength();
      triangle_area_Y             = tri.Area();
      triangle_straight_Y         = tri.StraightLineTest(img);
      triangle_empty_area_ratio_Y = tri.EmptyAreaRatio();
      triangle_empty_area_Y       = tri.EmptyArea();

      break;
    }

    default: {
      throw llcv_err("exception");
      break;
    }

    }
    
    
    return;
  }

  void SelTriangleStudy::FillPolygon(const Polygon& poly, const size_t idx, const size_t plane) {

    switch(plane) {
      
    case 0: {
      auto& polygon_number_defects_U           = _polygon_number_defects_U_v[idx];
      auto& polygon_number_defects_no_start_U  = _polygon_number_defects_no_start_U_v[idx];
      auto& polygon_largest_defect_U           = _polygon_largest_defect_U_v[idx];
      auto& polygon_smallest_defect_U          = _polygon_smallest_defect_U_v[idx];
      auto& polygon_largest_defect_no_start_U  = _polygon_largest_defect_no_start_U_v[idx];
      auto& polygon_smallest_defect_no_start_U = _polygon_smallest_defect_no_start_U_v[idx];
      auto& polygon_empty_area_ratio_U         = _polygon_empty_area_ratio_U_v[idx];
      auto& polygon_empty_area_U               = _polygon_empty_area_U_v[idx];
      auto& polygon_pocket_area_U              = _polygon_pocket_area_U_v[idx];
      auto& polygon_pocket_area_no_start_U     = _polygon_pocket_area_no_start_U_v[idx];

      polygon_number_defects_U           = poly.NumberDefects(_min_defect_sz);
      polygon_number_defects_no_start_U  = poly.NumberDefectsNoStart(_min_defect_sz);
      polygon_largest_defect_U           = poly.LargestDefect();
      polygon_smallest_defect_U          = poly.SmallestDefect();
      polygon_largest_defect_no_start_U  = poly.LargestDefectNoStart();
      polygon_smallest_defect_no_start_U = poly.SmallestDefectNoStart();
      polygon_empty_area_ratio_U         = poly.EmptyAreaRatio();
      polygon_empty_area_U               = poly.EmptyArea();
      polygon_pocket_area_U              = poly.PocketArea();
      polygon_pocket_area_no_start_U     = poly.PocketAreaNoStart();

      break;
    }

    case 1: {

      auto& polygon_number_defects_V           = _polygon_number_defects_V_v[idx];
      auto& polygon_number_defects_no_start_V  = _polygon_number_defects_no_start_V_v[idx];
      auto& polygon_largest_defect_V           = _polygon_largest_defect_V_v[idx];
      auto& polygon_smallest_defect_V          = _polygon_smallest_defect_V_v[idx];
      auto& polygon_largest_defect_no_start_V  = _polygon_largest_defect_no_start_V_v[idx];
      auto& polygon_smallest_defect_no_start_V = _polygon_smallest_defect_no_start_V_v[idx];
      auto& polygon_empty_area_ratio_V         = _polygon_empty_area_ratio_V_v[idx];
      auto& polygon_empty_area_V               = _polygon_empty_area_V_v[idx];
      auto& polygon_pocket_area_V              = _polygon_pocket_area_V_v[idx];
      auto& polygon_pocket_area_no_start_V     = _polygon_pocket_area_no_start_V_v[idx];

      polygon_number_defects_V           = poly.NumberDefects(_min_defect_sz);
      polygon_number_defects_no_start_V  = poly.NumberDefectsNoStart(_min_defect_sz);
      polygon_largest_defect_V           = poly.LargestDefect();
      polygon_smallest_defect_V          = poly.SmallestDefect();
      polygon_largest_defect_no_start_V  = poly.LargestDefectNoStart();
      polygon_smallest_defect_no_start_V = poly.SmallestDefectNoStart();
      polygon_empty_area_ratio_V         = poly.EmptyAreaRatio();
      polygon_empty_area_V               = poly.EmptyArea();
      polygon_pocket_area_V              = poly.PocketArea();
      polygon_pocket_area_no_start_V     = poly.PocketAreaNoStart();

      break;
    }

    case 2: {

      auto& polygon_number_defects_Y           = _polygon_number_defects_Y_v[idx];
      auto& polygon_number_defects_no_start_Y  = _polygon_number_defects_no_start_Y_v[idx];
      auto& polygon_largest_defect_Y           = _polygon_largest_defect_Y_v[idx];
      auto& polygon_smallest_defect_Y          = _polygon_smallest_defect_Y_v[idx];
      auto& polygon_largest_defect_no_start_Y  = _polygon_largest_defect_no_start_Y_v[idx];
      auto& polygon_smallest_defect_no_start_Y = _polygon_smallest_defect_no_start_Y_v[idx];
      auto& polygon_empty_area_ratio_Y         = _polygon_empty_area_ratio_Y_v[idx];
      auto& polygon_empty_area_Y               = _polygon_empty_area_Y_v[idx];
      auto& polygon_pocket_area_Y              = _polygon_pocket_area_Y_v[idx];
      auto& polygon_pocket_area_no_start_Y     = _polygon_pocket_area_no_start_Y_v[idx];

      polygon_number_defects_Y           = poly.NumberDefects(_min_defect_sz);
      polygon_number_defects_no_start_Y  = poly.NumberDefectsNoStart(_min_defect_sz);
      polygon_largest_defect_Y           = poly.LargestDefect();
      polygon_smallest_defect_Y          = poly.SmallestDefect();
      polygon_largest_defect_no_start_Y  = poly.LargestDefectNoStart();
      polygon_smallest_defect_no_start_Y = poly.SmallestDefectNoStart();
      polygon_empty_area_ratio_Y         = poly.EmptyAreaRatio();
      polygon_empty_area_Y               = poly.EmptyArea();
      polygon_pocket_area_Y              = poly.PocketArea();
      polygon_pocket_area_no_start_Y     = poly.PocketAreaNoStart();      

      break;
    }

    default: {
      throw llcv_err("exception");
      break;
    }

    }


    return;
  }


  void SelTriangleStudy::ResizeOutput(size_t sz) {

    
    //
    // clear
    //

    _triangle_height_U_v.clear();
    _triangle_height_V_v.clear();
    _triangle_height_Y_v.clear();

    _triangle_base_U_v.clear();
    _triangle_base_V_v.clear();
    _triangle_base_Y_v.clear();

    _triangle_area_U_v.clear();
    _triangle_area_V_v.clear();
    _triangle_area_Y_v.clear();

    _triangle_straight_U_v.clear();
    _triangle_straight_V_v.clear();
    _triangle_straight_Y_v.clear();
      
    _triangle_empty_area_ratio_U_v.clear();
    _triangle_empty_area_ratio_V_v.clear();
    _triangle_empty_area_ratio_Y_v.clear();

    _triangle_empty_area_U_v.clear();
    _triangle_empty_area_V_v.clear();
    _triangle_empty_area_Y_v.clear();

    _triangle_brem_U_v.clear();
    _triangle_brem_V_v.clear();
    _triangle_brem_Y_v.clear();

    _polygon_number_defects_U_v.clear();
    _polygon_number_defects_V_v.clear();
    _polygon_number_defects_Y_v.clear();

    _polygon_number_defects_no_start_U_v.clear();
    _polygon_number_defects_no_start_V_v.clear();
    _polygon_number_defects_no_start_Y_v.clear();

    _polygon_largest_defect_U_v.clear();
    _polygon_largest_defect_V_v.clear();
    _polygon_largest_defect_Y_v.clear();
    
    _polygon_smallest_defect_U_v.clear();
    _polygon_smallest_defect_V_v.clear();
    _polygon_smallest_defect_Y_v.clear();
    
    _polygon_largest_defect_no_start_U_v.clear();
    _polygon_largest_defect_no_start_V_v.clear();
    _polygon_largest_defect_no_start_Y_v.clear();
      
    _polygon_smallest_defect_no_start_U_v.clear();
    _polygon_smallest_defect_no_start_V_v.clear();
    _polygon_smallest_defect_no_start_Y_v.clear();

    _polygon_empty_area_ratio_U_v.clear();
    _polygon_empty_area_ratio_V_v.clear();
    _polygon_empty_area_ratio_Y_v.clear();
    
    _polygon_empty_area_U_v.clear();
    _polygon_empty_area_V_v.clear();
    _polygon_empty_area_Y_v.clear();

    _polygon_pocket_area_U_v.clear();
    _polygon_pocket_area_V_v.clear();
    _polygon_pocket_area_Y_v.clear();

    _polygon_pocket_area_no_start_U_v.clear();
    _polygon_pocket_area_no_start_V_v.clear();
    _polygon_pocket_area_no_start_Y_v.clear();
    
    
    //
    // resize
    //

    _triangle_height_U_v.resize(sz,-1);
    _triangle_height_V_v.resize(sz,-1);
    _triangle_height_Y_v.resize(sz,-1);

    _triangle_base_U_v.resize(sz,-1);
    _triangle_base_V_v.resize(sz,-1);
    _triangle_base_Y_v.resize(sz,-1);

    _triangle_area_U_v.resize(sz,-1);
    _triangle_area_V_v.resize(sz,-1);
    _triangle_area_Y_v.resize(sz,-1);

    _triangle_straight_U_v.resize(sz,-1);
    _triangle_straight_V_v.resize(sz,-1);
    _triangle_straight_Y_v.resize(sz,-1);
      
    _triangle_empty_area_ratio_U_v.resize(sz,-1);
    _triangle_empty_area_ratio_V_v.resize(sz,-1);
    _triangle_empty_area_ratio_Y_v.resize(sz,-1);

    _triangle_empty_area_U_v.resize(sz,-1);
    _triangle_empty_area_V_v.resize(sz,-1);
    _triangle_empty_area_Y_v.resize(sz,-1);

    _triangle_brem_U_v.resize(sz,-1);
    _triangle_brem_V_v.resize(sz,-1);
    _triangle_brem_Y_v.resize(sz,-1);

    _polygon_number_defects_U_v.resize(sz,-1);
    _polygon_number_defects_V_v.resize(sz,-1);
    _polygon_number_defects_Y_v.resize(sz,-1);

    _polygon_number_defects_no_start_U_v.resize(sz,-1);
    _polygon_number_defects_no_start_V_v.resize(sz,-1);
    _polygon_number_defects_no_start_Y_v.resize(sz,-1);

    _polygon_largest_defect_U_v.resize(sz,-1);
    _polygon_largest_defect_V_v.resize(sz,-1);
    _polygon_largest_defect_Y_v.resize(sz,-1);
    
    _polygon_smallest_defect_U_v.resize(sz,-1);
    _polygon_smallest_defect_V_v.resize(sz,-1);
    _polygon_smallest_defect_Y_v.resize(sz,-1);
    
    _polygon_largest_defect_no_start_U_v.resize(sz,-1);
    _polygon_largest_defect_no_start_V_v.resize(sz,-1);
    _polygon_largest_defect_no_start_Y_v.resize(sz,-1);
      
    _polygon_smallest_defect_no_start_U_v.resize(sz,-1);
    _polygon_smallest_defect_no_start_V_v.resize(sz,-1);
    _polygon_smallest_defect_no_start_Y_v.resize(sz,-1);

    _polygon_empty_area_ratio_U_v.resize(sz,-1);
    _polygon_empty_area_ratio_V_v.resize(sz,-1);
    _polygon_empty_area_ratio_Y_v.resize(sz,-1);
    
    _polygon_empty_area_U_v.resize(sz,-1);
    _polygon_empty_area_V_v.resize(sz,-1);
    _polygon_empty_area_Y_v.resize(sz,-1);

    _polygon_pocket_area_U_v.resize(sz,-1);
    _polygon_pocket_area_V_v.resize(sz,-1);
    _polygon_pocket_area_Y_v.resize(sz,-1);

    _polygon_pocket_area_no_start_U_v.resize(sz,-1);
    _polygon_pocket_area_no_start_V_v.resize(sz,-1);
    _polygon_pocket_area_no_start_Y_v.resize(sz,-1);

  }


}


#endif

