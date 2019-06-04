#ifndef __SHOWERTOOLS_CXX__
#define __SHOWERTOOLS_CXX__

#include "ShowerTools.h"
#include "LArUtil/GeometryHelper.h"
#include "InterTool_Util/InterImageUtils.h"
#include <cassert>
#include <numeric>
#include "Geo2D/Core/Geo2D.h"
#include "LLCVBase/llcv_err.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/VectorAnalysis.h"

namespace llcv {

  // adapted from https://goo.gl/Q7MBpE

  void ShowerTools::ReconstructAngle(const std::vector<larcv::Image2D*>& img_v,
				     const std::array<cv::Mat,3>& aimg_v,
				     Object2DCollection& obj_col) {

    auto geomH = larutil::GeometryHelper::GetME();
    
    // planes with largest number of hits used to get 3D direction
    std::vector<int> planeHits(3,0);
    std::vector<larutil::Point2D> planeDir(3);

    for(const auto& obj : obj_col) {

      const auto pl = obj._plane;
      
      assert(pl < 3);

      larutil::Point2D weightedDir;
      weightedDir.w = 0;
      weightedDir.t = 0;
      float Qtot = 0;
      
      int nhits = 0;
      
      for(const auto& poly : obj._polygon_v) {
	auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(aimg_v[pl],poly.Contour(),-1,false));
	
	for (const auto& nz_pt  : nz_pt_v){
	  float charge = MatToImage2DPixel(nz_pt,aimg_v[pl],*(img_v[pl]));

	  if (charge<=0) continue;
	  weightedDir.w += (nz_pt.y - obj.Start().y) * charge;
	  weightedDir.t += (nz_pt.x - obj.Start().x) * charge;

	  Qtot += charge;
	  nhits++;
	}
      }
     

      if (Qtot==0)  {
        LLCV_WARNING() << "Qtot=0!" << std::endl;
        Qtot=1; 
      }
      weightedDir.w /= Qtot;
      weightedDir.t /= Qtot;

      planeHits[pl] = nhits;
      planeDir[pl]  = weightedDir;
    }

    int pl_max = larlite::data::kINVALID_INT;
    int pl_mid = larlite::data::kINVALID_INT;
    int pl_min = larlite::data::kINVALID_INT;

    int n_max  = -1.0*larlite::data::kINVALID_INT;
    int n_min  =      larlite::data::kINVALID_INT;

    for (size_t pl=0; pl < planeHits.size(); pl++){
      if (planeHits[pl] > n_max){
	pl_max = pl;
	n_max  = planeHits[pl];
      }
      if (planeHits[pl] < n_min){
	pl_min = pl;
	n_min  = planeHits[pl];
      }
    }

    assert(pl_max != larlite::data::kINVALID_INT);
    assert(pl_min != larlite::data::kINVALID_INT);

    // find the medium plane
    for(int pp=0; pp<3; ++pp) {
      if (pp == pl_max) continue;
      if (pp == pl_min) continue;
      pl_mid = pp;
    }

    assert(pl_mid != larlite::data::kINVALID_INT);

    float slope_max, slope_mid;
    float angle_max, angle_mid;

    slope_mid = planeDir[pl_mid].t / planeDir[pl_mid].w;
    slope_max = planeDir[pl_max].t / planeDir[pl_max].w;
    
    angle_mid = std::atan(slope_mid);
    angle_mid = std::atan2( planeDir[pl_mid].t , planeDir[pl_mid].w );
    
    angle_max = std::atan(slope_max);
    angle_max = std::atan2( planeDir[pl_max].t , planeDir[pl_max].w );
    
    double theta, phi;
    geomH->Get3DAxisN(pl_max, pl_mid,
		      angle_max, angle_mid,
		      phi, theta);
    
    obj_col.SetTheta(theta);
    obj_col.SetPhi(phi);
    
  }
  
  // adapted from https://goo.gl/1pCDEa

  void ShowerTools::ReconstructLength(const std::vector<larcv::Image2D*>& img_v,
				      const std::array<cv::Mat,3>& aimg_v,
				      Object2DCollection& obj_col) {
    
    auto geomH = larutil::GeometryHelper::GetME();

    // output
    std::array<double,3> length_v;

    // geometry
    std::array<double,3> plane_f_v;
    std::array<TVector2,3> line_dir_v;

    //initialize
    for(auto& v : length_v) 
      v = -1.0*::larlite::data::kINVALID_DOUBLE;
    
    plane_f_v = length_v;
    
    // calcluate line
    float alpha = 5;

    TVector3 dir3D(obj_col.dX(),obj_col.dY(),obj_col.dZ());

    for(size_t plane=0; plane<3; ++plane) {
      plane_f_v[plane] = geomH->Project_3DLine_OnPlane(dir3D, plane).Mag();
      
      const auto& startPoint2D = obj_col.front().Start();
      TVector3 secondPoint3D = obj_col.Start() + alpha * dir3D;
      
      int px_x, px_y;
      ProjectMat(img_v[plane]->meta(),
		 secondPoint3D.X(),secondPoint3D.Y(),secondPoint3D.Z(),
		 px_x, px_y);
      
      geo2d::Vector<float> secondPoint2D(px_y,px_x);
      LLCV_DEBUG() << "@plane=" << plane << " start=" << startPoint2D << " end=" << secondPoint2D << std::endl;

      TVector2 dir(secondPoint2D.y - startPoint2D.y, secondPoint2D.x - startPoint2D.x);
      dir *= 1.0 / dir.Mod();
      line_dir_v[plane] = dir;
    }

    // calculate the length
    for(const auto& obj : obj_col) {

      const auto pl = obj._plane;

      const auto& dr_w = line_dir_v[pl].X();
      const auto& dr_t = line_dir_v[pl].Y();
      
      std::vector<std::pair<float,float> > dist_v;

      float qsum = 0;
      float d2D = 0;
      
      for(const auto& poly : obj._polygon_v) {
	auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(aimg_v[pl],poly.Contour(),-1,false));
	for (const auto& nz_pt  : nz_pt_v) {
	  float charge = MatToImage2DPixel(nz_pt,aimg_v[pl],*(img_v[pl]));
	  if (charge<=0) continue;
	  qsum += charge;
	}
      }

      for(const auto& poly : obj._polygon_v) {
	auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(aimg_v[pl],poly.Contour(),-1,false));
     
	for(size_t nid=0; nid < nz_pt_v.size(); ++nid) {
	  auto nz_pt = nz_pt_v[nid];

	  float charge = MatToImage2DPixel(nz_pt,aimg_v[pl],*(img_v[pl]));

	  if (charge<=0) continue;
	  
	  float ptw = (nz_pt.y - obj.Start().y)*0.3;
	  float ptt = (nz_pt.x - obj.Start().x)*0.3;
	
	  // calculate distance along the line
	  d2D  = (ptw * dr_w) + (ptt * dr_t);
	  d2D  = std::abs(d2D);
	  dist_v.emplace_back(std::make_pair(d2D, charge / qsum));
	}
      } // end "hit" loop
      
      std::sort(std::begin(dist_v),std::end(dist_v),
		[](const std::pair<float,float>& lhs, const std::pair<float,float>& rhs)
		{ return lhs.first < rhs.first; });
      
      qsum = 0;
      d2D = 0;
      float _qfraction = 1.0;
      for(const auto& dist_pair : dist_v) {
	d2D   = dist_pair.first;
	qsum += dist_pair.second;
	if (qsum > _qfraction) break;
      }
      
      auto f = plane_f_v.at(pl);
      
      double length = d2D / f;

      LLCV_DEBUG() << "@pl=" << pl << " f=" << f << " length=" << length << std::endl;
      length_v[pl] = length;
    }
    
    float sum = 0.0;
    float nplanes = 0.0;

    for(size_t plane=0; plane<3; ++plane) {
      const auto length = length_v[plane];
      
      if (length == -1.0*::larlite::data::kINVALID_DOUBLE)
	continue;

      auto& obj2d = obj_col.PlaneObjectRW(plane);
      obj2d._length = length;

      nplanes += 1;
      sum += (float)length;
    }
    
    sum /= nplanes;

    obj_col.SetLength(sum);
  }


  // adapted from https://goo.gl/fb894P


  void ShowerTools::ReconstructdQdx(const std::vector<larcv::Image2D*>& img_v,
				    const std::array<cv::Mat,3>& aimg_v,
				    Object2DCollection& obj_col,
				    double dtrunk) {
    
    auto geomHelper = larutil::GeometryHelper::GetME();
    
    TVector3 dir3D(obj_col.dX(),obj_col.dY(),obj_col.dZ());

    std::array<double,3> plane_f_v;
    std::array<double,3> pitch_v;
    
    for(size_t plane=0; plane<3; ++plane) {
      plane_f_v[plane] = geomHelper->Project_3DLine_OnPlane(dir3D, plane).Mag();
      pitch_v[plane]   = geomHelper->GetPitch(dir3D,plane);
    }
    
    // loop through planes
    for(auto& obj : obj_col) {

      const auto pl = obj._plane;

      // grab the 2D start point of the cluster
      const auto& start2D = obj.Start();
      
      double f      = plane_f_v[pl];
      double pitch  = pitch_v[pl];

      std::vector<double> dqdx_v(3*dtrunk, 0.);
      
      // loop through hits and find those within some radial distance of the start point
      // loop over hits
      size_t nhits = 0;
      for(const auto& poly : obj._polygon_v) {
	auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(aimg_v[pl],poly.Contour(),-1,false));

	for(size_t nid=0; nid < nz_pt_v.size(); ++nid) {
	  const auto& nz_pt = nz_pt_v[nid];

	  double d2D = std::sqrt( std::pow((nz_pt.y - start2D.y)*0.3, 2) + 
				  std::pow((nz_pt.x - start2D.x)*0.3, 2) );
	  
	  double d3D = d2D / f;
	  
	  if (d3D >= dtrunk) continue;
	  
	  auto charge = MatToImage2DPixel(nz_pt,aimg_v[pl],*(img_v[pl]));

	  if (charge <= 0) continue;
	  
	  dqdx_v.at((size_t)(d3D * 3)) += charge;

	  nhits++;
	} // end "hit"
      } // loop over all polygons

      static std::vector<double> dqdx_nonzero_v;
      dqdx_nonzero_v.clear();
      dqdx_nonzero_v.reserve(nhits);

      for (const auto& dqdx : dqdx_v) {
	if (dqdx) dqdx_nonzero_v.push_back(dqdx); 
      }
      
      double dqdx = -1.0*larlite::data::kINVALID_DOUBLE;
      
      if (dqdx_nonzero_v.empty()) dqdx = 0.;
      else {

	size_t half_idx = (size_t)(((float)dqdx_nonzero_v.size()) / 2.0);

	std::nth_element(dqdx_nonzero_v.begin(),
			 dqdx_nonzero_v.begin() + half_idx,
			 dqdx_nonzero_v.end());

	dqdx  = dqdx_nonzero_v.at(half_idx);
	dqdx /= pitch;
      }
      
      obj._dqdx = dqdx;
      
    } // end plane
    
    return;
  }
  
  void ShowerTools::ReconstructdQdxProfile(const std::vector<larcv::Image2D*>& img_v,
					   const std::array<cv::Mat,3>& aimg_v,
					   Object2DCollection& obj_col) {
    
    auto geomHelper = larutil::GeometryHelper::GetME();
    
    TVector3 dir3D(obj_col.ddX(),obj_col.ddY(),obj_col.ddZ());
    
    std::array<geo2d::Line<float>,3> plane_r_v;
    std::array<bool,3> valid_v;
    std::array<float,3> plane_f_v;
    
    for(size_t plane=0; plane<3; ++plane) {

      int pt1_x;
      int pt1_y;

      int pt2_x;
      int pt2_y;
      
      // project the vertex
      ProjectMat(img_v[plane]->meta(),
		 obj_col.Start().X(),
		 obj_col.Start().Y(),
		 obj_col.Start().Z(),
    		 pt1_y, pt1_x);
      
      // project the next point
      ProjectMat(img_v[plane]->meta(),
		 obj_col.Start().X() + 10*dir3D.X(),
		 obj_col.Start().Y() + 10*dir3D.Y(),
		 obj_col.Start().Z() + 10*dir3D.Z(),
    		 pt2_y, pt2_x);

      geo2d::Vector<float> pt1(pt1_x,pt1_y);
      geo2d::Vector<float> pt2(pt2_x,pt2_y);
      
      auto ldir = pt2-pt1;
      auto ldir_len = geo2d::length(ldir);
      if (ldir_len!=0)
	ldir /= geo2d::length(ldir);
      
      valid_v[plane] = true;
      if (ldir.x==0 and ldir.y==0)
	valid_v[plane] = false;

      // std::cout << "ldir=(" << ldir.x << "," << ldir.y << ") valid=" << valid_v[plane] << std::endl;
      plane_r_v[plane] = geo2d::Line<float>(pt1,ldir);
      plane_f_v[plane] = geomHelper->Project_3DLine_OnPlane(dir3D, plane).Mag();
    }
    
    // loop through planes
    for(auto& obj : obj_col) {

      const auto pl = obj._plane;

      if (!valid_v[pl]) continue;

      // get the 2D start point of the cluster
      const auto& start2D = obj.Start();
      
      const auto& line = plane_r_v[pl];
      float f = plane_f_v[pl];
      
      std::vector<float> q_v;
      std::vector<float> dx_v;

      larocv::GEO2D_Contour_t nz_pt_v;

      // put in the near vertex region
      for(const auto& nz_pt_line : obj._vtx_pt_v)
	nz_pt_v.push_back(nz_pt_line);

      // put in the polygons
      for(const auto& poly : obj._polygon_v) {
	auto nz_v = larocv::FindNonZero(larocv::MaskImage(aimg_v[pl],poly.Contour(),-1,false));
	for(auto& nz : nz_v) nz_pt_v.emplace_back(std::move(nz));
      }

      // loop over ``hits''
      for(const auto& nz_pt : nz_pt_v) {
	  
	geo2d::Vector<float> pt_nz(nz_pt.x,nz_pt.y);
	geo2d::Vector<float> pt_li;
	  
	geo2d::ClosestPoint(line, pt_nz, pt_li, pt_nz);
	  
	float q = MatToImage2DPixel(nz_pt,aimg_v[pl],*(img_v[pl]));
	  
	q_v.emplace_back(q);
	  
	auto dx = geo2d::dist(pt_li*0.3,start2D*0.3);
	  
	dx_v.emplace_back(dx);
	  
      } // end "hit"


      if (dx_v.empty()) continue;

      // sort the dx vector
      std::vector<size_t> idx_v(dx_v.size());
      std::iota(idx_v.begin(), idx_v.end(), 0);

      std::sort(idx_v.begin(), idx_v.end(),
		[&dx_v](size_t i1, size_t i2) { return dx_v[i1] < dx_v[i2]; });
      
      float min_dx = 0;
      float max_dx = dx_v.at(idx_v.back());
      max_dx /= f;

      float ddx = 0.3;
      size_t xlo = 0;
      size_t xhi = (size_t)((max_dx / ddx)+0.5)+1;

      if (xhi == 0) continue;

      std::vector<float> dqdx_v(xhi,0.0);
      std::vector<float> ddx_v(xhi,0.0);
      
      for(size_t ix=1; ix < ddx_v.size(); ++ix)
	ddx_v[ix] = ddx_v[ix-1] + ddx;
      
      for(auto idx : idx_v) {
	auto q  = q_v.at(idx);
	auto dx = dx_v.at(idx);

	dx /= f;
	
	int bin = (size_t)((dx / ddx)+0.5);

	if (bin >= (int)ddx_v.size() or bin < 0) {
	  std::stringstream ss;
	  ss << "@bin=" << bin << " & ddx_v sz=" << ddx_v.size() << "!";
	  throw llcv_err(ss.str());
	}
	
	// correct for pitch

	dqdx_v.at(bin) += q  * (f / ddx);
      }
      
      obj._dqdx_v = std::move(dqdx_v);
      obj._dx_v   = std::move(ddx_v);
      
      obj._dqdx_step = ddx;
      obj._dqdx_pitch = f;

      obj._dir = line.dir;

    } // end plane
    
    return;
  }


  void ShowerTools::TruncatedQdxProfile(Object2DCollection& obj_col, const float ftsigma) {
    std::vector<float> tdqdx_v;
    
    for(auto& obj : obj_col) {
      tdqdx_v.clear();
      _TruncMean.CalcTruncMeanProfile(obj._dx_v, obj._dqdx_v, tdqdx_v, ftsigma);
      obj._tdqdx_v = std::move(tdqdx_v);
    }

    return;
  }
  

  std::array<float,3> ShowerTools::ComputePCA(std::vector<std::array<float,3> > pts_v, const Object2DCollection& obj_col) {
  
    std::array<float,3> ret_v;
    ret_v[0] = -1*larocv::kINVALID_FLOAT;
    ret_v[1] = -1*larocv::kINVALID_FLOAT;
    ret_v[2] = -1*larocv::kINVALID_FLOAT;

    if (pts_v.empty()) return ret_v;

    cv::Mat vertex_mat(pts_v.size(), 3, CV_32FC1);

    for(size_t pid=0; pid < pts_v.size(); ++pid) {
      vertex_mat.at<float>(pid, 0) = pts_v[pid][0];
      vertex_mat.at<float>(pid, 1) = pts_v[pid][1];
      vertex_mat.at<float>(pid, 2) = pts_v[pid][2];
    }

    cv::PCA pca_ana(vertex_mat, cv::Mat(), CV_PCA_DATA_AS_ROW, 0);

    std::array<std::array<float,3>, 3> mean_vv;
    std::array<std::array<float,3>, 3> eigen_vv;
    
    for(size_t pid=0; pid<3; ++pid) {
      for(size_t plane=0; plane<3; ++plane) {
	mean_vv[pid][plane]  = pca_ana.mean.at<float>(pid,plane);
	eigen_vv[pid][plane] = pca_ana.eigenvectors.at<float>(pid,plane);
      }
      
      float eigen_len = larocv::Norm(eigen_vv[pid]);
      for(size_t plane=0; plane<3; ++plane)
	eigen_vv[pid][plane] /= eigen_len;
    }

    std::array<float,3> mean_dir_v;

    mean_dir_v[0] = mean_vv[0][0] - obj_col.Start().X();
    mean_dir_v[1] = mean_vv[0][1] - obj_col.Start().Y();
    mean_dir_v[2] = mean_vv[0][2] - obj_col.Start().Z();

    auto mean_dir_len = larocv::Norm(mean_dir_v);

    for(size_t plane=0; plane<3; ++plane)
      mean_dir_v[plane] /= mean_dir_len;

    // implement direction handling, negative dot product, flip the sign on eigen
    float dot_product = larocv::Dot(mean_dir_v,eigen_vv[0]);

    if (dot_product < 0) {
      eigen_vv[0][0] *= -1;
      eigen_vv[0][1] *= -1;
      eigen_vv[0][2] *= -1;
    }

    return mean_dir_v;
  }  

  
}
#endif
