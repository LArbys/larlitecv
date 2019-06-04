#ifndef __CONTOURSCAN_CXX__
#define __CONTOURSCAN_CXX__

#include "ContourScan.h"
#include "LLCVBase/llcv_err.h"

namespace llcv {
  
  void ContourScan::Reset() {
    for(auto& ctor : _ctor_v)
      ctor.clear();

    _scan_v.clear();
    _end_pt_v.clear();
    _end_pt_plane_v.clear();
    return;
  }
  
  void ContourScan::RegisterContour(const cv::Mat& img, const larocv::GEO2D_Contour_t& ctor, const size_t plane, const float rad) {
    auto mask = larocv::MaskImage(img,ctor,-1,false);

    if (rad > 0)
      mask = larocv::MaskImage(img,geo2d::Circle<float>(200,200,rad),-1,false);

    _ctor_v[plane] = larocv::FindNonZero(mask);
    return;
  }
  
  void ContourScan::AddPixels(const larocv::GEO2D_Contour_t& pt_v, const size_t plane) {
    for(const auto& pt : pt_v)
      _ctor_v[plane].push_back(pt);
  }

  void ContourScan::Scan(const std::array<cv::Mat,3>& in_img_v,
			 const std::array<cv::Mat,3>& dead_img_v,
			 std::array<cv::Mat,3>& out_img_v) {

    larocv::data::Vertex3D res;

    for(const auto& plane_comb : _plane_comb_v) {
      auto plane0 = plane_comb[0];
      auto plane1 = plane_comb[1];
      auto plane2 = plane_comb[2];

      for(size_t pid0=0; pid0<_ctor_v[plane0].size(); ++pid0) {
	for(size_t pid1=0; pid1<_ctor_v[plane1].size(); ++pid1) {

	  if (!_geo.YZPoint(_ctor_v[plane0][pid0],plane0,
			    _ctor_v[plane1][pid1],plane1,
			    res)) continue;
	  

	  int px = (int)res.vtx2d_v[plane2].pt.x;
	  int py = (int)res.vtx2d_v[plane2].pt.y;

	  int img_pt   = (int)in_img_v[plane2].at<uchar>(py,px);
	  
	  int not_dead = (int)dead_img_v[plane2].at<uchar>(py,px);
					
	  if (img_pt<10 and not_dead) continue;
	  
	  out_img_v[plane0].at<uchar>(_ctor_v[plane0][pid0].y,
				      _ctor_v[plane0][pid0].x)  = (uchar) 255;
	  
	  out_img_v[plane1].at<uchar>(_ctor_v[plane1][pid1].y,
				      _ctor_v[plane1][pid1].x)  = (uchar) 255;

	  	  
	} // end plane0
      } // end plane1
    } // end plane combo

    
    return;
  }

  void ContourScan::Scan() {
    
    for(const auto&  plane_comb : _plane_comb_2d_v) {
      
      auto plane1 = plane_comb.first;
      auto plane2 = plane_comb.second;
      
      const auto& nz_pt1_v = _ctor_v[plane1];
      const auto& nz_pt2_v = _ctor_v[plane2];

      // std::cout << "@(plane1,plane2)=" 
      // 		<< "(" << plane1 << "," << plane2 << ")" 
      // 		<< "sz=(" << nz_pt1_v.size() << "," << nz_pt2_v.size() << ")" << std::endl;

      if (nz_pt1_v.empty()) continue;
      if (nz_pt2_v.empty()) continue;

      _scan_v.reserve(_scan_v.size() + nz_pt1_v.size()*nz_pt2_v.size());

      static larocv::data::Vertex3D res;
      for(const auto& nz_pt1 : nz_pt1_v) {
	for(const auto& nz_pt2 : nz_pt2_v) {
	  if (!_geo.YZPoint(nz_pt1,plane1,
			    nz_pt2,plane2,
			    res)) continue;
	  _scan_v.emplace_back(res);
	} // end pt2
      } // end pt1

      //std::cout << "... _scan_v sz=" << _scan_v.size() << std::endl;
    } // end plane_comb
    
    return;
  }
  

  std::vector<std::array<float,3> > ContourScan::Voxelize(const float dx, const float dy, const float dz) const {

    std::vector<std::array<float, 3> > ret_v;
    
    // get the minimum x,y,z
    double min_x = larocv::kINVALID_DOUBLE;
    double min_y = larocv::kINVALID_DOUBLE;
    double min_z = larocv::kINVALID_DOUBLE;

    double max_x = -1*larocv::kINVALID_DOUBLE;
    double max_y = -1*larocv::kINVALID_DOUBLE;
    double max_z = -1*larocv::kINVALID_DOUBLE;
    
    if (_scan_v.empty()) return ret_v;

    for(const auto& scan : _scan_v) {
      min_x = std::min(min_x,scan.x);
      min_y = std::min(min_y,scan.y);
      min_z = std::min(min_z,scan.z);
      max_x = std::max(max_x,scan.x);
      max_y = std::max(max_y,scan.y);
      max_z = std::max(max_z,scan.z);
    }
    
    //std::cout << "min=(" << min_x << "," << min_y << "," << min_z << ")" << std::endl;
    //std::cout << "max=(" << max_x << "," << max_y << "," << max_z << ")" << std::endl;

    int nbins_x = (int)((max_x - min_x) / dx + 0.5) + 1;
    int nbins_y = (int)((max_y - min_y) / dy + 0.5) + 1;
    int nbins_z = (int)((max_z - min_z) / dz + 0.5) + 1;
    
    std::vector<std::vector<std::vector<int> > > voxel_vvv;
    voxel_vvv.resize(nbins_x);
    for(auto& voxel_vv : voxel_vvv) {
      voxel_vv.resize(nbins_y);
      for(auto& voxel_v : voxel_vv) {
	voxel_v.resize(nbins_z,0);
      }
    }

    for(const auto& scan : _scan_v) {
      int binx = (int)(((scan.x - min_x) / dx) + 0.5);
      int biny = (int)(((scan.y - min_y) / dy) + 0.5);
      int binz = (int)(((scan.z - min_z) / dz) + 0.5);

      if (binx >= nbins_x or biny >= nbins_y or binz >= nbins_z) {
	std::stringstream ss;
	ss << "@(binx,biny,binz)="
	   << "(" << binx << "," << biny << "," << binz << ") "
	   << ">= " 
	   << "(" << nbins_x << "," << nbins_y << "," << nbins_z << ")";
	throw llcv_err(ss.str());
      }
      voxel_vvv.at(binx).at(biny).at(binz) = 1;
    }

    
    for(int bx=0; bx<nbins_x; ++bx) {
      for(int by=0; by<nbins_y; ++by) {
	for(int bz=0; bz<nbins_z; ++bz) {

	  if (!voxel_vvv.at(bx).at(by).at(bz)) 
	    continue;
	  
	  float pt_x, pt_y, pt_z;
	  
	  pt_x = (min_x + bx*dx) + (dx / 0.5);
	  pt_y = (min_y + by*dy) + (dy / 0.5);
	  pt_z = (min_z + bz*dz) + (dz / 0.5);
	  
	  std::array<float,3> pt_v;
	  pt_v[0] = pt_x;
	  pt_v[1] = pt_y;
	  pt_v[2] = pt_z;

	  ret_v.emplace_back(std::move(pt_v));
	}
      }
    }
    

    return ret_v;
  }

  void ContourScan::RegisterEndPoint(const geo2d::Vector<float>& pt,
				     const size_t plane) {
    
    _end_pt_v.push_back(pt);
    _end_pt_plane_v.push_back(plane);
    
    return;
  }

  std::array<float,3> ContourScan::EndPoint() const {

    std::array<float,3> ret_v;
    ret_v[0] = -1*larocv::kINVALID_FLOAT;
    ret_v[1] = -1*larocv::kINVALID_FLOAT;
    ret_v[2] = -1*larocv::kINVALID_FLOAT;
    
    if (_end_pt_v.size() == 2) {
      
      larocv::data::Vertex3D res;
      auto compat = _geo.YZPoint(_end_pt_v.front(),_end_pt_plane_v.front(),
				 _end_pt_v.back() ,_end_pt_plane_v.back(),
				 res);
    
      if (compat) {
	ret_v[0] = res.x;
	ret_v[1] = res.y;
	ret_v[2] = res.z;
      }

    } // end size == 2
    else if (_end_pt_v.size() == 3) {
    
      std::vector<larocv::data::Vertex3D> res_v;
      res_v.reserve(3);
      
      for(size_t pid1=0; pid1 < _end_pt_v.size(); ++pid1) {
	for(size_t pid2=pid1+1; pid2 < _end_pt_v.size(); ++pid2) {
	  
	  larocv::data::Vertex3D res;
	  auto compat = _geo.YZPoint(_end_pt_v[pid1], _end_pt_plane_v[pid1],
				     _end_pt_v[pid2], _end_pt_plane_v[pid2],
				     res);
	  
	  if (compat)
	    res_v.emplace_back(std::move(res));
	  
	}
      }

      // average
      ret_v[0] = 0;
      ret_v[1] = 0;
      ret_v[2] = 0;
      for(const auto& res : res_v) {
	ret_v[0] += res.x;
	ret_v[1] += res.y;
	ret_v[2] += res.z;
      }

      if (!res_v.empty()) {
	ret_v[0] /= (float)res_v.size();
	ret_v[1] /= (float)res_v.size();
	ret_v[2] /= (float)res_v.size();
      }

    } // end size == 3;
    
    
    return ret_v;
  }

}

#endif
