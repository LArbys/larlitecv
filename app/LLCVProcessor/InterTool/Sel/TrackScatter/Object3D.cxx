#ifndef __OBJECT3D_CXX__
#define __OBJECT3D_CXX__

#include "Object3D.h"

#include <array>

#include "LArOpenCV/ImageCluster/AlgoFunction/VectorAnalysis.h"

#include "TVector3.h"

#ifndef __CLING__
#ifndef __CINT__
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#endif
#endif

#include <cassert>

namespace llcv {
  
  Object3D::Object3D(const larocv::data::Vertex3D& s,
		     const std::vector<const larocv::data::Vertex3D*> & v) {
    
    _start = s;

    _pts_v.clear();
    _pts_v.resize(v.size());

    for(size_t vid=0; vid<v.size(); ++vid) 
      _pts_v[vid] = *(v[vid]);

    FillPCA();
  }


  void Object3D::FillPCA() {
  
    cv::Mat vertex_mat(_pts_v.size(), 3, CV_32FC1);
    
    for(size_t vtxid=0; vtxid<_pts_v.size(); ++vtxid) {
      vertex_mat.at<float>(vtxid, 0) = _pts_v[vtxid].x;
      vertex_mat.at<float>(vtxid, 1) = _pts_v[vtxid].y;
      vertex_mat.at<float>(vtxid, 2) = _pts_v[vtxid].z;
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

    mean_dir_v[0] = mean_vv[0][0] - _start.x;
    mean_dir_v[1] = mean_vv[0][1] - _start.y;
    mean_dir_v[2] = mean_vv[0][2] - _start.z;

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

    FillOOBB(mean_vv[0],eigen_vv);

    return;
  }


  void Object3D::FillOOBB(const std::array<float,3> mean_v, const std::array<std::array<float,3>, 3> eigen_vv) {
    _center.x = mean_v[0];
    _center.y = mean_v[1];
    _center.z = mean_v[2];

    float min_dist0 = -1.0*larocv::kINVALID_FLOAT;
    float min_dist1 = -1.0*larocv::kINVALID_FLOAT;
    float min_dist2 = -1.0*larocv::kINVALID_FLOAT;

    std::array<float,3> far_pt0 = {{0,0,0}};
    std::array<float,3> far_pt1 = {{0,0,0}};
    std::array<float,3> far_pt2 = {{0,0,0}};

    auto start_v = ToVector(_start);

    auto vtx0 = larocv::Sum(start_v,eigen_vv[0]);
    auto vtx1 = larocv::Sum(start_v,eigen_vv[1]);
    auto vtx2 = larocv::Sum(start_v,eigen_vv[2]);

    std::array<float,3> far_pt_tmp0,far_pt_tmp1,far_pt_tmp2;
    float dist_tmp0,dist_tmp1,dist_tmp2;

    _deviation_v.clear();
    _deviation_v.resize(_pts_v.size(),-1);

    for(size_t pid=0; pid<_pts_v.size(); ++pid) {
      const auto& pt =  _pts_v[pid];
      auto& deviation = _deviation_v[pid];

      auto apt = ToVector(pt);

      far_pt_tmp0 = larocv::ClosestPoint(start_v,vtx0,apt);
      far_pt_tmp1 = larocv::ClosestPoint(start_v,vtx1,apt);
      far_pt_tmp2 = larocv::ClosestPoint(start_v,vtx2,apt);

      dist_tmp0   = larocv::Distance(start_v,far_pt_tmp0);
      dist_tmp1   = larocv::Distance(start_v,far_pt_tmp1);
      dist_tmp2   = larocv::Distance(start_v,far_pt_tmp2);
      
      deviation = larocv::Distance(apt,far_pt_tmp0);

      if (dist_tmp0 > min_dist0) {
	min_dist0 = dist_tmp0;
	far_pt0 = far_pt_tmp0;
      }

      if (dist_tmp1 > min_dist1) {
	min_dist1 = dist_tmp1;
	far_pt1 = far_pt_tmp1;
      }

      if (dist_tmp2 > min_dist2) {
	min_dist2 = dist_tmp2;
	far_pt2 = far_pt_tmp2;
      }

    }

    _end.x = far_pt0[0];
    _end.y = far_pt0[1];
    _end.z = far_pt0[2];

    _length = larocv::Distance(far_pt0, start_v);

    _width1 = 2*larocv::Distance(far_pt1, start_v);
    _width2 = 2*larocv::Distance(far_pt2, start_v);

    far_pt1 = larocv::Sum(far_pt1,larocv::Scale(_length,eigen_vv.front()));
    far_pt2 = larocv::Sum(far_pt2,larocv::Scale(_length,eigen_vv.front()));

    _edge1.x = far_pt1[0];
    _edge1.y = far_pt1[1];
    _edge1.z = far_pt1[2];

    _edge2.x = far_pt2[0];
    _edge2.y = far_pt2[1];
    _edge2.z = far_pt2[2];

    auto eigen_dir = TVector3(eigen_vv.front()[0],
			      eigen_vv.front()[1],
			      eigen_vv.front()[2]);

    _theta = eigen_dir.Theta();
    _phi   = eigen_dir.Phi();

    float cos;

    cos  = 0;
    cos  = _width1*_width1;
    cos /= (2*larocv::QuadDiff2(start_v,far_pt1));
    cos  = 1 - cos;

    _oangle1 = std::acos(cos);

    cos  = 0;
    cos  = _width2*_width2;
    cos /= (2*larocv::QuadDiff2(start_v,far_pt2));
    cos  = 1 - cos;

    _oangle2 = std::acos(cos);

    return;

  }
  
  std::array<float,3> Object3D::ToVector(const larocv::data::Vertex3D& vtx) const
  { return {{(float)vtx.x,(float)vtx.y,(float)vtx.z}}; }

  std::array<float,3> Object3D::ToVector(const TVector3& pt) const
  { return {{(float)pt.X(),(float)pt.Y(),(float)pt.Z()}}; }

  std::vector<float> Object3D::TrackDeviation(const larlite::track& trk) const {

    std::vector<float> res_v;
    res_v.resize(_pts_v.size(),-1);

    // for each segment compute the closest distance to segment
    for(size_t pid=0; pid<_pts_v.size(); ++pid) {
      auto apt = ToVector(_pts_v[pid]);
      float dist = larocv::kINVALID_FLOAT;
      for(size_t pid=0; pid< trk.NumberTrajectoryPoints()-1; ++pid) {
	auto pt1 = ToVector(trk.LocationAtPoint(pid));
	auto pt2 = ToVector(trk.LocationAtPoint(pid+1));
	auto line_pt = larocv::ClosestPoint(pt1,pt2,apt);
	auto dist_tmp= larocv::Distance(line_pt,apt);
	if (dist_tmp < dist)
	  dist = dist_tmp;
      }
      res_v[pid] = dist;
    }
    
    return res_v;
  }

}

#endif
