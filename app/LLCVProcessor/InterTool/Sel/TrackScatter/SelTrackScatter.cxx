#ifndef __SELTRACKSCATTER_CXX__
#define __SELTRACKSCATTER_CXX__

#include "SelTrackScatter.h"

#include "InterTool_Util/InterImageUtils.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoClass/PixelChunk.h"

#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <array>
#include <algorithm>
#include <unordered_set>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

namespace llcv {

  void SelTrackScatter::Configure (const larcv::PSet &pset) {

    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    
    _cropx = 400;
    _cropy = 400;

    auto lrdb_pset = pset.get<larcv::PSet>("LongRangeDB");
    auto srdb_pset = pset.get<larcv::PSet>("ShortRangeDB");

    _LR_DefectBreaker.Configure(lrdb_pset.get<int>("MinDefectSize"),
				lrdb_pset.get<int>("NHullEdgePoints"),
				lrdb_pset.get<int>("NAllowedBreaks"));

    _SR_DefectBreaker.Configure(srdb_pset.get<int>("MinDefectSize"),
				srdb_pset.get<int>("NHullEdgePoints"),
				srdb_pset.get<int>("NAllowedBreaks"));
    
    _debug            = pset.get<bool>("Debug");
    _fill2d           = pset.get<bool>("Fill2D");
    _allow_dead_image = pset.get<bool>("AllowDeadImage");
    _skeletonize      = pset.get<bool>("Skeletonize");
    _sub_skeleton     = pset.get<bool>("SubSkeleton");

    _DBSCAN.Configure(3,5);
    
    _twatch.Stop();

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelTrackScatter::Initialize() {
    LLCV_DEBUG() << "start" << std::endl;

    _outtree = new TTree("scatter","");
    AttachRSEV(_outtree);

    if(_debug) {
      _outtree->Branch("track_x_vv", &_track_x_vv);
      _outtree->Branch("track_y_vv", &_track_y_vv);
      _outtree->Branch("track_z_vv", &_track_z_vv);
      
      _outtree->Branch("shower_x_vv", &_shower_x_vv);
      _outtree->Branch("shower_y_vv", &_shower_y_vv);
      _outtree->Branch("shower_z_vv", &_shower_z_vv);

      _outtree->Branch("shower_skel_x_vv", &_shower_skel_x_vv);
      _outtree->Branch("shower_skel_y_vv", &_shower_skel_y_vv);
      _outtree->Branch("shower_skel_z_vv", &_shower_skel_z_vv);

      _outtree->Branch("shower_p0_x_vv", &_shower_p0_x_vv);
      _outtree->Branch("shower_p0_y_vv", &_shower_p0_y_vv);
      _outtree->Branch("shower_p0_z_vv", &_shower_p0_z_vv);
      
      _outtree->Branch("shower_p1_x_vv", &_shower_p1_x_vv);
      _outtree->Branch("shower_p1_y_vv", &_shower_p1_y_vv);
      _outtree->Branch("shower_p1_z_vv", &_shower_p1_z_vv);

      _outtree->Branch("shower_p2_x_vv", &_shower_p2_x_vv);
      _outtree->Branch("shower_p2_y_vv", &_shower_p2_y_vv);
      _outtree->Branch("shower_p2_z_vv", &_shower_p2_z_vv);
      
      _outtree->Branch("shower_start_x_vv" , &_shower_start_x_vv);
      _outtree->Branch("shower_start_y_vv" , &_shower_start_y_vv);
      _outtree->Branch("shower_start_z_vv" , &_shower_start_z_vv);
      
      _outtree->Branch("shower_center_x_vv", &_shower_center_x_vv);
      _outtree->Branch("shower_center_y_vv", &_shower_center_y_vv);
      _outtree->Branch("shower_center_z_vv", &_shower_center_z_vv);
      
      _outtree->Branch("shower_end_x_vv"   , &_shower_end_x_vv);
      _outtree->Branch("shower_end_y_vv"   , &_shower_end_y_vv);
      _outtree->Branch("shower_end_z_vv"   , &_shower_end_z_vv);
      
      _outtree->Branch("shower_edge1_x_vv"   , &_shower_edge1_x_vv);
      _outtree->Branch("shower_edge1_y_vv"   , &_shower_edge1_y_vv);
      _outtree->Branch("shower_edge1_z_vv"   , &_shower_edge1_z_vv);
    
      _outtree->Branch("shower_edge2_x_vv"   , &_shower_edge2_x_vv);
      _outtree->Branch("shower_edge2_y_vv"   , &_shower_edge2_y_vv);
      _outtree->Branch("shower_edge2_z_vv"   , &_shower_edge2_z_vv);

      _outtree->Branch("shower_pca_dev_vv", &_shower_pca_dev_vv);
      _outtree->Branch("shower_trk_dev_vv", &_shower_trk_dev_vv);

      _outtree->Branch("shower_cid_vv", &_shower_cid_vv);
    }

    _outtree->Branch("shower3D_n_points_v", &_shower3D_n_points_v);

    _outtree->Branch("shower3D_length_v", &_shower3D_length_v);
    _outtree->Branch("shower3D_width_v",  &_shower3D_width_v);
    _outtree->Branch("shower3D_width1_v", &_shower3D_width1_v);
    _outtree->Branch("shower3D_width2_v", &_shower3D_width2_v);

    _outtree->Branch("shower3D_theta_v", &_shower3D_theta_v);
    _outtree->Branch("shower3D_phi_v",   &_shower3D_phi_v);

    _outtree->Branch("shower3D_opening_v",  &_shower3D_opening_v);
    _outtree->Branch("shower3D_opening1_v", &_shower3D_opening1_v);
    _outtree->Branch("shower3D_opening2_v", &_shower3D_opening2_v);
    
    _outtree->Branch("shower3D_pca_mean_dev_v",         &_shower3D_pca_mean_dev_v);
    _outtree->Branch("shower3D_start_pca_mean_dev_v",   &_shower3D_start_pca_mean_dev_v);
    _outtree->Branch("shower3D_middle_pca_mean_dev_v",  &_shower3D_middle_pca_mean_dev_v);
    _outtree->Branch("shower3D_end_pca_mean_dev_v",     &_shower3D_end_pca_mean_dev_v);

    _outtree->Branch("shower3D_track_mean_dev_v",        &_shower3D_track_mean_dev_v);      
    _outtree->Branch("shower3D_start_track_mean_dev_v",  &_shower3D_start_track_mean_dev_v);
    _outtree->Branch("shower3D_middle_track_mean_dev_v", &_shower3D_middle_track_mean_dev_v);
    _outtree->Branch("shower3D_end_track_mean_dev_v",    &_shower3D_end_track_mean_dev_v);

    _outtree->Branch("shower3D_n_clusters_v", &_shower3D_n_clusters_v);
    
    _outtree->Branch("shower3D_cluster_n_points_vv", &_shower3D_cluster_n_points_vv);

    _outtree->Branch("shower3D_cluster_length_vv" , &_shower3D_cluster_length_vv);
    _outtree->Branch("shower3D_cluster_width_vv"  , &_shower3D_cluster_width_vv);
    _outtree->Branch("shower3D_cluster_width1_vv" , &_shower3D_cluster_width1_vv);
    _outtree->Branch("shower3D_cluster_width2_vv" , &_shower3D_cluster_width2_vv);
    
    _outtree->Branch("shower3D_cluster_theta_vv" , &_shower3D_cluster_theta_vv);
    _outtree->Branch("shower3D_cluster_phi_vv"   , &_shower3D_cluster_phi_vv);
    
    _outtree->Branch("shower3D_cluster_opening_vv"   , &_shower3D_cluster_opening_vv);
    _outtree->Branch("shower3D_cluster_opening1_vv"  , &_shower3D_cluster_opening1_vv);
    _outtree->Branch("shower3D_cluster_opening2_vv"  , &_shower3D_cluster_opening2_vv);
    
    _outtree->Branch("shower3D_cluster_distance_vv"  , &_shower3D_cluster_distance_vv);

    if(_fill2d) {
      _outtree->Branch("shower2D_n_clusters_U_v" , &_shower2D_n_clusters_U_v);
      _outtree->Branch("shower2D_n_clusters_V_v" , &_shower2D_n_clusters_V_v);
      _outtree->Branch("shower2D_n_clusters_Y_v" , &_shower2D_n_clusters_Y_v);

      _outtree->Branch("shower2D_area_U_vv"   , &_shower2D_area_U_vv);
      _outtree->Branch("shower2D_length_U_vv" , &_shower2D_length_U_vv);
      _outtree->Branch("shower2D_width_U_vv"  , &_shower2D_width_U_vv);
      _outtree->Branch("shower2D_npixel_U_vv" , &_shower2D_npixel_U_vv);
      _outtree->Branch("shower2D_qsum_U_vv"   , &_shower2D_qsum_U_vv);

      _outtree->Branch("shower2D_area_V_vv"   , &_shower2D_area_V_vv);
      _outtree->Branch("shower2D_length_V_vv" , &_shower2D_length_V_vv);
      _outtree->Branch("shower2D_width_V_vv"  , &_shower2D_width_V_vv);
      _outtree->Branch("shower2D_npixel_V_vv" , &_shower2D_npixel_V_vv);
      _outtree->Branch("shower2D_qsum_V_vv"   , &_shower2D_qsum_V_vv);
    
      _outtree->Branch("shower2D_area_Y_vv"   , &_shower2D_area_Y_vv);
      _outtree->Branch("shower2D_length_Y_vv" , &_shower2D_length_Y_vv);
      _outtree->Branch("shower2D_width_Y_vv"  , &_shower2D_width_Y_vv);
      _outtree->Branch("shower2D_npixel_Y_vv" , &_shower2D_npixel_Y_vv);
      _outtree->Branch("shower2D_qsum_Y_vv"   , &_shower2D_qsum_Y_vv);

      _outtree->Branch("shower2D_n_defects_U_vv"  , &_shower2D_n_defects_U_vv);
      _outtree->Branch("shower2D_n_defects_V_vv"  , &_shower2D_n_defects_V_vv);
      _outtree->Branch("shower2D_n_defects_Y_vv"  , &_shower2D_n_defects_Y_vv);
    }

    LLCV_DEBUG() << "end" << std::endl;
  }

  double SelTrackScatter::Select() {
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

    for(auto const& meta : meta_v)
      _PixelScan3D.SetPlaneInfo(*meta);
    
    //std::array<cv::Mat,3> white_mat_v;
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
      // white_mat_v[plane] = larocv::BlankImage(*mat);
    }
    

    //
    // Project the vertex onto the plane
    //
    GEO2D_Contour_t vertex_ctor;
    vertex_ctor.resize(3);

    for(size_t plane=0; plane<3; ++plane) {
      int px_x = kINVALID_INT;
      int px_y = kINVALID_INT;

      ProjectMat(img_v[plane]->meta(),
    		 Data().Vertex()->X(),Data().Vertex()->Y(),Data().Vertex()->Z(),
    		 px_x, px_y);

      vertex_ctor[plane] = cv::Point_<int>(px_y,px_x);
    }
    

    //
    // Project reconstructed tracks onto plane
    //
    auto track_v = Data().Tracks();
    
    std::vector<std::vector<GEO2D_Contour_t> > track_ctor_vv;
    track_ctor_vv.resize(track_v.size());

    for(size_t tid=0; tid<track_v.size(); ++tid) {
      const auto track = track_v[tid];
      auto& track_ctor_v = track_ctor_vv[tid];
      track_ctor_v.resize(3);

      for(auto& v : track_ctor_v) 
    	v.reserve(track->NumberTrajectoryPoints());

      for(size_t pid=0; pid< track->NumberTrajectoryPoints(); ++pid) {
	
    	const auto& pt = track->LocationAtPoint(pid);
	
    	for(size_t plane=0; plane<3; ++plane) {
	  
    	  auto& track_ctor = track_ctor_v[plane];
	  
    	  int px_x = kINVALID_INT;
    	  int px_y = kINVALID_INT;
	  
    	  ProjectMat(img_v[plane]->meta(),
    		     pt.X(),pt.Y(),pt.Z(),
    		     px_x, px_y);
	  
    	  if (px_x < 0) continue;
    	  if (px_y < 0) continue;

    	  if (px_x >= (int)_cropx) continue;
    	  if (px_y >= (int)_cropy) continue;
	  
	  if(_debug) {
	    auto& mat3d = mat3d_v[plane];
	  
	    if (tid==0) mat3d.at<cv::Vec3b>(px_x,px_y) = {0,0,255};
	    if (tid==1) mat3d.at<cv::Vec3b>(px_x,px_y) = {0,255,0};
	    if (tid==2) mat3d.at<cv::Vec3b>(px_x,px_y) = {255,0,0};
	  }

    	  track_ctor.emplace_back(px_y,px_x);
	  
    	} // end plane
      } // end track point
    } // end track
    

    //
    // generate clusters, associate to tracks
    //
    std::vector<larocv::GEO2D_ContourArray_t> ctor_vv;
    ctor_vv.resize(3);

    std::vector<std::vector<std::vector<size_t> > > tass_vvv;
    tass_vvv.resize(3);

    for(size_t plane=0; plane<3; ++plane) {
      const auto& mat = thresh_mat_v[plane];
      auto& ctor_v = ctor_vv[plane];
      auto& tass_vv= tass_vvv[plane];

      ctor_v  = FindAndMaskVertex(mat,vertex_ctor[plane]);
      tass_vv = AssociateToTracks(ctor_v, track_ctor_vv, plane);
      
      if(!_debug) continue;

      auto& img3d = mat3d_v[plane];
      for(size_t cid=0; cid<ctor_v.size(); ++cid) {
      	cv::Scalar color(std::rand() % 255,std::rand() % 255,std::rand() % 255);
      	if (!ContainsTrack(tass_vv[cid])) continue;
	cv::drawContours(img3d,ctor_v,cid,color);
      }

    }

    //
    // generate 3D consistent points using radial approximation
    // copmuter 3D object paramters
    //
    
    larocv::data::Vertex3D vtx3d;
    vtx3d.x = Data().Vertex()->X();
    vtx3d.y = Data().Vertex()->Y();
    vtx3d.z = Data().Vertex()->Z();
  
    float _track_threshold = 0.5;

    ResizeOutput(track_v.size());

    for(size_t tid=0; tid<track_v.size(); ++tid) {

      if(_debug) {
	auto& track_x_v = _track_x_vv[tid];
	auto& track_y_v = _track_y_vv[tid];
	auto& track_z_v = _track_z_vv[tid];

	for(size_t pid=0; pid< track_v[tid]->NumberTrajectoryPoints(); ++pid) {
	  const auto& pt = track_v[tid]->LocationAtPoint(pid);
	  track_x_v.push_back(pt.X());
	  track_y_v.push_back(pt.Y());
	  track_z_v.push_back(pt.Z());
	}
      }

      LLCV_DEBUG() << "@tid=" << tid << std::endl;
      
      std::array<cv::Mat,3> timg_v;

      for(size_t plane=0; plane<3; ++plane) {
	
    	larocv::GEO2D_ContourArray_t veto_v;

    	const auto& ctor_v = ctor_vv[plane];
    	const auto& tass_vv= tass_vvv[plane];
	
    	for(size_t cid=0; cid<ctor_v.size(); ++cid) {
    	  const auto& tass_v = tass_vv[cid];

    	  if (TrackFraction(tass_v,tid) < _track_threshold) 
    	    continue;

    	  veto_v.emplace_back(ctor_v[cid]);
    	}
	
    	const auto& thresh_mat = thresh_mat_v[plane];
    	timg_v[plane] = larocv::MaskImage(thresh_mat,veto_v,-1,false);
      }
      
      //
      // Estimate the track angle
      //
      auto theta_phi = TrackAngle(*track_v[tid]);
      auto track_len = TrackLength(*track_v[tid]);

      float pi8 = 3.14159 / 8.0;
      float pi4 = 3.14159 / 4.0;
      
      float rad_min = 0.5;
      float rad_max = track_len + 5;
      float rad_step = 0.5;
      
      float theta_min = theta_phi.first - pi8;
      float theta_max = theta_phi.first + pi8;

      float phi_min = theta_phi.second - pi4;
      float phi_max = theta_phi.second + pi4;

      LLCV_DEBUG() << "Scanning Spheres" << std::endl;

      //
      // Build 3D object
      //
      _twatch.Start();

      _PixelScan3D.Reconfigure(rad_min, rad_max, rad_step,
			       theta_min, theta_max,
			       phi_min, phi_max);
      
      
      //auto reg_v = _PixelScan3D.SphereScan3D(white_mat_v,vtx3d,3);
      std::vector<larocv::data::Vertex3D> reg_v;
      std::vector<larocv::data::Vertex3D> skel_v;
      
      if(_allow_dead_image)
	reg_v = _PixelScan3D.SphereScan3D(timg_v,dead_v,vtx3d);
      else
	reg_v = _PixelScan3D.SphereScan3D(timg_v,vtx3d);
      
      _twatch.Stop();
      LLCV_DEBUG() << reg_v.size() << " spheres scanned in " 
		   << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;

      if(_skeletonize and !reg_v.empty()) {
	_twatch.Start();
	_Skeletonize.Initialize(reg_v,0.5);
	_twatch.Stop();
	
	LLCV_DEBUG() << "Initialize skeleton in " 
		     << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;
	
	_twatch.Start();
	skel_v = _Skeletonize.Run();
	_twatch.Stop();
	
	LLCV_DEBUG() << "Skeletonize to " << skel_v.size() << " pts in "
		     << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;
      }

      if(_sub_skeleton)
	reg_v = skel_v;

      if (!reg_v.empty()) {
	LLCV_DEBUG() << "Objectify" << std::endl;
	_twatch.Start();
	Object3D obj(vtx3d,reg_v);
	_twatch.Stop();
	LLCV_DEBUG() << "...objectified in " << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;

	if(_debug) {
	  //
	  // Shower3D start
	  //
	  auto& shower3D_start_x_v = _shower_start_x_vv[tid];
	  auto& shower3D_start_y_v = _shower_start_y_vv[tid];
	  auto& shower3D_start_z_v = _shower_start_z_vv[tid];

	  shower3D_start_x_v.push_back(obj.Start().x);
	  shower3D_start_y_v.push_back(obj.Start().y);
	  shower3D_start_z_v.push_back(obj.Start().z);

	  //
	  // Shower3D Center
	  //
	  auto& shower3D_center_x_v = _shower_center_x_vv[tid];
	  auto& shower3D_center_y_v = _shower_center_y_vv[tid];
	  auto& shower3D_center_z_v = _shower_center_z_vv[tid];

	  shower3D_center_x_v.push_back(obj.Center().x);
	  shower3D_center_y_v.push_back(obj.Center().y);
	  shower3D_center_z_v.push_back(obj.Center().z);

	  //
	  // Shower3D End
	  //
	  auto& shower3D_end_x_v = _shower_end_x_vv[tid];
	  auto& shower3D_end_y_v = _shower_end_y_vv[tid];
	  auto& shower3D_end_z_v = _shower_end_z_vv[tid];

	  shower3D_end_x_v.push_back(obj.End().x);
	  shower3D_end_y_v.push_back(obj.End().y);
	  shower3D_end_z_v.push_back(obj.End().z);

	  //
	  // Shower3D Edge1
	  //
	  auto& shower3D_edge1_x_v = _shower_edge1_x_vv[tid];
	  auto& shower3D_edge1_y_v = _shower_edge1_y_vv[tid];
	  auto& shower3D_edge1_z_v = _shower_edge1_z_vv[tid];

	  shower3D_edge1_x_v.push_back(obj.Edge1().x);
	  shower3D_edge1_y_v.push_back(obj.Edge1().y);
	  shower3D_edge1_z_v.push_back(obj.Edge1().z);

	  //
	  // Shower3D Edge2
	  //
	  auto& shower3D_edge2_x_v = _shower_edge2_x_vv[tid];
	  auto& shower3D_edge2_y_v = _shower_edge2_y_vv[tid];
	  auto& shower3D_edge2_z_v = _shower_edge2_z_vv[tid];
	
	  shower3D_edge2_x_v.push_back(obj.Edge2().x);
	  shower3D_edge2_y_v.push_back(obj.Edge2().y);
	  shower3D_edge2_z_v.push_back(obj.Edge2().z);
	}

	//
	// Shower3D dimensions
	//
	_shower3D_n_points_v[tid] = obj.Points().size();
	_shower3D_length_v[tid]   = obj.Length();
	_shower3D_width_v[tid]    = obj.Width();
	_shower3D_width1_v[tid]   = obj.Width1();
	_shower3D_width2_v[tid]   = obj.Width2();
      
	_shower3D_theta_v[tid] = obj.Theta();
	_shower3D_phi_v[tid]   = obj.Phi();
      
	_shower3D_opening_v[tid]  = obj.Opening();
	_shower3D_opening1_v[tid] = obj.Opening1();
	_shower3D_opening2_v[tid] = obj.Opening2();


	//
	// Compute average deviation from PCA and track
	//
	auto& shower_pca_dev_v = _shower_pca_dev_vv[tid];
	auto& shower_trk_dev_v = _shower_trk_dev_vv[tid];

	shower_pca_dev_v = obj.PCADeviation();
	shower_trk_dev_v = obj.TrackDeviation(*(track_v[tid]));

	auto pca_dev = ComputeMeans(shower_pca_dev_v);
	auto trk_dev = ComputeMeans(shower_trk_dev_v);

	_shower3D_pca_mean_dev_v[tid]        = pca_dev[0];
	_shower3D_start_pca_mean_dev_v[tid]  = pca_dev[1];
	_shower3D_middle_pca_mean_dev_v[tid] = pca_dev[2];
	_shower3D_end_pca_mean_dev_v[tid]    = pca_dev[3];
	
	_shower3D_track_mean_dev_v[tid]        = trk_dev[0];
	_shower3D_start_track_mean_dev_v[tid]  = trk_dev[1];
	_shower3D_middle_track_mean_dev_v[tid] = trk_dev[2];
	_shower3D_end_track_mean_dev_v[tid]    = trk_dev[3];

	//
	// Cluster using DBSCAN
	//
	LLCV_DEBUG() << "Clustering..." << std::endl;

	auto& shower_cid_v = _shower_cid_vv[tid];

	_twatch.Start();
	shower_cid_v = Cluster(obj);
	_twatch.Stop();

	LLCV_DEBUG() << "... clustered in " 
		     << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;

	int n_clusters = CountClusters(shower_cid_v);
	_shower3D_n_clusters_v[tid] = n_clusters;

	LLCV_DEBUG() << "returned " << n_clusters << " clusters" << std::endl;

	auto& shower_x_v = _shower_x_vv[tid];
	auto& shower_y_v = _shower_y_vv[tid];
	auto& shower_z_v = _shower_z_vv[tid];

	std::vector<std::vector<const larocv::data::Vertex3D*> > pts_cluster_vv;
	pts_cluster_vv.resize(n_clusters);

	for(auto& v : pts_cluster_vv) 
	  v.reserve(obj.Points().size());

	auto& shower_skel_x_v = _shower_skel_x_vv[tid];
	auto& shower_skel_y_v = _shower_skel_y_vv[tid];
	auto& shower_skel_z_v = _shower_skel_z_vv[tid];

	for(const auto& skel : skel_v) {
	  shower_skel_x_v.push_back(skel.x);
	  shower_skel_y_v.push_back(skel.y);
	  shower_skel_z_v.push_back(skel.z);
	}

	auto p0_v = _PixelScan3D.ProjectAndDistance(timg_v[0],0,obj.Points(),shower_trk_dev_v);
	auto p1_v = _PixelScan3D.ProjectAndDistance(timg_v[1],1,obj.Points(),shower_trk_dev_v);
	auto p2_v = _PixelScan3D.ProjectAndDistance(timg_v[2],2,obj.Points(),shower_trk_dev_v);

	auto& shower_p0_x_v = _shower_p0_x_vv[tid];
	auto& shower_p0_y_v = _shower_p0_y_vv[tid];
	auto& shower_p0_z_v = _shower_p0_z_vv[tid];
	
	auto& shower_p1_x_v = _shower_p1_x_vv[tid];
	auto& shower_p1_y_v = _shower_p1_y_vv[tid];
	auto& shower_p1_z_v = _shower_p1_z_vv[tid];
	
	auto& shower_p2_x_v = _shower_p2_x_vv[tid];
	auto& shower_p2_y_v = _shower_p2_y_vv[tid];
	auto& shower_p2_z_v = _shower_p2_z_vv[tid];

	for(const auto& p0 : p0_v) {
	  shower_p0_x_v.push_back(p0->x);
	  shower_p0_y_v.push_back(p0->y);
	  shower_p0_z_v.push_back(p0->z);
	}

	for(const auto& p1 : p1_v) {
	  shower_p1_x_v.push_back(p1->x);
	  shower_p1_y_v.push_back(p1->y);
	  shower_p1_z_v.push_back(p1->z);
	}

	for(const auto& p2 : p2_v) {
	  shower_p2_x_v.push_back(p2->x);
	  shower_p2_y_v.push_back(p2->y);
	  shower_p2_z_v.push_back(p2->z);
	}

	//
	// Make the clusters
	//
	int nfail = 0;
	for(size_t pid=0; pid < obj.Points().size(); ++pid) {
	  const auto& pt = obj.Points()[pid];
	  
	  int cid = shower_cid_v[pid];
	  
	  if (cid == -3) {
	    nfail++;
	    continue;
	  }


	  const auto pt_ptr = &pt;
	  assert(pt_ptr);

	  pts_cluster_vv[cid].push_back(pt_ptr);
	  
	  if(!_debug) continue;

	  shower_x_v.push_back(pt.x);
	  shower_y_v.push_back(pt.y);
	  shower_z_v.push_back(pt.z);

	  if (_sub_skeleton) continue;

	  for(size_t plane=0; plane<3; ++plane) {
	    auto& mat3d = mat3d_v[plane];
	    int px_x = pt.vtx2d_v[plane].pt.x;
	    int px_y = pt.vtx2d_v[plane].pt.y;
	    if (px_x < 0) continue;
	    if (px_y < 0) continue;
	    if (px_x >= mat3d.rows) continue;
	    if (px_y >= mat3d.cols) continue;
	    
	    mat3d.at<cv::Vec3b>(px_y,px_x) = {255,255,0};
	  }
	} // end sphere point

	LLCV_DEBUG() << "nfail=(" << nfail << "/" << obj.Points().size() << ")" << std::endl;

	//
	// analyze cluster objects
	//
	auto& shower3D_cluster_n_points_v = _shower3D_cluster_n_points_vv[tid];

	auto& shower3D_cluster_length_v = _shower3D_cluster_length_vv[tid];
	auto& shower3D_cluster_width_v = _shower3D_cluster_width_vv[tid];
	auto& shower3D_cluster_width1_v = _shower3D_cluster_width1_vv[tid];
	auto& shower3D_cluster_width2_v = _shower3D_cluster_width2_vv[tid];
	
	auto& shower3D_cluster_theta_v = _shower3D_cluster_theta_vv[tid];
	auto& shower3D_cluster_phi_v = _shower3D_cluster_phi_vv[tid];
	
	auto& shower3D_cluster_opening_v = _shower3D_cluster_opening_vv[tid];
	auto& shower3D_cluster_opening1_v = _shower3D_cluster_opening1_vv[tid];
	auto& shower3D_cluster_opening2_v = _shower3D_cluster_opening2_vv[tid];
	
	auto& shower3D_cluster_distance_v = _shower3D_cluster_distance_vv[tid];

	shower3D_cluster_n_points_v.resize(n_clusters);

	shower3D_cluster_length_v.resize(n_clusters);
	shower3D_cluster_width_v.resize(n_clusters);
	shower3D_cluster_width1_v.resize(n_clusters);
	shower3D_cluster_width2_v.resize(n_clusters);
	
	shower3D_cluster_theta_v.resize(n_clusters);
	shower3D_cluster_phi_v.resize(n_clusters);
	
	shower3D_cluster_opening_v.resize(n_clusters);
	shower3D_cluster_opening1_v.resize(n_clusters);
	shower3D_cluster_opening2_v.resize(n_clusters);
	
	shower3D_cluster_distance_v.resize(n_clusters);

	for(size_t cid=0; cid<pts_cluster_vv.size(); ++cid) {
	  
	  LLCV_DEBUG() << "@cid=" << cid << std::endl;
	  const auto& pts_cluster_v = pts_cluster_vv[cid];
	  _twatch.Start();
	  Object3D cluster(vtx3d,pts_cluster_v);
	  _twatch.Stop();
	  LLCV_DEBUG() << "... " << std::setprecision(6) << _twatch.RealTime() << "s" << std::endl;

	  shower3D_cluster_n_points_v[cid] = cluster.Points().size();
	  
	  shower3D_cluster_length_v[cid] = cluster.Length();
	  shower3D_cluster_width_v[cid] = cluster.Width();
	  shower3D_cluster_width1_v[cid] = cluster.Width1();
	  shower3D_cluster_width2_v[cid] = cluster.Width2();
	  
	  shower3D_cluster_theta_v[cid] = cluster.Theta();
	  shower3D_cluster_phi_v[cid] = cluster.Phi();
	  
	  shower3D_cluster_opening_v[cid] = cluster.Opening();
	  shower3D_cluster_opening1_v[cid] = cluster.Opening1();
	  shower3D_cluster_opening2_v[cid] = cluster.Opening2();
	}

	if(!_fill2d) continue;
	if(_sub_skeleton) continue;

	//
	// Analyze 2D clusters
	//
      
	// make an image which has pixels only from this track
	std::vector<cv::Mat>  track_mat_v;
	track_mat_v.resize(3);
      
	std::vector<larocv::GEO2D_ContourArray_t> track_ctor_vv;
	track_ctor_vv.resize(3);
      
	for(size_t plane=0; plane<3; ++plane) {
	  auto track_mat = larocv::BlankImage(thresh_mat_v[plane],0);

	  for(const auto& pt : obj.Points()) {
	    int px_x = pt.vtx2d_v[plane].pt.x;
	    int px_y = pt.vtx2d_v[plane].pt.y;
	    if (px_x < 0) continue;
	    track_mat.at<uchar>(px_y,px_x) = (uchar)10;
	  }
	
	  auto& track_ctor_v = track_ctor_vv[plane];

	  track_ctor_v = larocv::FindContours(track_mat);
	  track_mat_v[plane] = track_mat;
	}
      
	//
	// Fill parameters
	//
	auto& shower2D_n_clusters_U = _shower2D_n_clusters_U_v[tid];
	auto& shower2D_n_clusters_V = _shower2D_n_clusters_V_v[tid];
	auto& shower2D_n_clusters_Y = _shower2D_n_clusters_Y_v[tid];

	auto& shower2D_area_U_v    = _shower2D_area_U_vv[tid];
	auto& shower2D_length_U_v  = _shower2D_length_U_vv[tid];
	auto& shower2D_width_U_v   = _shower2D_width_U_vv[tid];
	auto& shower2D_npixel_U_v  = _shower2D_npixel_U_vv[tid];
	auto& shower2D_qsum_U_v    = _shower2D_qsum_U_vv[tid];	

	auto& shower2D_area_V_v    = _shower2D_area_V_vv[tid];
	auto& shower2D_length_V_v  = _shower2D_length_V_vv[tid];
	auto& shower2D_width_V_v   = _shower2D_width_V_vv[tid];
	auto& shower2D_npixel_V_v  = _shower2D_npixel_V_vv[tid];
	auto& shower2D_qsum_V_v    = _shower2D_qsum_V_vv[tid];

	auto& shower2D_area_Y_v    = _shower2D_area_Y_vv[tid];
	auto& shower2D_length_Y_v  = _shower2D_length_Y_vv[tid];
	auto& shower2D_width_Y_v   = _shower2D_width_Y_vv[tid];
	auto& shower2D_npixel_Y_v  = _shower2D_npixel_Y_vv[tid];
	auto& shower2D_qsum_Y_v    = _shower2D_qsum_Y_vv[tid];

	auto& shower2D_n_defects_U_v = _shower2D_n_defects_U_vv[tid];
	auto& shower2D_n_defects_V_v = _shower2D_n_defects_V_vv[tid];
	auto& shower2D_n_defects_Y_v = _shower2D_n_defects_Y_vv[tid];

	for(size_t plane=0; plane<3; ++plane) {
	  const auto& track_mat    = track_mat_v[plane];
	  const auto& track_ctor_v = track_ctor_vv[plane];

	  int nctor = track_ctor_v.size();
	  LLCV_DEBUG() << "a plane=" << plane << " nctor=" << nctor << std::endl;

	  if(plane==0) { 
	    shower2D_n_clusters_U = nctor;

	    shower2D_area_U_v.resize(nctor);	    
	    shower2D_length_U_v.resize(nctor);
	    shower2D_width_U_v.resize(nctor);
	    shower2D_npixel_U_v.resize(nctor);
	    shower2D_qsum_U_v.resize(nctor);

	    for(int cid=0; cid<nctor; ++cid) {
	      const auto& ctor = track_ctor_v[cid];
	      larocv::PixelChunk pc(ctor,*(mat_v[plane]));
	      shower2D_area_U_v[cid]   = pc.area;
	      shower2D_length_U_v[cid] = pc.length;
	      shower2D_width_U_v[cid]  = pc.width ;
	      shower2D_npixel_U_v[cid] = pc.npixel;
	      shower2D_qsum_U_v[cid]   = pc.qsum;
	    }

	  }

	  if(plane==1) { 
	    shower2D_n_clusters_V = nctor;

	    shower2D_area_V_v.resize(nctor);	    
	    shower2D_length_V_v.resize(nctor);
	    shower2D_width_V_v.resize(nctor);
	    shower2D_npixel_V_v.resize(nctor);
	    shower2D_qsum_V_v.resize(nctor);

	    for(int cid=0; cid<nctor; ++cid) {
	      const auto& ctor = track_ctor_v[cid];
	      larocv::PixelChunk pc(ctor,*(mat_v[plane]));
	      shower2D_area_V_v[cid]   = pc.area;
	      shower2D_length_V_v[cid] = pc.length;
	      shower2D_width_V_v[cid]  = pc.width ;
	      shower2D_npixel_V_v[cid] = pc.npixel;
	      shower2D_qsum_V_v[cid]   = pc.qsum;
	    }

	  }

	  if(plane==2) { 
	    shower2D_n_clusters_Y = nctor;	
	
	    shower2D_area_Y_v.resize(nctor);	    
	    shower2D_length_Y_v.resize(nctor);
	    shower2D_width_Y_v.resize(nctor);
	    shower2D_npixel_Y_v.resize(nctor);
	    shower2D_qsum_Y_v.resize(nctor);

	    for(int cid=0; cid<nctor; ++cid) {
	      const auto& ctor = track_ctor_v[cid];
	      larocv::PixelChunk pc(ctor,*(mat_v[plane]));
	      shower2D_area_Y_v[cid]   = pc.area;
	      shower2D_length_Y_v[cid] = pc.length;
	      shower2D_width_Y_v[cid]  = pc.width ;
	      shower2D_npixel_Y_v[cid] = pc.npixel;
	      shower2D_qsum_Y_v[cid]   = pc.qsum;
	    }
	  }
	} // end plane

      } //end non-empty scan

      LLCV_DEBUG() << "done!" << std::endl;
    } // end this track
    

    if(_debug) {
      for(size_t plane=0; plane<3; ++plane) {
	auto& img3d = mat3d_v[plane];
	std::stringstream ss; 
	ss << "png/plane_img_" << Run() << "_" << SubRun() << "_" << Event() << "_" << VertexID() << "_" << plane << ".png";
	cv::imwrite(ss.str(),img3d);
      }
    }

    _outtree->Fill();
    
    LLCV_DEBUG() << "end" << std::endl;
    return 0.0;
  }


  larocv::GEO2D_ContourArray_t SelTrackScatter::FindAndMaskVertex(const cv::Mat& mat,
								  const cv::Point_<int> vertex) {
    larocv::GEO2D_ContourArray_t ctor_v;

    auto vtx_circ = geo2d::Circle<float>(vertex.x, vertex.y, 4);
    auto mat_mask = larocv::MaskImage(mat,vtx_circ);

    auto ctor_tmp_v = larocv::FindContours(mat_mask);

    // Split off the small contours
    larocv::GEO2D_ContourArray_t small_ctor_v, large_ctor_v;
    small_ctor_v.reserve((int) ((float)ctor_tmp_v.size() / 2.0));
    large_ctor_v.reserve((int) ((float)ctor_tmp_v.size() / 2.0));
    
    for(auto& ctor_tmp : ctor_tmp_v) {
      if (ctor_tmp.size()<=2) 
	small_ctor_v.emplace_back(std::move(ctor_tmp));
      else 
	large_ctor_v.emplace_back(std::move(ctor_tmp));
    }
    
    std::swap(ctor_tmp_v,large_ctor_v);

    larocv::GEO2D_ContourArray_t ctor_tmp_split_v;
    ctor_tmp_split_v.reserve(ctor_tmp_v.size());

    for(size_t cid=0; cid< ctor_tmp_v.size(); ++cid) {
      auto split_v = _LR_DefectBreaker.SplitContour(ctor_tmp_v[cid]);

      for(auto& split : split_v)
	ctor_tmp_split_v.emplace_back(std::move(split));
    }
    
    size_t sz_estimate = ctor_tmp_split_v.size() + small_ctor_v.size();
    ctor_v.reserve(sz_estimate);

    for(auto& ctor_tmp_split : ctor_tmp_split_v)
      ctor_v.emplace_back(std::move(ctor_tmp_split));

    for(auto& small_ctor : small_ctor_v)
      ctor_v.emplace_back(std::move(small_ctor));    
    
    return ctor_v;
  }


  larocv::GEO2D_ContourArray_t SelTrackScatter::FindAndBreakVertex(const cv::Mat& mat,
								   const cv::Point_<int> vertex) {

    larocv::GEO2D_ContourArray_t ctor_v;
    
    auto ctor_tmp_v = larocv::FindContours(mat);

    // Split off the small contours
    larocv::GEO2D_ContourArray_t small_ctor_v, large_ctor_v;
    small_ctor_v.reserve((int) ((float)ctor_tmp_v.size() / 2.0));
    large_ctor_v.reserve((int) ((float)ctor_tmp_v.size() / 2.0));
    
    for(auto& ctor_tmp : ctor_tmp_v) {
      if (ctor_tmp.size()<=2) 
	small_ctor_v.emplace_back(std::move(ctor_tmp));
      else 
	large_ctor_v.emplace_back(std::move(ctor_tmp));
    }
      
    std::swap(ctor_tmp_v,large_ctor_v);

    // Find the contour which contains the vertex
    auto vid = larocv::FindContainingContour(ctor_tmp_v,vertex);
    if (vid==kINVALID_SIZE) 
      throw llcv_err("Could not find containing contour");

    // Break this contour -- long range
    auto blr_ctor_tmp_v = _LR_DefectBreaker.SplitContour(ctor_tmp_v[vid]);

    // Find the broken contour with vertex
    auto blr_vid = larocv::FindContainingContour(blr_ctor_tmp_v,vertex);
    if (blr_vid==kINVALID_SIZE) 
      throw llcv_err("Could not find containing contour");

    // Break this contour -- short range
    auto bsr_ctor_tmp_v = _SR_DefectBreaker.SplitContour(blr_ctor_tmp_v[blr_vid]);

    size_t sz_estimate = ctor_tmp_v.size() + blr_ctor_tmp_v.size() + bsr_ctor_tmp_v.size() + small_ctor_v.size();
    ctor_v.reserve(sz_estimate);

    for(size_t cid = 0; cid < ctor_tmp_v.size(); ++cid) {
      if (cid == vid) continue;
      ctor_v.emplace_back(std::move(ctor_tmp_v[cid]));
    }

    for(size_t cid = 0; cid < blr_ctor_tmp_v.size(); ++cid) {
      if (cid == blr_vid) continue;
      ctor_v.emplace_back(std::move(blr_ctor_tmp_v[cid]));
    }
      
    for(auto& bsr_ctor_tmp : bsr_ctor_tmp_v)
      ctor_v.emplace_back(std::move(bsr_ctor_tmp));

    for(auto& small_ctor : small_ctor_v)
      ctor_v.emplace_back(std::move(small_ctor));

    return ctor_v;

  }
  

  std::vector<std::vector<size_t> > SelTrackScatter::AssociateToTracks(const larocv::GEO2D_ContourArray_t& ctor_v,
								       const std::vector<larocv::GEO2D_ContourArray_t>& track_ctor_vv,
								       const size_t plane) {
    LLCV_DEBUG() << "start" << std::endl;
    // Save the contours where projected pixel is inside
    std::vector<std::vector<size_t> > tin_vv;
    tin_vv.resize(ctor_v.size());
    for(auto& v : tin_vv) 
      v.resize(track_ctor_vv.size(),0);
    
    for(size_t cid = 0; cid < ctor_v.size(); ++cid) {
      const auto& ctor = ctor_v[cid];
      auto& tin_v = tin_vv[cid];

      for(size_t tid = 0; tid < track_ctor_vv.size(); ++tid) {
	const auto& track_ctor = track_ctor_vv[tid][plane];

	for(const auto& pt : track_ctor) {
	  double dist;
	  auto inside = larocv::PointPolygonTest(ctor,pt,dist);	    
	  if (inside or dist<=1) {
	    tin_v[tid] += 1;
	  }
	}
      } // end track
    } // end contour
    
    LLCV_DEBUG() << "end" << std::endl;
    return tin_vv;
  }

  std::pair<float,float> SelTrackScatter::TrackAngle(const larlite::track& track) {
    
    std::pair<float,float> res;
    
    const auto& start_pt= track.Vertex();
    const auto& end_pt  = track.End();
    auto dir = end_pt - start_pt;

    res.first  = std::acos( dir.Z() / dir.Mag());
    res.second = std::atan2(dir.X(),dir.Y()); // wtf
    
    return res;
  }

  float SelTrackScatter::TrackLength(const larlite::track& track) {
    
    float res;
    
    const auto& start_pt= track.Vertex();
    const auto& end_pt  = track.End();
    auto dir = end_pt - start_pt;

    res = dir.Mag();

    return res;
  }

  bool SelTrackScatter::ContainsTrack(const std::vector<size_t>& tin_v) {
    
    for(const auto tin : tin_v) {
      if (tin > 0) return true;
    }

    return false;
  }

  float SelTrackScatter::TrackFraction(const std::vector<size_t>& tin_v, size_t tid) {
    
    float ret = 0;
    float sum = 0;
    for(const auto tin : tin_v)
      sum += (float) tin;
    
    if (sum==0) return ret;
    
    ret = ((float)tin_v[tid]) / sum;
    
    return ret;
  }
  
  void SelTrackScatter::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _outtree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }

  std::vector<int> SelTrackScatter::Cluster(const Object3D& obj) {
    std::vector<int> res_v;
    res_v.resize(obj.Points().size(),-1);
    
    static std::vector<Point> pts_v;
    pts_v.clear();
    pts_v.resize(obj.Points().size());
    
    for(size_t pid=0; pid<obj.Points().size(); ++pid) {
      pts_v[pid].x = obj.Points()[pid].x;
      pts_v[pid].y = obj.Points()[pid].y;
      pts_v[pid].z = obj.Points()[pid].z;
    }
    
    _DBSCAN.Reset(pts_v);
    _DBSCAN.run();

    for(size_t pid=0; pid<obj.Points().size(); ++pid) 
      res_v[pid] = _DBSCAN.Points()[pid].clusterID - 1;

    return res_v;
  }


  std::array<float,4> SelTrackScatter::ComputeMeans(const std::vector<float>& data_v) {

    std::array<float,4> res_v = {{-1,-1,-1,-1}};

    float npts = (float)data_v.size();
    
    if (data_v.size() <= 3) {
      for(size_t did=0; did<data_v.size(); ++did) 
	{ res_v[did] = data_v[did]; }
      return res_v;
    }

    std::array<size_t,4> range_v;
    range_v[0] = 0;
    range_v[1] = npts/3;
    range_v[2] = 2*npts/3;
    range_v[3] = npts;
    

    res_v[0] = Average(data_v,0,data_v.size());

    for(size_t rid=0; rid<3; ++rid) 
      res_v[rid+1] = Average(data_v,range_v[rid],range_v[rid+1]);
    
    return res_v;
  }

  float SelTrackScatter::Average(const std::vector<float>& data_v, size_t start, size_t end) {
    float res = 0.0;
    
    for(size_t id=start; id<end; ++id)
      res += data_v[id];
    
    res /= ((float) (end - start));
    
    return res;
  }


  int SelTrackScatter::CountClusters(const std::vector<int>& cid_v) {
    int res = 0;
    static std::unordered_set<int> s;
    s.clear();

    for(auto cid : cid_v) {
      if (cid == -3) continue;
      s.insert(cid);
    }
    
    res = s.size();

    return res;
  }


  void SelTrackScatter::ResizeOutput(size_t sz) {
    
    _track_x_vv.clear();
    _track_y_vv.clear();
    _track_z_vv.clear();

    _shower_x_vv.clear();
    _shower_y_vv.clear();
    _shower_z_vv.clear();

    _shower_skel_x_vv.clear();
    _shower_skel_y_vv.clear();
    _shower_skel_z_vv.clear();

    _shower_p0_x_vv.clear();
    _shower_p0_y_vv.clear();
    _shower_p0_z_vv.clear();

    _shower_p1_x_vv.clear();
    _shower_p1_y_vv.clear();
    _shower_p1_z_vv.clear();

    _shower_p2_x_vv.clear();
    _shower_p2_y_vv.clear();
    _shower_p2_z_vv.clear();

    _shower_start_x_vv.clear();
    _shower_start_y_vv.clear();
    _shower_start_z_vv.clear();

    _shower_end_x_vv.clear();
    _shower_end_y_vv.clear();
    _shower_end_z_vv.clear();

    _shower_center_x_vv.clear();
    _shower_center_y_vv.clear();
    _shower_center_z_vv.clear();

    _shower_cid_vv.clear();
    
    _shower_pca_dev_vv.clear();
    _shower_trk_dev_vv.clear();

    _shower_edge1_x_vv.clear();
    _shower_edge1_y_vv.clear();
    _shower_edge1_z_vv.clear();

    _shower_edge2_x_vv.clear();
    _shower_edge2_y_vv.clear();
    _shower_edge2_z_vv.clear();

    _shower3D_n_points_v.clear();
    
    _shower3D_length_v.clear();
    _shower3D_width_v.clear();
    _shower3D_width1_v.clear();
    _shower3D_width2_v.clear();

    _shower3D_theta_v.clear();
    _shower3D_phi_v.clear();

    _shower3D_opening_v.clear();
    _shower3D_opening1_v.clear();
    _shower3D_opening2_v.clear();

    _shower3D_start_pca_mean_dev_v.clear();
    _shower3D_middle_pca_mean_dev_v.clear();
    _shower3D_end_pca_mean_dev_v.clear();

    _shower3D_start_track_mean_dev_v.clear();
    _shower3D_middle_track_mean_dev_v.clear();
    _shower3D_end_track_mean_dev_v.clear();

    _shower3D_n_clusters_v.clear();

    _shower3D_cluster_n_points_vv.clear();

    _shower3D_cluster_length_vv.clear();
    _shower3D_cluster_width_vv.clear();
    _shower3D_cluster_width1_vv.clear();
    _shower3D_cluster_width2_vv.clear();
    
    _shower3D_cluster_theta_vv.clear();
    _shower3D_cluster_phi_vv.clear();
    
    _shower3D_cluster_opening_vv.clear();
    _shower3D_cluster_opening1_vv.clear();
    _shower3D_cluster_opening2_vv.clear();

    _shower3D_cluster_distance_vv.clear();
    
    _shower2D_n_clusters_U_v.clear();
    _shower2D_n_clusters_V_v.clear();
    _shower2D_n_clusters_Y_v.clear();
    
    _shower2D_area_U_vv.clear();
    _shower2D_length_U_vv.clear();
    _shower2D_width_U_vv.clear();
    _shower2D_npixel_U_vv.clear();
    _shower2D_qsum_U_vv.clear();

    _shower2D_area_V_vv.clear();
    _shower2D_length_V_vv.clear();
    _shower2D_width_V_vv.clear();
    _shower2D_npixel_V_vv.clear();
    _shower2D_qsum_V_vv.clear();

    _shower2D_area_Y_vv.clear();
    _shower2D_length_Y_vv.clear();
    _shower2D_width_Y_vv.clear();
    _shower2D_npixel_Y_vv.clear();
    _shower2D_qsum_Y_vv.clear();

    _shower2D_n_defects_U_vv.clear();
    _shower2D_n_defects_V_vv.clear();
    _shower2D_n_defects_Y_vv.clear();
    
    _track_x_vv.resize(sz);
    _track_y_vv.resize(sz);
    _track_z_vv.resize(sz);

    _shower_x_vv.resize(sz);
    _shower_y_vv.resize(sz);
    _shower_z_vv.resize(sz);

    _shower_p0_x_vv.resize(sz);
    _shower_p0_y_vv.resize(sz);
    _shower_p0_z_vv.resize(sz);

    _shower_p1_x_vv.resize(sz);
    _shower_p1_y_vv.resize(sz);
    _shower_p1_z_vv.resize(sz);

    _shower_p2_x_vv.resize(sz);
    _shower_p2_y_vv.resize(sz);
    _shower_p2_z_vv.resize(sz);

    _shower_skel_x_vv.resize(sz);
    _shower_skel_y_vv.resize(sz);
    _shower_skel_z_vv.resize(sz);

    _shower_start_x_vv.resize(sz);
    _shower_start_y_vv.resize(sz);
    _shower_start_z_vv.resize(sz);

    _shower_end_x_vv.resize(sz);
    _shower_end_y_vv.resize(sz);
    _shower_end_z_vv.resize(sz);

    _shower_center_x_vv.resize(sz);
    _shower_center_y_vv.resize(sz);
    _shower_center_z_vv.resize(sz);

    _shower_cid_vv.resize(sz);
    
    _shower_pca_dev_vv.resize(sz);
    _shower_trk_dev_vv.resize(sz);

    _shower_edge1_x_vv.resize(sz);
    _shower_edge1_y_vv.resize(sz);
    _shower_edge1_z_vv.resize(sz);

    _shower_edge2_x_vv.resize(sz);
    _shower_edge2_y_vv.resize(sz);
    _shower_edge2_z_vv.resize(sz);

    _shower3D_n_points_v.resize(sz);

    _shower3D_length_v.resize(sz);
    _shower3D_width_v.resize(sz);
    _shower3D_width1_v.resize(sz);
    _shower3D_width2_v.resize(sz);

    _shower3D_theta_v.resize(sz);
    _shower3D_phi_v.resize(sz);

    _shower3D_opening_v.resize(sz);
    _shower3D_opening1_v.resize(sz);
    _shower3D_opening2_v.resize(sz);

    _shower3D_pca_mean_dev_v.resize(sz);
    _shower3D_start_pca_mean_dev_v.resize(sz);
    _shower3D_middle_pca_mean_dev_v.resize(sz);
    _shower3D_end_pca_mean_dev_v.resize(sz);

    _shower3D_track_mean_dev_v.resize(sz);
    _shower3D_start_track_mean_dev_v.resize(sz);
    _shower3D_middle_track_mean_dev_v.resize(sz);
    _shower3D_end_track_mean_dev_v.resize(sz);

    _shower3D_n_clusters_v.resize(sz);

    _shower3D_cluster_n_points_vv.resize(sz);

    _shower3D_cluster_length_vv.resize(sz);
    _shower3D_cluster_width_vv.resize(sz);
    _shower3D_cluster_width1_vv.resize(sz);
    _shower3D_cluster_width2_vv.resize(sz);
    
    _shower3D_cluster_theta_vv.resize(sz);
    _shower3D_cluster_phi_vv.resize(sz);
    
    _shower3D_cluster_opening_vv.resize(sz);
    _shower3D_cluster_opening1_vv.resize(sz);
    _shower3D_cluster_opening2_vv.resize(sz);

    _shower3D_cluster_distance_vv.resize(sz);
    
    _shower2D_n_clusters_U_v.resize(sz);
    _shower2D_n_clusters_V_v.resize(sz);
    _shower2D_n_clusters_Y_v.resize(sz);

    _shower2D_area_U_vv.resize(sz);
    _shower2D_length_U_vv.resize(sz);
    _shower2D_width_U_vv.resize(sz);
    _shower2D_npixel_U_vv.resize(sz);
    _shower2D_qsum_U_vv.resize(sz);

    _shower2D_area_V_vv.resize(sz);
    _shower2D_length_V_vv.resize(sz);
    _shower2D_width_V_vv.resize(sz);
    _shower2D_npixel_V_vv.resize(sz);
    _shower2D_qsum_V_vv.resize(sz);

    _shower2D_area_Y_vv.resize(sz);
    _shower2D_length_Y_vv.resize(sz);
    _shower2D_width_Y_vv.resize(sz);
    _shower2D_npixel_Y_vv.resize(sz);
    _shower2D_qsum_Y_vv.resize(sz);
    
    _shower2D_n_defects_U_vv.resize(sz);
    _shower2D_n_defects_V_vv.resize(sz);
    _shower2D_n_defects_Y_vv.resize(sz);
  }


}


#endif

