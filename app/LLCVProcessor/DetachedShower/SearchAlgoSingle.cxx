#ifndef __SEARCHALGOSINGLE_CXX__
#define __SEARCHALGOSINGLE_CXX__

#include "SearchAlgoSingle.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoData/AlgoDataUtils.h"

namespace llcv {
  
  void SearchAlgoSingle::Configure(const larcv::PSet &pset) {

    larocv::logger::get_shared().set((larocv::msg::Level_t)2);
    _shower_frac   = pset.get<float>("ShowerFrac",0.8);
    _shower_size   = pset.get<float>("ShowerSize",10);
    _mask_vertex   = pset.get<bool>("MaskVertex");
    _PixelScan3D.set_verbosity((larocv::msg::Level_t)this->logger().level());
    _PixelScan3D.Configure();

    return;
  }

  std::vector<llcv::DetachedCandidate>
  SearchAlgoSingle::_Search_(const larocv::data::Vertex3D& vtx3d,
			     std::vector<cv::Mat>& adc_mat_v,
			     std::vector<cv::Mat>& shr_mat_v,
			     const std::vector<larocv::ImageMeta>& meta_v) { 

    LLCV_DEBUG() << "start" << std::endl;
    
    //
    // setup the return
    // 
    std::vector<llcv::DetachedCandidate> res_v;
    res_v.clear();


    //
    // threshold the image
    //
    for(auto& img : adc_mat_v)
      img = larocv::Threshold(img,10,255);    

    for(auto& img : shr_mat_v)
      img = larocv::Threshold(img,10,255);    

    for(auto const& meta : meta_v)
      _PixelScan3D.SetPlaneInfo(meta);

    //
    // find adc contours per plane
    //
    std::vector<larocv::GEO2D_ContourArray_t> ctor_vv;
    ctor_vv.resize(3);

    for(size_t plane=0; plane<3; ++plane)
      ctor_vv[plane] = larocv::FindContours(adc_mat_v.at(plane)); 

    
    // shower image which will mask small cluster
    std::array<cv::Mat,3> simg_v; 

    // id of contour which is closest to vertex ( 10 px threshold )
    std::array<size_t, 3> vtx_ctor_v; 

    // contours in the shower image
    std::array<larocv::GEO2D_ContourArray_t,3> sctor_vv;

    // candidate contours to perform 3D scan
    std::array<larocv::GEO2D_ContourArray_t,3> actor_vv;

    LLCV_DEBUG() << "@vtx3d=(" << vtx3d.x << "," << vtx3d.y << "," << vtx3d.z << ")" << std::endl;
    
    //
    // find the contour closest to the vertex -- may or may not exist
    //
    for(size_t plane=0; plane<3; ++plane) {
      LLCV_DEBUG() << "@plane=" << plane << std::endl;
      
      simg_v[plane] = shr_mat_v[plane].clone();

      const auto& ctor_v = ctor_vv.at(plane);
      const auto& vtx2d = vtx3d.vtx2d_v.at(plane);
      
      LLCV_DEBUG() << "2d pt=" << vtx2d.pt << std::endl;
      
      if (vtx2d.pt.x < 0 or vtx2d.pt.x == kINVALID_FLOAT) continue;
      if (vtx2d.pt.y < 0 or vtx2d.pt.y == kINVALID_FLOAT) continue;
      
      double distance = kINVALID_DOUBLE;
      auto id = larocv::FindContainingContour(ctor_v, vtx2d.pt,distance);

      if (id == kINVALID_SIZE) {
	LLCV_DEBUG() << "Could not be found..." << std::endl;
	vtx_ctor_v[plane] = kINVALID_SIZE;
	continue;      
      }

      if (distance < -10) {
	LLCV_DEBUG() << "Too far away..." << std::endl;
	vtx_ctor_v[plane] = kINVALID_SIZE;
	continue;
      }

      vtx_ctor_v[plane] = id;
      
      LLCV_DEBUG() << "Identified vertex contour @ id=" << id << std::endl;
    }

    //
    // find shower contours
    //
    for(size_t plane=0; plane<3; ++plane)
      sctor_vv[plane] = larocv::FindContours(simg_v[plane]); 
    
    //
    // clear possible detached shower contours
    //
    for(auto& v : actor_vv) v.clear();
    
    //
    // filter contours by size, area, shower fraction
    //
    for(size_t plane = 0; plane<3; ++plane) {

      const auto& sctor_v  = sctor_vv[plane];
      const auto& ctor_v   = ctor_vv[plane];

      auto& actor_v = actor_vv[plane];
      actor_v.reserve(ctor_v.size());
      
      // loop over ADC contours
      for(size_t aid = 0 ; aid < ctor_v.size(); ++aid) {
	LLCV_DEBUG() << "@aid=" << aid << std::endl;
	const auto& ctor = ctor_v[aid];

	// it's the vertex contour
	if ( (aid == vtx_ctor_v[plane]) && (_mask_vertex))  {
	  LLCV_DEBUG() << "skip vtx contour @aid=" << aid << std::endl;
	  simg_v[plane] = larocv::MaskImage(simg_v[plane],ctor,0,true);
	  continue; 
	}

	// it's too small to be a second shower 
	if(larocv::ContourArea(ctor) < _shower_size) {
	  LLCV_DEBUG() << "too small" << std::endl;
	  simg_v[plane] = larocv::MaskImage(simg_v[plane],ctor,0,true);
	  continue; 
	}
	
	// check if this adc contours is in the shower image
	auto sid = larocv::FindContainingContour(sctor_v,ctor);
	if (sid == kINVALID_SIZE) { 
	  LLCV_DEBUG() << "no simg correspond" << std::endl;
	  continue;
	}

	// it's not shower enough
	double frac = larocv::AreaRatio(sctor_v.at(sid),ctor);
	if (frac < _shower_frac)  {
	  LLCV_DEBUG() << "not enough shower frac" << std::endl;
	  simg_v[plane] = larocv::MaskImage(simg_v[plane],ctor,0,true); 
	  continue;
	}

	// candidate shower cluster
	actor_v.emplace_back(ctor);
      }
    }

    //
    // not enough contours, continue
    //
    size_t cnt=0;
    for(size_t plane=0; plane<3; ++plane) {
      const auto& actor_v = actor_vv[plane];
      if (!actor_v.empty()) cnt +=1; 
      LLCV_DEBUG() << "@plane=" << plane << " actor sz=" << actor_v.size() << std::endl;
    }

    //
    // There must be at least 1 on two planes to match
    //
    if (cnt<2) {
      LLCV_DEBUG() << "... no candidate" << std::endl;
      return res_v;
    }
    
    //
    // Register regions for the scan -- spheres of fixed dimension
    //
    LLCV_DEBUG() << "Scan" << std::endl;
    auto reg_vv = _PixelScan3D.RegionScan3D(simg_v,vtx3d);
      
    //
    // Associate the candidate contours to valid 3D points on sphere
    //
    LLCV_DEBUG() << "Associate" << std::endl;
    auto ass_vv = _PixelScan3D.AssociateContours(reg_vv,actor_vv);
      
    //
    // count the number of consistent 3D points per contour ID
    //

    // vector of (trip)lets -- each triplet is a match of contours acrossp planes
    std::vector<std::array<size_t,3>> trip_v;
    
    // counter per triplet of number of consistent 3D points
    std::vector<size_t> trip_cnt_v;
    
    // vector of 3D points per triplet to play with
    std::vector<std::vector<const larocv::data::Vertex3D*> > trip_vtx_ptr_vv;


    bool found = false;
    for(size_t assvid =0; assvid < ass_vv.size(); ++assvid) {
      const auto& ass_v = ass_vv[assvid];
      for(size_t assid =0; assid < ass_v.size(); ++assid) {
	const auto& ass = ass_v[assid];
	const auto& reg = reg_vv.at(assvid).at(assid);

	found = false;
	for(size_t tid=0; tid<trip_v.size(); ++tid) {
	  if(!CompareAsses(ass,trip_v[tid])) continue;
	  found = true;
	  trip_cnt_v[tid] += 1;
	  trip_vtx_ptr_vv[tid].emplace_back(&reg);
	  break;
	}

	if (found) continue;
	trip_v.emplace_back(ass);
	trip_cnt_v.emplace_back(1);
	trip_vtx_ptr_vv.emplace_back(std::vector<const larocv::data::Vertex3D*>(1,&reg));
      }
    }

    LLCV_DEBUG() << "trip_v sz=    " << trip_v.size() << std::endl;
    LLCV_DEBUG() << "trip_cnt_v sz=" << trip_cnt_v.size() << std::endl;
    for(size_t tid=0; tid< trip_v.size(); ++tid) {
      LLCV_DEBUG() << "@tid=" << tid 
		   << " {" << trip_v[tid][0] << "," << trip_v[tid][1]  << "," << trip_v[tid][2] << "} = " 
		   << trip_cnt_v[tid] << std::endl;
    }
      
    
    //
    // no candidate could be associated across planes
    //
    if (trip_v.empty()) {
      LLCV_DEBUG() << "... no associated particles across planes" << std::endl;
      return res_v;
    }

    //
    // pick the most 3D consistent triplet
    //
    size_t maxid = std::distance(trip_cnt_v.begin(), std::max_element(trip_cnt_v.begin(), trip_cnt_v.end()));
    LLCV_DEBUG() << "Selected max element @pos=" << maxid << std::endl;
    const auto& trip           = trip_v.at(maxid);
    const auto& trip_vtx_ptr_v = trip_vtx_ptr_vv.at(maxid);
    
    //
    // make a vector of space points (was to do PCA, but forget it)
    //
    std::vector<larocv::data::SpacePt> sps_v;
    sps_v.resize(trip_vtx_ptr_v.size());
    for(size_t sid=0; sid<sps_v.size(); ++sid) {
      auto& sps = sps_v[sid];
      const auto vtx = trip_vtx_ptr_v[sid];
      sps.pt = *vtx;
    }

    //
    // get the radial point closest to the vertex in 3D to call a start point
    //
    double min_dist = kINVALID_DOUBLE;
    const larocv::data::SpacePt* min_sp = nullptr;
    for(const auto& sp : sps_v) {
      auto d = larocv::data::Distance(sp,vtx3d);
      if (d<min_dist) {
	min_dist = d;
	min_sp = &sp;
      }
    }
    
    LLCV_DEBUG() << "Second shower candidate idenfitied" << std::endl;
    

    //
    // make a single shower since that's what this module does
    //
    res_v.resize(1);
    auto& res = res_v.front();
    
    // set the start point
    res.origin = vtx3d;
    if(min_sp) res.start  = min_sp->pt;
    
    // fill the output
    for(int plane=0; plane<3; ++plane) {
      LLCV_DEBUG() << "@plane=" << plane << std::endl;
      DetachedCluster dc;
      auto aid = trip.at(plane);
      if (aid==kINVALID_SIZE) {
	LLCV_DEBUG() << "...skip!" << std::endl;
	continue;
      }
      dc.ctor = actor_vv.at(plane).at(aid);
      dc.start_x = kINVALID_FLOAT;
      dc.start_y = kINVALID_FLOAT;
      dc.plane = plane;
      res.Move(std::move(dc),plane);
    }
    
    LLCV_DEBUG() << "end" << std::endl;
    // send the output back to the driver to make pgraph + pixel2d
    return res_v;
  }

}
#endif
