#ifndef __SELCOSMICID_CXX__
#define __SELCOSMICID_CXX__

#include "SelCosmicID.h"

#include "InterTool_Util/InterImageUtils.h"

#include "LArOpenCV/ImageCluster/AlgoFunction/Contour2DAnalysis.h"
#include "LArOpenCV/ImageCluster/AlgoFunction/ImagePatchAnalysis.h"

#include <cassert>
#include <cstdlib>
#include <sstream>

namespace llcv {

  void SelCosmicID::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));

    LLCV_DEBUG() << "start" << std::endl;

    _max_radius = pset.get<float>("MaxRadius");
    _min_radius = pset.get<float>("MinRadius");
    _radius_step= pset.get<float>("RadiusStep");

    assert(_max_radius>_min_radius);
    
    _radius_v.reserve((_max_radius-_min_radius) / _radius_step);
    for(float rad=_min_radius; rad<_max_radius; rad+=_radius_step)
      _radius_v.push_back(rad);

    ResetTree();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelCosmicID::Initialize() {
    _fout->cd();
    _out_tree = new TTree("CosmicID","");

    AttachRSEV(_out_tree);
    
    _out_tree->Branch("pt_xing_vv", &_pt_xing_vv);
    _out_tree->Branch("connected_vv", &_connected_vv);
    _out_tree->Branch("distance_vv",&_distance_vv);

    return;
  }

  double SelCosmicID::Select() {
    llcv::logger::get_shared().set((msg::Level_t)2);

    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "=======================" << std::endl;
    
    size_t cropx = (size_t) (_max_radius * 1.5);
    size_t cropy = cropx;

    auto mat_v = Image().Image<cv::Mat>(kImageADC,cropx,cropy);
    auto img_v = Image().Image<larcv::Image2D>(kImageADC,cropx,cropy);
    
    const auto& vertex = Data().Vertex();

    for(size_t plane=0; plane<3; ++plane) {

      LLCV_DEBUG() << "@plane=" << plane << std::endl;

      const auto& meta = img_v.at(plane)->meta();
      int row = kINVALID_INT;
      int col = kINVALID_INT;
      ProjectMat(meta,
		 vertex->X(),vertex->Y(),vertex->Z(),
		 row, col);
      
      std::stringstream ss;

      // ss.str("");
      // ss << VertexID() << "plane_" << plane << "_crop.png";
      // cv::imwrite(ss.str().c_str(),*(mat_v.at(plane)));

      auto mat = larocv::Threshold(*(mat_v.at(plane)),10,255);
      
      // ss.str("");
      // ss << VertexID() << "plane_" << plane << "_crop_thresh.png";
      // cv::imwrite(ss.str().c_str(),mat);

      auto ctor_v = larocv::FindContours(mat);

      geo2d::Vector<float> pt(row,col);

      LLCV_DEBUG() << "@pt=" << pt << std::endl;

      auto ctor_id = larocv::FindContainingContour(ctor_v,pt);
      if (ctor_id == kINVALID_SIZE) {
	LLCV_DEBUG() << "... point not associated with any contour" << std::endl;
	continue;
      }
      
      
      LLCV_DEBUG() << "@ctor_id=" << ctor_id << std::endl;
      
      const auto& ctor = ctor_v.at(ctor_id);

      mat = larocv::MaskImage(mat,ctor,0,false);
	 
      // ss.str("");
      // ss << VertexID() << "plane_" << plane << "_crop_thresh_mask.png";
      // cv::imwrite(ss.str().c_str(),mat);

      auto temp_xs_vv = larocv::OnCircleGroupsOnCircleArray(mat,pt,_radius_v);
      
      for (size_t r_idx = 0; r_idx < _radius_v.size(); ++r_idx) {
	const auto  radius = _radius_v[r_idx];
	const auto& temp_xs_v = temp_xs_vv[r_idx];
	
	_pt_xing_vv[plane][r_idx] = temp_xs_v.size();

	LLCV_DEBUG() << "@r_idx=" << r_idx << " radius=" << radius << " temp_xs_v sz=" << temp_xs_v.size() << std::endl;
	if (temp_xs_v.size()==2) {
	  LLCV_DEBUG() << "Testing pt1=" << temp_xs_v.front() << " pt2=" << temp_xs_v.back() << std::endl;
	  auto connected = larocv::Connected(mat,temp_xs_v.front(),temp_xs_v.back(),2);
	  LLCV_DEBUG() << "Connected=" << connected << std::endl;
	  _connected_vv[plane][r_idx] = (int)connected;
	  _distance_vv[plane][r_idx] = geo2d::dist(temp_xs_v.front(),temp_xs_v.back());
	}
      }
    }    
    LLCV_DEBUG() << "=======================" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;

    _out_tree->Fill();
    ResetTree();
    return 0.0;
  }
  
  void SelCosmicID::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;

    _fout->cd();
    _out_tree->Write();
    
    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelCosmicID::ResetTree() {
    _pt_xing_vv.clear();
    _connected_vv.clear();
    _distance_vv.clear();

    _pt_xing_vv.resize(3);
    _connected_vv.resize(3);
    _distance_vv.resize(3);

    for(size_t plane=0; plane<3; ++plane) {
      _pt_xing_vv[plane].resize(_radius_v.size(),-1);
      _connected_vv[plane].resize(_radius_v.size(),-1);
      _distance_vv[plane].resize(_radius_v.size(),-1);
    }

  }

}




#endif
