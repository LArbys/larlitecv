#ifndef __SEARCHALGOBASE_CXX__
#define __SEARCHALGOBASE_CXX__

#include "SearchAlgoBase.h"

namespace llcv {

  const std::string& SearchAlgoBase::Name() const { return _name; }

  std::vector<llcv::DetachedCandidate>
  SearchAlgoBase::Search(larocv::data::Vertex3D& vtx3d,
			 const std::vector<std::tuple<cv::Mat,larocv::ImageMeta> >& adc_mat_meta_v,
			 const std::vector<std::tuple<cv::Mat,larocv::ImageMeta> >& shr_mat_meta_v) {
    
    std::vector<cv::Mat> adc_mat_v(3);
    std::vector<cv::Mat> shr_mat_v(3);
    std::vector<larocv::ImageMeta> adc_meta_v(3);

    vtx3d.vtx2d_v.resize(3);
    
    size_t nvalid = 0;
    for(size_t plane=0; plane<3; ++plane) {
      LLCV_DEBUG() << "project @plane=" << plane << std::endl;
      const auto& adc_mat_meta = adc_mat_meta_v[plane];
      const auto& shr_mat_meta = shr_mat_meta_v[plane];
      
      // copy the image
      adc_mat_v[plane]  = std::get<0>(adc_mat_meta).clone();
      adc_meta_v[plane] = std::get<1>(adc_mat_meta);

      shr_mat_v[plane]  = std::get<0>(shr_mat_meta).clone();

      // set the meta
      _geo.ResetPlaneInfo(adc_meta_v[plane]);

      // set the 2d point in the plane
      auto& vtx2d  = vtx3d.vtx2d_v[plane];

      float x = kINVALID_FLOAT;
      float y = kINVALID_FLOAT;
      
      try { x = _geo.x2col (vtx3d.x,plane); } 
      //catch(const larcv::larbys& what) { 
      catch(...) { 
      	LLCV_DEBUG() << "skip @ plane=" << plane << std::endl;
      	continue;
      }
      
      try { y = _geo.yz2row(vtx3d.y,vtx3d.z,plane); }
      //catch(const larcv::larbys& what) { 
      catch(...) { 
      	LLCV_DEBUG() << "skip @ plane=" << plane << std::endl;
      	continue; 
      }
      
      LLCV_DEBUG() << " @ vtx3d " << &vtx3d << " (x,y,z)=("<<vtx3d.x<<","<<vtx3d.y<<","<<vtx3d.z<<")"<<std::endl;
      LLCV_DEBUG() << " => (x,y)=(col,row)=("<<x<<","<<y<<")"<<std::endl;
	
      vtx2d.pt.x = x;
      vtx2d.pt.y = y;
      
      nvalid += 1;
    }

    // search
    LLCV_DEBUG() << "_Search_" << std::endl;

    if (nvalid<2) return std::vector<llcv::DetachedCandidate>();

    return _Search_(vtx3d,adc_mat_v,shr_mat_v,adc_meta_v);
  }
  
}
#endif
