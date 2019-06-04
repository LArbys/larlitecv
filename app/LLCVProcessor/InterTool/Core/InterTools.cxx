#ifndef __INTERTOOLS_CXX__
#define __INTERTOOLS_CXX__

#include "InterTools.h"

namespace llcv {
  
  void MaskImage(const std::vector<larcv::PGraph>& pgraph_v,
  		 const std::map<larcv::PlaneID_t, std::vector<larcv::Pixel2DCluster> >& pix_m,
  		 const std::map<larcv::PlaneID_t, std::vector<larcv::ImageMeta> >&   pix_meta_m,
  		 std::vector<larcv::Image2D>& img_v) {
    
    for(const auto& pgraph : pgraph_v) {
      auto const& roi_v         = pgraph.ParticleArray();
      auto const& cluster_idx_v = pgraph.ClusterIndexArray();

      for(size_t roid=0; roid < roi_v.size(); ++roid) {
	
  	const auto& roi = roi_v[roid];
  	auto cidx = cluster_idx_v.at(roid);

  	for(size_t plane=0; plane<3; ++plane) {

  	  auto iter_pix = pix_m.find(plane);
  	  if(iter_pix == pix_m.end())
  	    continue;

  	  auto iter_pix_meta = pix_meta_m.find(plane);
  	  if(iter_pix_meta == pix_meta_m.end())
  	    continue;
	  
  	  auto const& pix_v      = (*iter_pix).second;
  	  auto const& pix_meta_v = (*iter_pix_meta).second;
	  
  	  auto const& pix      = pix_v.at(cidx);
  	  auto const& pix_meta = pix_meta_v.at(cidx);
	  
  	  auto& plane_img = img_v.at(plane);

  	  for(const auto& px : pix) {
  	    auto posx = pix_meta.pos_x(px.Y());
  	    auto posy = pix_meta.pos_y(px.X());
  	    auto row = plane_img.meta().row(posy);
  	    auto col = plane_img.meta().col(posx);
  	    plane_img.set_pixel(row,col,0);
  	  }
  	} // end plane
      } // end roi
    } // end pgraph
  }

}

#endif
