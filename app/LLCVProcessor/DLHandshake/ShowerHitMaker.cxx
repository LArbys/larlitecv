#ifndef SHOWERHITMAKER_CXX
#define SHOWERHITMAKER_CXX

#include "ShowerHitMaker.h"
#include "DataFormat/hit.h"
#include <cassert>

namespace llcv {

  void ShowerHitMaker::configure(const larcv::PSet&) {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
  }

  void ShowerHitMaker::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
  }

  bool ShowerHitMaker::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {
    LLCV_DEBUG() << "start" << std::endl;

    const auto ev_pixel2d  = (larcv::EventPixel2D*) mgr.get_data(larcv::kProductPixel2D,"dlshr");
    const auto ev_hit_in   = (larlite::event_hit*)  sto.get_data(larlite::data::kHit,"gaushit");
    auto ev_hit_out        = (larlite::event_hit*)  sto.get_data(larlite::data::kHit,"dlshr");
    
    std::vector<larcv::Image2D> shr_img_v;
    auto roi_exists = MakeShowerImage(ev_pixel2d,shr_img_v);
    if(!roi_exists) {
      LLCV_INFO() << "no ROI exists..." << std::endl;
      return true;
    }

    for(size_t plane=0; plane<3; ++plane) {
      
      const auto& shr_img = shr_img_v[plane];
      const auto& meta    = shr_img.meta();

      int rows = meta.rows();
      int cols = meta.cols();

      for (auto const& hit : (*ev_hit_in)) {
	    
	if (hit.WireID().Plane != plane) continue;
	    
	auto time = hit.PeakTime() + 2400;
	auto wire = hit.WireID().Wire;
      
	auto xpixel = (wire - meta.min_x()) / meta.pixel_width();
	auto ypixel = (meta.max_y() - time) / meta.pixel_height();
	    
	int xx = (int)(xpixel+0.5);
	int yy = (int)(ypixel+0.5);
	    
	if ( (xx < 0) || (yy < 0) ) continue;
	if ( (xx >= cols) || (yy >= rows) ) continue;
	
	auto px_value  = shr_img.pixel(yy,xx);

	if (px_value>0) ev_hit_out->push_back(hit);
	
      } // end this hit
    } // end this plane
    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  
  bool ShowerHitMaker::MakeShowerImage(const larcv::EventPixel2D* ev_pixel2d,
				       std::vector<larcv::Image2D>& shr_img_v) {
				       
    shr_img_v.clear();
    shr_img_v.resize(3);

    const auto& pixel_m = ev_pixel2d->Pixel2DArray();
    const auto& meta_m  = ev_pixel2d->MetaArray();

    if (pixel_m.empty()) return false;
    if (meta_m.empty()) return false;
    
    for(size_t plane=0; plane<3; ++plane) {

      auto meta_itr = meta_m.find(plane);
      if (meta_itr == meta_m.end()) return false;
      const auto& meta = (*meta_itr).second;

      auto& shr_image = shr_img_v[plane];
      shr_image = larcv::Image2D(meta);
      shr_image.paint(0.0);

      size_t npx_x = meta.cols();
      size_t npx_y = meta.rows();

      auto pixel_itr = pixel_m.find(plane);
      if (pixel_itr == pixel_m.end()) continue;
      const auto& pixel_v = (*pixel_itr).second;

      for(const auto& pixel : pixel_v) {
	if(pixel.X() >= npx_x || pixel.Y() >= npx_y) {
	  LLCV_WARNING() << "Ignoring shower pixel (row,col) = ("
			  << pixel.Y() << "," << pixel.X() << ")"
			 << " as it is out of bounds (ncol=" << npx_x << ", nrow=" << npx_y << ")" << std::endl;
	  continue;
	}
	shr_image.set_pixel(pixel.Y(),pixel.X(),pixel.Intensity());
      } // end pixel
    } // end plane
	
    return true;
  }

  void ShowerHitMaker::finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
  }
}

#endif
	
