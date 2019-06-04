#ifndef MATCHUTILS_CXX
#define MATCHUTILS_CXX

#include "MatchUtils.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include <numeric>

namespace llcv {
  
  void
  Project3D(const larcv::ImageMeta& meta,
	    double parent_x,
	    double parent_y,
	    double parent_z,
	    double parent_t,
	    uint plane,
	    double& xpixel, double& ypixel) 
  {
    
    auto geohelp = larutil::GeometryHelper::GetME();
    auto larpro  = larutil::LArProperties::GetME();

    auto vtx_2d = geohelp->Point_3Dto2D(parent_x, parent_y, parent_z, plane );
    
    double x_compression  = (double)meta.width()  / (double)meta.cols();
    double y_compression  = (double)meta.height() / (double)meta.rows();
    xpixel = (vtx_2d.w/geohelp->WireToCm() - meta.tl().x) / x_compression;
    ypixel = (((parent_x/larpro->DriftVelocity() + parent_t/1000.)*2+3200)-meta.br().y)/y_compression;
  }
  
  void
  mask_image(larcv::Image2D& target, const larcv::Image2D& ref)
  {
    if(target.meta() != ref.meta()) 
      throw larcv::larbys("Cannot mask images w/ different meta");

    auto meta = target.meta();
    std::vector<float> data = target.move();
    auto const& ref_vec = ref.as_vector();

    for(size_t i=0; i<data.size(); ++i) { if(ref_vec[i]>0) data[i]=0; }	

    target.move(std::move(data));
  }
  
  float TestPixelType(int row, int col,
		      const larcv::Image2D& adc_img, const larcv::Image2D& pgraph_img,
		      bool ignore_zero) {

    float res = -2;

    int maxrow = pgraph_img.meta().rows() - 1;
    int maxcol = pgraph_img.meta().cols() - 1;

    static std::vector<int> pixel_type_v;
    static std::vector<int> row_v(3);
    static std::vector<int> col_v(3);

    row_v[0] = row;
    row_v[1] = row+1;
    row_v[2] = row-1;
    col_v[0] = col;
    col_v[1] = col+1;
    col_v[2] = col-1;

    for(auto& v : row_v) { 
      if(v<0)      v=0;
      if(v>maxrow) v=maxrow;
    }

    for(auto& v : col_v) { 
      if(v<0)      v=0;
      if(v>maxcol) v=maxcol;
    }
    
    pixel_type_v.clear();

    bool valid_px = false;
    
    for(size_t rid=0; rid<3; ++rid) {
      for(size_t cid=0; cid<3; ++cid) {
	
	if (adc_img.pixel(row_v[rid],col_v[cid]) == 0)
	  continue;

	int px = (int) pgraph_img.pixel(row_v[rid],col_v[cid]);

	valid_px = true;

	if (px < 0) continue;
	if (px == 0 and ignore_zero) continue;

	if (px >= (int) pixel_type_v.size())
	  pixel_type_v.resize(px+1,0);

	pixel_type_v[px] += 1;

      }
    }

    if (pixel_type_v.empty()) {
      if (!valid_px) return -2;
      else           return -1;
    }
    
    static std::vector<int>::iterator res_iter;

    res_iter = std::max_element(std::begin(pixel_type_v), std::end(pixel_type_v));
    res      = (float) std::distance(std::begin(pixel_type_v), res_iter);
    
    return res;
  }

}

#endif
