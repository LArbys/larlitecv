#ifndef __PATH_PIXEL_CHARGE_TYPES_H__
#define __PATH_PIXEL_CHARGE_TYPES_H__

#include <vector>

// larcv
#include "DataFormat/Image2D.h"

namespace larlitecv {

  class PixelQPt {

    // unusable default constructor
    PixelQPt() {};

  public:
    
    PixelQPt( const std::vector<double>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
	      const int pixel_radius, const float threshold );
    virtual ~PixelQPt() {};

    const std::vector<double>& getPos() const { return m_pos3d; };
    const std::vector<double>& getPlaneCharge() const { return m_planeq; };
    int row() const { return m_imagerow; };
    const std::vector<int>& cols() const { return m_imagecols; };
    //const larcv::ImageMeta& meta() const ;
    const std::vector<int> getGoodPixelsPerPlane() const { return m_goodpixels_per_plane; };
    const std::vector<int> getBadPixelsPerPlane() const { return m_badchpixels_per_plane; };

  protected:

    std::vector<double> m_pos3d;
    std::vector<double> m_planeq;
    int m_imagerow;
    std::vector<int> m_imagecols;
    std::vector<int> m_goodpixels_per_plane;
    std::vector<int> m_badchpixels_per_plane;
    //larcv::ImageMeta meta;

  };
  
}

#endif
