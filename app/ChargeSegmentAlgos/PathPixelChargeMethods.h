#ifndef __PATH_PIXEL_CHARGE_METHODS__
#define __PATH_PIXEL_CHARGE_METHODS__

#include <vector>
#include "PathPixelChargeTypes.h"
#include "DataFormat/Image2D.h"

namespace larlitecv {

  // This function returns a list of 3D points and the pixel charge around those points
  // The charge is simply summed in the neighborhood of the projected point: pixels are not unique
  // Do not sum the charge together, you will over count
  std::vector< PixelQPt > getPixelQPts( const std::vector< std::vector<double> >& path,
					const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					const int pixel_radius, const float threshold, const float max_stepsize );

  double getTrackTotalPixelCharge( const std::vector< std::vector<double> >& path,
				   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
				   const int pixel_radius, const float threshold, const float max_stepsize );
  
}

#endif
