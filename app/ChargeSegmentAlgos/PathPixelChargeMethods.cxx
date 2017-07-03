#include "PathPixelChargeMethods.h"

namespace larlitecv {

  std::vector< PixelQPt > getPixelQPts( const std::vector< std::vector<double> >& path,
					const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					const int pixel_radius, const float threshold, const float max_stepsize ) {
    std::vector< PixelQPt > pathq;

    return pathq;
  }

  double getTrackTotalPixelCharge( const std::vector< std::vector<double> >& path,
				   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
				   const int pixel_radius, const float threshold, const float max_stepsize ) {
    return 0;
  }
  
}
