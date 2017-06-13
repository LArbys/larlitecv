#ifndef __PATH_2_PIXELS_H__
#define __PATH_2_PIXELS_H__

#include <vector>

#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2DCluster.h"

namespace larlitecv {

  // w/ badch
  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImages( const std::vector< std::vector<double> >& path3d,
							       const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
							       const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
							       const float stepsize );
  // w/o badch
  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImagesNoBadCh( const std::vector< std::vector<double> >& path3d,
								      const std::vector<larcv::Image2D>& imgs, 
								      const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
								      const float stepsize );

}

#endif
