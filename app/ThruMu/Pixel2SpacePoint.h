#ifndef __PIXEL_2_SPACE_POINT_H__
#define __PIXEL_2_SPACE_POINT_H__

#include <vector>
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2D.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundarySpacePoint.h"

namespace larlitecv {

  BoundarySpacePoint Pixel2SpacePoint( const std::vector<larcv::Pixel2D>& pixels, const BoundaryEnd_t endtype, const larcv::ImageMeta& meta );

}

#endif
