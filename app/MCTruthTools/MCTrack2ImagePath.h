#ifndef __MCTrack2ImagePath_h__
#define __MCTrack2ImagePath_h__

#include <vector>

// larlite
#include "DataFormat/mctrack.h"

// larcv
#include "DataFormat/Image2D.h"
#include "SCE/SpaceChargeMicroBooNE.h"

namespace larlitecv {
  
  std::vector< std::vector<double> > mctrack2tyz( const larlite::mctrack& truthtrack,
						  const float trig_time,
						  bool returnxpos=false,
						  larlitecv::SpaceChargeMicroBooNE* psce=NULL );
  
  std::vector< std::vector<int> > tyzpath2imagepath( const std::vector< std::vector<double> >& path_tyz_sce, const std::vector<larcv::Image2D>& img_v );
  
  std::vector< std::vector<int> > mctrack2imagepath( const std::vector<larcv::Image2D>& img_v,
						     const larlite::mctrack& truthtrack,
						     const float trig_time,
						     larlitecv::SpaceChargeMicroBooNE* psce=NULL );
  
}

#endif
