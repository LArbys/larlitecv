#ifndef __BM_TRACK_CLUSTER2D__
#define __BM_TRACK_CLUSTER2D__

/** 
    Needed a data type to store 2D track clusters.
    
    Will make the process of modifying and merging them easier.
*/

#include "DataFormat/Pixel2DCluster.h"
#include "BoundaryEndPt.h"

namespace larlitecv {

  class BMTrackCluster2D {
    
  public:
    BMTrackCluster2D() {
      start_pix_idx = -1;
      end_pix_idx = -1;
    };
    virtual ~BMTrackCluster2D() {};

    BoundaryEndPt start;
    BoundaryEndPt end;
    larcv::Pixel2DCluster pixelpath;
    int plane;
    int start_pix_idx;
    int end_pix_idx;

  };

}

#endif