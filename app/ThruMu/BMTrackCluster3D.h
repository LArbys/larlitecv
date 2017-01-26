#ifndef __BM_TRACK_CLUSTER_3D__
#define __BM_TRACK_CLUSTER_3D__

#include <vector>
#include "BoundaryMuonTaggerTypes.h"
#include "BoundaryEndPt.h"
#include "BoundarySpacePoint.h"
#include "BMTrackCluster2D.h"

namespace larlitecv {

  class BMTrackCluster3D {

  public:
    BMTrackCluster3D();
    virtual ~BMTrackCluster3D();

    int row_start;
    int row_end;
    float tick_start;
    float tick_end;
    std::vector<int> start_wire;
    std::vector<int> end_wire;
    std::vector<double> start3D;
    std::vector<double> end3D;
    BoundarySpacePoint start_endpts;
    BoundarySpacePoint end_endpts;
    BoundaryEnd_t start_type;
    BoundaryEnd_t end_type;
    std::vector< BMTrackCluster2D > plane_paths;
    std::vector< std::vector<double> > path3d;
    int track2d_index;
    
    //std::vector< float > node_triangle_area;

    bool operator< (const BMTrackCluster3D& rhs ) const {
      if ( tick_start<rhs.tick_start ) return true;
      if ( tick_start==rhs.tick_start && tick_end<rhs.tick_end ) return true;
      return false;
    };
    
    bool operator== (const BMTrackCluster3D& rhs ) const  {
      if ( tick_start==rhs.tick_start && tick_end==rhs.tick_end )
	return true;
      return false;
    };
    
  };


}

#endif
