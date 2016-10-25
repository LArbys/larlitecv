#ifndef __BM_TRACK_CLUSTER_3D__
#define __BM_TRACK_CLUSTER_3D__

#include <vector>
#include "BoundaryEndPt.h"
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
    std::vector<float> start3D;
    std::vector<float> end3D;
    std::vector< BoundaryEndPt > start_endpts;
    std::vector< BoundaryEndPt > end_endpts;
    BoundaryEndPt::BoundaryEnd_t start_type;
    BoundaryEndPt::BoundaryEnd_t end_type;
    std::vector< BMTrackCluster2D > plane_paths;
    std::vector< std::vector<float> > path3d;
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
