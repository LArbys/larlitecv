#ifndef __BM_TRACK_CLUSTER_3D__
#define __BM_TRACK_CLUSTER_3D__

#include <vector>
#include "BoundaryEndPt.h"

namespace larlitecv {

  class BMTrackCluster3D {

  public:
    BMTrackCluster3D();
    virtual ~BMTrackCluster3D();

    std::vector<int> trackidx; // index of tracks in list of TrackClusters2D
    int row_start;
    int row_end;
    float tick_start;
    float tick_end;
    std::vector<int> start_wire;
    std::vector<int> end_wire;
    std::vector<float> start3D;
    std::vector<float> end3D;
    BoundaryEndPt::BoundaryEnd_t start_type;
    BoundaryEndPt::BoundaryEnd_t end_type;
    float score;

    bool operator< (const BMTrackCluster3D& rhs ) const {
      if ( trackidx[0]<rhs.trackidx[0] ) return true;
      else if ( trackidx[0]==rhs.trackidx[0] ) {
	if ( trackidx[1]<rhs.trackidx[1] ) 
	  return true;
	else if ( trackidx[1]==rhs.trackidx[1] && trackidx[2]<rhs.trackidx[2] ) 
	  return true;
      }
      return false;
    };
    
    bool operator== (const BMTrackCluster3D& rhs ) const  {
      if ( trackidx[0]==rhs.trackidx[0] 
	   && trackidx[1]==rhs.trackidx[1] 
	   && trackidx[2]==rhs.trackidx[2]  )
	return true;
      return false;
    };

  };


}

#endif
