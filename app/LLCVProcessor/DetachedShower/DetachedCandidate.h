#ifndef __DETACHEDCANDIDATE_H__
#define __DETACHEDCANDIDATE_H__

// ll
#include "DataFormat/cluster.h"

// locv
#include "LArOpenCV/ImageCluster/Base/ImageClusterTypes.h"
#include "LArOpenCV/ImageCluster/AlgoData/Vertex.h"

// cpp
#include <stdexcept>

namespace llcv {
  
  //
  // candidate cluster detached from vertex
  // defined by a bounding contour
  //
  class DetachedCluster {
  public:
    larocv::GEO2D_Contour_t ctor;
    float start_x;
    float start_y;
    int plane;
  };
  
  //
  // a grouping of upto 3 one-per-plane DetachedClusters
  //
  class DetachedCandidate {
  public:
    DetachedCandidate() { _dcluster_v.resize(3); }
    ~DetachedCandidate() {}

    void Insert(const llcv::DetachedCluster& dc,int plane);
    void Move(llcv::DetachedCluster&& dc,int plane);
    
    const std::vector<llcv::DetachedCluster>& CandidateClusters() const 
    { return _dcluster_v; }
    
    const llcv::DetachedCluster& CandidateCluster(int plane) const
    { if (plane >=3 or plane<0) throw std::runtime_error("Invalid plane"); return _dcluster_v[plane]; }

    larocv::data::Vertex3D origin;
    larocv::data::Vertex3D start;
    
  private:
    std::vector<llcv::DetachedCluster> _dcluster_v;

    
  };

}

#endif
