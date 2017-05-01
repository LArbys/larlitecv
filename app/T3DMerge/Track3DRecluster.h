#ifndef __TRACK3D_RECLUSTER_H__
#define __TRACK3D_RECLUSTER_H__

#include <vector>

namespace larlitecv {

  class Track3DRecluster {
  public:
    Track3DRecluster() {};
    virtual ~Track3DRecluster() {};

    void addPath( const std::vector< std::vector<double> >& path );
    void addPath( const std::vector< std::vector<float> >& path );    

    std::vector< T3DCluster > m_tracks;
    
  };
    
  
}

#endif
