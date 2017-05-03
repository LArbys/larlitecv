#ifndef __TRACK3D_RECLUSTER_H__
#define __TRACK3D_RECLUSTER_H__

#include <vector>
#include <set>

#include "TRandom3.h"

#include "GeoAlgo/GeoAlgo.h"

#include "T3DCluster.h"

namespace larlitecv {

  class Track3DRecluster {
  public:
    Track3DRecluster();
    virtual ~Track3DRecluster() {};

    void addPath( const std::vector< std::vector<double> >& path );
    void addPath( const std::vector< std::vector<float> >& path );
    std::vector< T3DCluster > recluster();
    bool mergeEnds( T3DCluster& fa, T3DCluster& fb, float max_dist, float min_cos );
    void setVerbosity( int v ) { m_verbosity=v; };

  protected:

    struct SegmentOverlap_t {
      std::vector<int> indices;
      std::set<int>    othertrackmatches;
    };
      
    
    std::vector< T3DCluster > m_tracks;

    bool ReclusterPair( const T3DCluster& tracka, const T3DCluster& trackb, std::vector<T3DCluster>& tracks_v );
    std::vector< SegmentOverlap_t > getOverlapSegmentsOfBonA( const T3DCluster& tracka, const T3DCluster& trackb, std::vector<int>& overlap_info );    

    geoalgo::GeoAlgo m_geoalgo;

    TRandom3 m_rand;
    int m_verbosity;
    int m_min_segoverlap;
    double m_ann_dist;
    double m_max_end_dist;
    double m_min_end_cos;    
  };
    
  
}

#endif
