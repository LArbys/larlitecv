#ifndef __T3DPCMerge_H__
#define __T3DPCMerge_H__

/* ----------------------------------------------------------------
   This routine merges using first principal component of clusters

   ---------------------------------------------------------------*/

#include <vector>
#include "T3DCluster.h"

#include "GeoAlgo/GeoAlgo.h"

namespace larlitecv {

  class T3DPCMerge {
  public:

    T3DPCMerge();
    virtual ~T3DPCMerge() {};

    std::vector<T3DCluster> merge( const std::vector<T3DCluster>& tracks );


    // support functions
    std::vector<T3DCluster> endPointMerge( const std::vector<T3DCluster>& tracks );
    bool shouldWeEndPointMerge( const T3DCluster& ta, const T3DCluster& tb, double& closest_dist, std::vector<int>& whichends );
    void mergeTracks( T3DCluster& ta, T3DCluster& tb, const std::vector<int>& whichends );
    void setVerbosity( int v ) { m_verbose = v; };
    void configEndPointMerge( double max_endmerge_dist, double min_pcacos, int max_iterations ) {
      m_max_endmerge_dist = max_endmerge_dist;
      m_min_pcacos = min_pcacos;
      m_max_iterations = max_iterations;
    };
    
  protected:

    double m_max_endmerge_dist;
    double m_min_pcacos;
    int m_max_iterations;
    int m_verbose;
    ::larlite::geoalgo::GeoAlgo m_geoalgo;

  };


}

#endif
