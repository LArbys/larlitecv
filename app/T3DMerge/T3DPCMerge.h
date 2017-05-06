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

  protected:

    double m_max_endmerge_dist;
    geoalgo::GeoAlgo m_geoalgo;

  };


}

#endif
