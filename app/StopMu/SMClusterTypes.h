#ifndef __SMCLUSTERTYPES__
#define __SMCLUSTERTYPES__

// Here we define some types used by the StopMuCluster algorithm

// LArCV
#include "dbscan/DBSCANAlgo.h"

namespace larlitecv {

	// input pixels and clusters stored together
  struct untagged_cluster_info_t {
    ::dbscan::dbPoints pixels;
    ::dbscan::dbscanOutput output;
    std::vector<dbscan::ClusterExtrema> extrema_v;
  };

}

#endif