#ifndef __SMCLUSTERTYPES__
#define __SMCLUSTERTYPES__

// Here we define some types used by the StopMuCluster algorithm

#include <vector>
#include <map>

// LArCV
#include "dbscan/DBSCANAlgo.h"

namespace larlitecv {

	typedef int ClusterIndex_t;
	// link between clusters
	struct ClusterLink_t {
		std::vector<ClusterIndex_t> indices;
		std::vector<dbscan::ClusterExtrema::Extrema_t> extrema;
		float dist;
		ClusterLink_t( ClusterIndex_t idx_a, ClusterIndex_t idx_b, dbscan::ClusterExtrema::Extrema_t ex_a, dbscan::ClusterExtrema::Extrema_t ex_b, float d ) {
			indices.resize(2);
			indices[0] = idx_a;
			indices[1] = idx_b;
			extrema.resize(2);
			extrema[0] = ex_a;
			extrema[1] = ex_b;
			dist = d;
		};
	};

	// input pixels and clusters stored together
  struct untagged_cluster_info_t {
    ::dbscan::dbPoints pixels;
    ::dbscan::dbscanOutput output;
    std::vector<dbscan::ClusterExtrema> extrema_v;
    std::vector<ClusterLink_t> link_v;
    std::map< ClusterIndex_t, std::vector<ClusterLink_t> > links;
    void makeLink( ClusterIndex_t a, ClusterIndex_t b, dbscan::ClusterExtrema::Extrema_t ex_a, dbscan::ClusterExtrema::Extrema_t ex_b, float dist );
    const std::vector<ClusterLink_t>& getLinks( ClusterIndex_t a);
  };

}

#endif