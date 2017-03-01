#include "SMClusterTypes.h"

namespace larlitecv {

	void untagged_cluster_info_t::makeLink( ClusterIndex_t a, ClusterIndex_t b, 
		dbscan::ClusterExtrema::Extrema_t ex_a, dbscan::ClusterExtrema::Extrema_t ex_b, float dist ) {

		ClusterLink_t link(a,b,ex_a,ex_b,dist);
		if ( links.find(a)==links.end() ) {
			links.insert( std::make_pair(a,std::vector<ClusterLink_t>() ) );
		}
		if ( links.find(b)==links.end() ) {
			links.insert( std::make_pair(b,std::vector<ClusterLink_t>() ) );
		}
		links[a].push_back( link );
		links[b].push_back( link );
		link_v.push_back( link );
	}

	const std::vector<ClusterLink_t>& untagged_cluster_info_t::getLinks( ClusterIndex_t a ) {
		if ( links.find(a)==links.end() ) {
			links.insert( std::make_pair(a,std::vector<ClusterLink_t>() ) );
		}
		return links[a];
	}

}