#include "ClusterGroupMatchingAlgo.h"

namespace larlitecv {

	void ClusterGroupMatchingAlgo::MatchClusterGroups( const std::vector<PlaneClusterGroups>& plane_groups ) {

		// setup internal data object
		AlgoData_t data;

		GenPreMatches( plane_groups, data );
		for (int i=0; i<data.prematch_combos_v.size(); i++) {
			PreMatchMetric_t& prematch = data.prematch_combos_v.at(i);
			std::cout << "prematch[" << i << "]: " << prematch.dtSpan << " " << prematch.dtEnd << std::endl;
		}

		// over minimal time overlap, break into time slices and calculate bounding polygons

	}

	void ClusterGroupMatchingAlgo::GenPreMatches( const std::vector<PlaneClusterGroups>& plane_groups, ClusterGroupMatchingAlgo::AlgoData_t& data ) {
		// make pre-matches
		// factorial fun!

		struct IndexCombo_t {
			std::vector<int> combo;
			const std::vector<int> nelems;
			IndexCombo_t( const std::vector<int>& size_per_dim ) 
			  : nelems(size_per_dim) {
				combo.resize( nelems.size(), 0 );
			};
			bool isLast() {
				for (size_t p=0; p<nelems.size(); p++) {
					if ( combo[p]!=nelems[p]-1 ) 
						return false;
				}
				return true;
			};
			bool next() {
				if ( isLast() ) return true; // prevent moving further
				for ( size_t p=0; p<nelems.size(); p++) {
					combo[p]++;
					if ( combo[p]<nelems[p] )
						break;
					else
						combo[p] = 0;
					// reset this value, move to next dim
				}
				return isLast();
			};
			int operator()(int i) {
				return combo[i];
			};
		};

		std::vector<int> ngroups_per_plane( plane_groups.size() );
		for ( size_t p=0; p<plane_groups.size(); p++)
			ngroups_per_plane[p] = plane_groups.at(p).size();
		IndexCombo_t combo( ngroups_per_plane );

		while ( !combo.isLast() ) {
			PreMatchMetric_t prematch( plane_groups.at(0).at(combo(0)), plane_groups.at(1).at(combo(1)), plane_groups.at(1).at(combo(2)) );
			data.prematch_combos_v.emplace_back( std::move(prematch) );
			combo.next();
		}

		std::sort( data.prematch_combos_v.begin(), data.prematch_combos_v.end() );
	}


}