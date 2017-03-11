#include "ClusterGroupMatchingAlgo.h"
#include "Combinator/Combinator.h"

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

		Combinator< ClusterGroup > combo( plane_groups );

		while ( !combo.isLast() ) {
			std::vector< const ClusterGroup* > groupcombo = combo.getCombo();
			PreMatchMetric_t prematch( *groupcombo[0], *groupcombo[1], *groupcombo[2] );
			data.prematch_combos_v.emplace_back( std::move(prematch) );
			combo.next();
		}

		std::sort( data.prematch_combos_v.begin(), data.prematch_combos_v.end() );
	}


}