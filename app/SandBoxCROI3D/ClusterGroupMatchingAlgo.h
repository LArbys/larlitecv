#ifndef __CLUSTERGROUPMATCHINGALGO_H__
#define __CLUSTERGROUPMATCHINGALGO_H__

#include "ClusterGroupAlgo.h"

namespace larlitecv {

	class ClusterGroupMatchingAlgo {
	public:
		ClusterGroupMatchingAlgo() {};
		virtual ~ClusterGroupMatchingAlgo() {};

		void MatchClusterGroups( const std::vector<PlaneClusterGroups>& plane_groups );

	protected:

		// triple for tracking pairs of indicices
		class GroupIndices_t : public std::vector<int> {
		public:
			GroupIndices_t(int idx1, int idx2, int idx3 ) {
				resize(3,0);
				at(0) = idx1;
				at(1) = idx2;
				at(2) = idx3;
			};
			virtual ~GroupIndices_t() {};
		};


		// We use this to metric to rank the matches.  Determines order of processing them.
		class PreMatchMetric_t {
		public:
			PreMatchMetric_t( const ClusterGroup& p1, const ClusterGroup& p2, const ClusterGroup& p3 ) { // assumes 3-plane match
		    dtSpan = 0;
				dtEnd = 0;

				m_clusters.push_back( &p1 );
				m_clusters.push_back( &p2 );
				m_clusters.push_back( &p3 );

    		for ( size_t p1=0; p1<m_clusters.size(); p1++ ) {
    			for ( size_t p2=p1+1; p2<m_clusters.size(); p2++ ) {
    				dtSpan += std::fabs( m_clusters.at(p1)->tick_width - m_clusters.at(p2)->tick_width );
    				dtEnd += std::fabs( m_clusters.at(p1)->tick_start - m_clusters.at(p2)->tick_start );
    				dtEnd += std::fabs( m_clusters.at(p1)->tick_end - m_clusters.at(p2)->tick_end );				
    			}
    		}				
			};
			virtual ~PreMatchMetric_t() {};

			bool operator()(const PreMatchMetric_t& lhs, const PreMatchMetric_t& rhs ) { //< sorting comparator
				if ( lhs.dtSpan<rhs.dtSpan)
					return true;
				else if ( lhs.dtSpan==rhs.dtSpan && lhs.dtEnd<rhs.dtEnd )
					return true;
				return false;
			};

			bool operator<( const PreMatchMetric_t& rhs ) { //< sorting comparator
				return (*this)( *this, rhs );
			};

			float dtSpan; // difference in time-spane
			float dtEnd;  // difference in end-points
			//float dQ;   // difference in charge (in future?)			
			std::vector< const ClusterGroup* > m_clusters;
		};

		struct AlgoData_t {
			std::vector<PreMatchMetric_t> prematch_combos_v;
		};

		// should use chain of command here
		void GenPreMatches( const std::vector<PlaneClusterGroups>& plane_groups, AlgoData_t& data );

	};


}

#endif