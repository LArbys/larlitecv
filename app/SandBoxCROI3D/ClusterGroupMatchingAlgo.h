#ifndef __CLUSTERGROUPMATCHINGALGO_H__
#define __CLUSTERGROUPMATCHINGALGO_H__

#include "ClusterGroupAlgo.h"
#include "ClusterGroupMatchingTypes.h"

namespace larlitecv {

	class ClusterGroupMatchingAlgo {
	public:
		ClusterGroupMatchingAlgo() {};
		virtual ~ClusterGroupMatchingAlgo() {};

		std::vector<ChargeVolume> MatchClusterGroups( const std::vector<larcv::Image2D>& untagged_v, const std::vector<PlaneClusterGroups>& plane_groups );

		void debugSetTargetCombo( const std::vector<int>& target ) { m_debug_targetcombo=target; };
		void debugUnsetTargetCombo() { m_debug_targetcombo.clear(); };

	protected:



		struct AlgoData_t {
			std::vector<PreMatchMetric_t> prematch_combos_v;
		};

		//
		std::vector<int> m_debug_targetcombo;

		// should use chain of command here
		void GenPreMatches( const std::vector<PlaneClusterGroups>& plane_groups, AlgoData_t& data );

		ChargeVolume GetIntersectionVolume( const std::vector<larcv::Image2D>& untagged_v, const PreMatchMetric_t& prematch );

		std::vector<int> GetCommonRowInterval( const PreMatchMetric_t& prematch );

		std::vector<WireInterval > GetWireInterval( const PreMatchMetric_t& prematch, const int row_start, const int row_end, const larcv::ImageMeta& meta );

		PointList_t GetIntersectionPoints( const std::vector<WireInterval>& plane_wire_intervals );

		Point_t GetIntersectionCentroid( const std::vector< WireInterval >& plane_wire_intervals );

		PointList_t GetBoundaryPoints( const PointList_t& crossingpts, const Point_t& centroid, const std::vector<WireInterval>& wireranges );

	  std::vector<WireInterval> RecalculateWireIntervalsFromBoundary( const PointList_t& yzboundary );

	  PointList_t EnforceTPCBounds( const PointList_t& yzboundary );

  	std::vector<float> SumContainedCharge( const PreMatchMetric_t& prematch, const std::vector<larcv::Image2D>& untagged_v, 
  		const std::vector<WireInterval>& overlap_intervals, const int row_start, const int row_end );

	};


}

#endif