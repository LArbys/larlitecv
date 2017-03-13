#ifndef __CLUSTERGROUPMATCHINGALGO_H__
#define __CLUSTERGROUPMATCHINGALGO_H__

#include "ClusterGroupAlgo.h"

namespace larlitecv {

	class ClusterGroupMatchingAlgo {
	public:
		ClusterGroupMatchingAlgo() {};
		virtual ~ClusterGroupMatchingAlgo() {};

		void MatchClusterGroups( const std::vector<larcv::Image2D>& untagged_v, const std::vector<PlaneClusterGroups>& plane_groups );

		void debugSetTargetCombo( const std::vector<int>& target ) { m_debug_targetcombo=target; };
		void debugUnsetTargetCombo() { m_debug_targetcombo.clear(); };

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
			PreMatchMetric_t( const ClusterGroup& p1, const ClusterGroup& p2, const ClusterGroup& p3 );
			virtual ~PreMatchMetric_t() {};

			bool operator()(const PreMatchMetric_t& lhs, const PreMatchMetric_t& rhs ) { //< sorting comparator
				if ( lhs.dtSpan + lhs.dtEnd < rhs.dtSpan+rhs.dtEnd ) 
					return true;
				return false;
				// if ( lhs.dtSpan<rhs.dtSpan)
				// 	return true;
				// else if ( lhs.dtSpan==rhs.dtSpan && lhs.dtEnd<rhs.dtEnd )
				// 	return true;
				// return false;
			};

			bool operator<( const PreMatchMetric_t& rhs ) { //< sorting comparator
				return (*this)( *this, rhs );
			};

			float dtSpan; // difference in time-spane
			float dtEnd;  // difference in end-points
			//float dQ;   // difference in charge (in future?)			
			std::vector< const ClusterGroup* > m_clusters;
			std::vector< int > m_index_combo;
		};

		// Intersection Points
		class Point_t : public std::vector<float>  {
		public:
			Point_t( float y, float z, float start, float end )
			 : tickstart(start), tickend(end) { 
			  resize(2); 
			  at(0)=y; 
			  at(1)=z; 
			};
			virtual ~Point_t() {};
			float tickstart;
			float tickend;
			const std::vector<float> as_vec() const {
				std::vector<float> vec(2);
				for (int i=0; i<2; i++) 
					vec[i] = (*this)[i];
				return vec;
			};
		};

		class WireInterval : public std::vector<int> {
		public:
			WireInterval() {
				resize(2);
				at(0) = -1;
				at(1) = -1;
			};
			WireInterval(int w1, int w2) {
				resize(2);
				at(0) = w1;
				at(1) = w2;
			};
			virtual ~WireInterval() {};
		};

		typedef std::vector<Point_t> PointList_t;
		typedef std::vector<PointList_t> Slices_t;

		class ChargeVolumes {

		};


		struct AlgoData_t {
			std::vector<PreMatchMetric_t> prematch_combos_v;
		};

		//
		std::vector<int> m_debug_targetcombo;

		// should use chain of command here
		void GenPreMatches( const std::vector<PlaneClusterGroups>& plane_groups, AlgoData_t& data );

		Slices_t GetIntersectionVolume( const std::vector<larcv::Image2D>& untagged_v, const PreMatchMetric_t& prematch );

		std::vector<int> GetCommonRowInterval( const PreMatchMetric_t& prematch );

		std::vector<WireInterval > GetWireInterval( const PreMatchMetric_t& prematch, const int row_start, const int row_end, const larcv::ImageMeta& meta );

		PointList_t GetIntersectionPoints( const std::vector<WireInterval>& plane_wire_intervals );

		Point_t GetIntersectionCentroid( const std::vector< WireInterval >& plane_wire_intervals );

		PointList_t GetBoundaryPoints( const PointList_t& crossingpts, const Point_t& centroid, const std::vector<WireInterval>& wireranges );

	  std::vector<WireInterval> RecalculateWireIntervalsFromBoundary( const PointList_t& yzboundary );

	};


}

#endif