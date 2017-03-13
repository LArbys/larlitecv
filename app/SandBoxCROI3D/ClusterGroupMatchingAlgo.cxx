#include "ClusterGroupMatchingAlgo.h"
#include "Combinator/Combinator.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {

	void ClusterGroupMatchingAlgo::MatchClusterGroups( const std::vector<larcv::Image2D>& untagged_v, const std::vector<PlaneClusterGroups>& plane_groups ) {

		// setup internal data object
		AlgoData_t data;

		std::vector<int> debug_target = {3,7,5};
		debugSetTargetCombo( debug_target );

		GenPreMatches( plane_groups, data );
		// for (size_t i=0; i<data.prematch_combos_v.size(); i++) {
			// PreMatchMetric_t& prematch = data.prematch_combos_v.at(i);
			// std::cout << "prematch[" << i << "]: " 
			  // << "(" << prematch.m_index_combo[0] << "," << prematch.m_index_combo[1] << "," << prematch.m_index_combo[2] << ") "
			  // << prematch.dtSpan << " " << prematch.dtEnd << std::endl;
		// }

		for (auto const& prematch : data.prematch_combos_v ) {
			//std::cout << "(" << m_debug_targetcombo[0] << "," << m_debug_targetcombo[1] << "," << m_debug_targetcombo[2] << ")" << std::endl;
			if ( m_debug_targetcombo.size()>0 && m_debug_targetcombo!=prematch.m_index_combo )
				continue;

			// over minimal time overlap, break into time slices and calculate bounding polygons
			Slices_t slices = GetIntersectionVolume( untagged_v, prematch );
			for ( auto const& slice : slices ) {
				std::cout << "slice: ";
				for ( auto const& pt : slice ) {
					std::cout << " (" << pt[0] << "," << pt[1] << ") ";
				}
				std::cout << std::endl;
			}
		}

	}

	void ClusterGroupMatchingAlgo::GenPreMatches( const std::vector<PlaneClusterGroups>& plane_groups, ClusterGroupMatchingAlgo::AlgoData_t& data ) {
		// make pre-matches
		// factorial fun!

		Combinator< ClusterGroup > combo( plane_groups );

		while ( !combo.isLast() ) {
			std::vector< const ClusterGroup* > groupcombo = combo.getCombo();
			PreMatchMetric_t prematch( *groupcombo[0], *groupcombo[1], *groupcombo[2] );
			prematch.m_index_combo = combo.getIndexCombo();
			data.prematch_combos_v.emplace_back( std::move(prematch) );
			combo.next();
		}

		std::sort( data.prematch_combos_v.begin(), data.prematch_combos_v.end() );
	}

	ClusterGroupMatchingAlgo::PreMatchMetric_t::PreMatchMetric_t( const ClusterGroup& p1, const ClusterGroup& p2, const ClusterGroup& p3 ) { // assumes 3-plane match

	  dtSpan = 0;
		dtEnd = 0;		

		m_clusters.push_back( &p1 );		
		m_clusters.push_back( &p2 );
		m_clusters.push_back( &p3 );

  	for ( size_t p1=0; p1<m_clusters.size(); p1++ ) {

  		for ( size_t p2=p1+1; p2<m_clusters.size(); p2++ ) {

  			dtSpan += std::fabs( m_clusters.at(p1)->tick_width - m_clusters.at(p2)->tick_width );
  			float diff1 = std::fabs( m_clusters.at(p1)->tick_start - m_clusters.at(p2)->tick_start ); // aligned
  			float diff2 = std::fabs( m_clusters.at(p1)->tick_start - m_clusters.at(p2)->tick_end );   // anti-aligned
  			if ( diff1<diff2 ) {
  				// aligned
  				dtEnd += diff1;
  				dtEnd += std::fabs( m_clusters.at(p1)->tick_end - m_clusters.at(p2)->tick_end );
  			}
  			else {
  				// anti-aligned
  				dtEnd += diff2;
  				dtEnd += std::fabs( m_clusters.at(p1)->tick_end - m_clusters.at(p2)->tick_start );
  			}
  		}
  	}				
	}

	ClusterGroupMatchingAlgo::Slices_t ClusterGroupMatchingAlgo::GetIntersectionVolume( const std::vector<larcv::Image2D>& untagged_v, 
		const ClusterGroupMatchingAlgo::PreMatchMetric_t& prematch ) {

		std::cout << "ClusterGroupMatchingAlgo::GetIntersectionVolume" << std::endl;

		const larcv::ImageMeta& meta = untagged_v.front().meta();

	  std::vector<int> row_interval = GetCommonRowInterval( prematch );
	  std::cout << " row interval: [" << row_interval[0] << "," << row_interval[1] << "] "
   	  << " ticks: [" << meta.pos_y(row_interval[1]) << "," << meta.pos_y(row_interval[0]) << "]"
   	  << std::endl;

   	int rows = abs( row_interval[1]-row_interval[0] );

   	int slice_size = 4;
   	int ntime_slices = rows/slice_size;
   	if ( rows%slice_size!=0 ) ntime_slices++;

   	Slices_t area_slices(ntime_slices);

   	for (int islice=0; islice<ntime_slices; islice++) {

   		// define the row slice
   		int row_start = row_interval[0] + slice_size*(islice);
   		int row_end   = row_start + slice_size;

   		// define the wire intervals over this slice
   		std::vector<WireInterval> wintervals = GetWireInterval( prematch, row_start, row_end, meta );

   		std::cout << " slice[" << islice << "] "
   			<< "tick=[" << meta.pos_y( row_end ) << "," << meta.pos_y(row_start) << "]"
   			<< " wire intervals: "
   			<< "[" << wintervals[0][0] << "," << wintervals[0][1] << "] "
   			<< "[" << wintervals[1][0] << "," << wintervals[1][1] << "] "
   			<< "[" << wintervals[2][0] << "," << wintervals[2][1] << "] ";

   		// get the (y,z) point from the center of each wire interval
   		Point_t centroid = GetIntersectionCentroid( wintervals );
   		std::cout << "centroid: (" << centroid[0] << "," << centroid[1] << ") ";
   		if ( centroid[0]<=-10000.0 && centroid[1]<=-10000.0 ) {
   			std::cout << " [bad slice]" << std::endl;
   			continue;
   		}

   		// get a list of all 12 intersection points made up by the intervals
   		PointList_t yzintersections = GetIntersectionPoints( wintervals );
   		std::cout << "nxsections=" << yzintersections.size() << " ";

   		// we filter out the intersection points along the (y,z) overlap polygon and form the boundary of the polygon
   		PointList_t yzboundary = GetBoundaryPoints( yzintersections, centroid, wintervals );
   		std::cout << "boundarypts=" << yzboundary.size() << " ";

   		// we calculate the area of this polygon

   		// we can use the polygon to narrow the wire-ranges
   		std::vector<WireInterval> overlap_intervals = RecalculateWireIntervalsFromBoundary( yzboundary );
   		std::cout << "overlap intervals: ";
   		for ( auto const& interval : overlap_intervals ) {
   			std::cout << "[" << interval[0] << "," << interval[1] << "] ";
   		}

   		// store the slice
   		area_slices.at(islice) = yzboundary;

   		std::cout << std::endl;
   	}

   	// volume is the set of slices
   	return area_slices;
	}

	std::vector<int> ClusterGroupMatchingAlgo::GetCommonRowInterval( const ClusterGroupMatchingAlgo::PreMatchMetric_t& prematch ) {
		// we get the minimal row interval common
		float row_start = -1.0;
		float row_end   = -1.0;

		// we simply scan the groups, getting the maximum start and minimum start
		for ( const auto& pgroup : prematch.m_clusters ) {
			if ( row_start<0 || pgroup->tick_start>row_start ) row_start = pgroup->tick_start;
			if ( row_end<0 || pgroup->tick_end<row_end ) row_end = pgroup->tick_end;
		}
		std::vector<int> interval = { (int)row_start, (int)row_end };
		return interval;
	}

	std::vector< ClusterGroupMatchingAlgo::WireInterval > ClusterGroupMatchingAlgo::GetWireInterval( 
		const ClusterGroupMatchingAlgo::PreMatchMetric_t& prematch, 
		const int row_start, const int row_end, const larcv::ImageMeta& meta ) {
		// we find the min and max col for a given row interval for each plane cluster
		std::vector< WireInterval > intervals;
		for ( size_t p=0; p<prematch.m_clusters.size(); p++) {
			int minwire = -1;
			int maxwire = -1;
		  const ClusterGroup& group = *(prematch.m_clusters.at(p));
			// for each cluster group, we cycle through the clusters and pick the pixels in range.
			// we can use the cluster extrema to rule some out
			for ( size_t ic=0; ic<group.m_clusters_v.size(); ic++ ) {

				const dbscan::ClusterExtrema& ex    = group.m_extrema_v.at(ic);
				const larcv::Pixel2DCluster& pixels = group.m_clusters_v.at(ic);

				bool start_is_inside = ex.bottommost()[1]<=row_start && row_start<=ex.topmost()[1];
				bool end_is_inside   = ex.bottommost()[1]<=row_end   && row_end<=ex.topmost()[1];
				bool stradles        = ex.bottommost()[1]>=row_start && row_end>=ex.topmost()[1];
				bool contained       = ex.bottommost()[1]<=row_start && row_end<=ex.topmost()[1];
				if ( start_is_inside || end_is_inside || stradles || contained ) {
					for ( auto const& pix : pixels ) {
						// is pixel within row interval
						if ( row_start<=(int)pix.Y() && (int)pix.Y()<=row_end ) {
							int wire = meta.pos_x( pix.X() );
							if ( minwire<0 || minwire>wire ) minwire = wire;
							if ( maxwire<0 || maxwire<wire ) maxwire = wire;
						}
					}
				}

			}
			WireInterval range(minwire,maxwire);
			intervals.emplace_back( std::move(range) );
		}
		return intervals;
	}

	ClusterGroupMatchingAlgo::PointList_t ClusterGroupMatchingAlgo::GetIntersectionPoints( 
		const std::vector< ClusterGroupMatchingAlgo::WireInterval >& plane_wire_intervals ) {

		PointList_t intersections_v;
		for (size_t p1=0; p1<plane_wire_intervals.size(); p1++) {
			for (size_t p2=p1+1; p2<plane_wire_intervals.size(); p2++) {
				const WireInterval& interval1 = plane_wire_intervals.at(p1);
				const WireInterval& interval2 = plane_wire_intervals.at(p2);
				for (int i=0; i<2; i++) {
					for (int j=0; j<2; j++) {

						// the interval value not found, so we skip it
						if ( interval1[i]<0 || interval2[j]<0 )
							continue;

						std::vector<float> intersection;
						int crosses = 0;
						larcv::UBWireTool::wireIntersection( p1, interval1[i], p2, interval2[j], intersection, crosses );
						Point_t pt( intersection[1], intersection[0], 0, 1 );
						intersections_v.emplace_back( std::move(pt) );
					}
				}
			}
		}

		return intersections_v;
	}

	ClusterGroupMatchingAlgo::Point_t ClusterGroupMatchingAlgo::GetIntersectionCentroid( const std::vector< ClusterGroupMatchingAlgo::WireInterval >& plane_wire_intervals ) {
		std::vector<int> wids;
		std::vector<int> goodplanes;

		for ( size_t p=0; p<plane_wire_intervals.size(); p++ ) {
			auto const& wints = plane_wire_intervals.at(p);
			int avewire = 0;
			int nwires = 0;
			for (int i=0; i<2; i++) {
				if ( wints[i]>=0 ){
					avewire += wints[i];
					nwires++;
				}
			}
			if (nwires>0 ) {
			  avewire /= nwires;
				wids.push_back( avewire );
				goodplanes.push_back(p);
			}
		}
		if ( goodplanes.size()==3 ) {
			// 3 wire intersection
	  	std::vector<float> inter;
	  	int crosses;
	  	double tri;
	  	larcv::UBWireTool::wireIntersection( wids, inter, tri, crosses );
	  	Point_t centroid( inter[1], inter[0], 0, 0 );
	  	return centroid;
	  }
	  else if ( goodplanes.size()==2 ) {
	  	// 2 wire intersection
	  	std::vector<float> inter;
	  	int crosses;
	  	larcv::UBWireTool::wireIntersection(  goodplanes[0], wids[0], goodplanes[1], wids[1], inter, crosses );
	  	Point_t centroid( inter[1], inter[0], 0, 0 );
	  	return centroid;
	  }

	  return Point_t(-1.0e6,-1.0e6,-1,-1);
	}

	ClusterGroupMatchingAlgo::PointList_t ClusterGroupMatchingAlgo::GetBoundaryPoints( const ClusterGroupMatchingAlgo::PointList_t& crossingpts, 
		const ClusterGroupMatchingAlgo::Point_t& centroid, const std::vector<ClusterGroupMatchingAlgo::WireInterval>& wireranges  ) {

		// we define an (y,z) point that we can sort by phi around centroid
		class BoundaryPt {
		public:
			BoundaryPt(const std::vector<float>& yz, const std::vector<float>& centroid ) 
			: m_yz(yz) {
				float r[2] = {0.0, 0.0};
				for (int i=0; i<2; i++)
					r[i] = yz[i] - centroid[i];
				m_phi = atan2( r[1], r[0] );
			};
			std::vector<float> m_yz;
			float m_phi;
			bool operator<( const BoundaryPt& rhs ) const {
				if ( m_phi<rhs.m_phi)
					return true;
				return false;
			};
		};

		float dyz[4][2] = { {0.6,0.0}, {0.0,0.6}, {-0.6,0.0}, {0.0,-0.6} };

		std::vector<BoundaryPt> boundary_v;
		for ( auto const& pt : crossingpts ) {
			// for each intersection, we bump the wire and find out if it is inside all wire intervals
			// it only has to be true for one of the four test points
			bool one_point_within_all_intervals = false;
			for (int i=0; i<4; i++) {
				Double_t test_xyz[3];
				test_xyz[0] = 0.;
				test_xyz[1] = pt[0]+dyz[i][0];
				test_xyz[2] = pt[1]+dyz[i][1];
				// must be in all intervals. so we search and stop if not in one of the plnes
				bool inrange = true;
				for (size_t p=0; p<wireranges.size(); p++) {
					float testwire = larutil::Geometry::GetME()->WireCoordinate( test_xyz, p );
					//std::cout << "  testpoint[" << i << "] plane=" << p << " (" << test_xyz[1] << "," << test_xyz[2] << ") "
					//	<< "wire=" << testwire << " range [" << wireranges[p][0] << "," << wireranges[p][1] << "]" << std::endl;
					if ( testwire<(float)wireranges[p][0] || testwire>(float)wireranges[p][1] ) {
						inrange = false;
						break;
					}
				}
				if ( inrange ) {
					// satisifed for all planes. we mark as passed and stop looking.
					one_point_within_all_intervals = true;
					break;
				}
			}

			if ( one_point_within_all_intervals ) {
				// this is a boundary point
				BoundaryPt bpt( pt.as_vec(), centroid.as_vec() );
				boundary_v.emplace_back( std::move(bpt) );
			}
		}//end of point loop

		// sort this
		std::sort( boundary_v.begin(), boundary_v.end() );

	  // make a point list
	  PointList_t boundary_pts;
	  for ( auto& bpt : boundary_v ) {
	  	Point_t pt( bpt.m_yz[0], bpt.m_yz[1], crossingpts.front().tickstart, crossingpts.front().tickend );
	  	boundary_pts.emplace_back( std::move(pt) );
	  }

		return boundary_pts;
	}

	std::vector<ClusterGroupMatchingAlgo::WireInterval> ClusterGroupMatchingAlgo::RecalculateWireIntervalsFromBoundary( const ClusterGroupMatchingAlgo::PointList_t& yzboundary ) {
		std::vector<WireInterval> output(3);

		for ( auto const& pt : yzboundary ) {
			Double_t xyz[3];
			xyz[0] = 0.;
			xyz[1] = pt[0];
			xyz[2] = pt[1];
			for ( size_t p=0; p<3; p++) {
				int wire = (int)larutil::Geometry::GetME()->WireCoordinate( xyz, p );
				if ( output[p][0]<0 || wire<output[p][0] )
					output[p][0] = wire;
				if ( output[p][1]<0 || wire>output[p][1] )
					output[p][1] = wire;
			}
		}
		return output;
	}


}