#include "ClusterGroupMatchingAlgo.h"
#include "Combinator/Combinator.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "TRandom.h"

namespace larlitecv {

	std::vector<ChargeVolume> ClusterGroupMatchingAlgo::MatchClusterGroups( const std::vector<larcv::Image2D>& untagged_v, 
		const std::vector<PlaneClusterGroups>& plane_groups ) {

		// setup internal data object
		AlgoData_t data;

		//std::vector<int> debug_target = {3,7,5}; // neutrino cluster
		//std::vector<int> debug_target = {3,0,3};
		//debugSetTargetCombo( debug_target );

		GenPreMatches( plane_groups, data );
		// for (size_t i=0; i<data.prematch_combos_v.size(); i++) {
			// PreMatchMetric_t& prematch = data.prematch_combos_v.at(i);
			// std::cout << "prematch[" << i << "]: " 
			  // << "(" << prematch.m_index_combo[0] << "," << prematch.m_index_combo[1] << "," << prematch.m_index_combo[2] << ") "
			  // << prematch.dtSpan << " " << prematch.dtEnd << std::endl;
		// }

		std::vector<ChargeVolume> vols;
		for (auto const& prematch : data.prematch_combos_v ) {
			//std::cout << "(" << m_debug_targetcombo[0] << "," << m_debug_targetcombo[1] << "," << m_debug_targetcombo[2] << ")" << std::endl;
			if ( m_debug_targetcombo.size()>0 && m_debug_targetcombo!=prematch.m_index_combo )
				continue;

			// over minimal time overlap, break into time slices and calculate bounding polygons
			//std::cout << "prematch index: (" << prematch.m_index_combo[0] << "," << prematch.m_index_combo[1] << "," << prematch.m_index_combo[2] << ")";
			ChargeVolume vol = GetIntersectionVolume( untagged_v, prematch );
		  vols.emplace_back( std::move(vol) );
		}

		std::sort( vols.begin(), vols.end() );

		// std::cout << "Charge Volumes: " << std::endl;
		// for ( auto const& vol : vols ) {
		// 	std::cout << " clgroup[" << vol._clustergroup_indices[0] << "," << vol._clustergroup_indices[1] << "," << vol._clustergroup_indices[2] << "] "
		// 	  << " numslices=" << vol.num_slices 
		// 	  << " goodslices=" << vol.num_good_slices 
		// 	  << " fracgood=" << vol.frac_good_slices 
		// 	  << " planecharge=[" << vol.plane_charge[0] << "," << vol.plane_charge[1] << "," << vol.plane_charge[2] << "]"
		// 	  << std::endl;
		// }

		return vols;
	}

	void ClusterGroupMatchingAlgo::GenPreMatches( const std::vector<PlaneClusterGroups>& plane_groups, ClusterGroupMatchingAlgo::AlgoData_t& data ) {
		// make pre-matches
		// factorial fun!

		Combinator< ClusterGroup > combo( plane_groups );

		do {
			std::vector< const ClusterGroup* > groupcombo = combo.getCombo();
			PreMatchMetric_t prematch( *groupcombo[0], *groupcombo[1], *groupcombo[2] );
			prematch.m_index_combo = combo.getIndexCombo();
			//std::cout << " (" << prematch.m_index_combo[0] << "," << prematch.m_index_combo[1] << "," << prematch.m_index_combo[2] << ")" << std::endl;
			data.prematch_combos_v.emplace_back( std::move(prematch) );
			combo.next();
		} while ( !combo.isLast() );

		std::sort( data.prematch_combos_v.begin(), data.prematch_combos_v.end() );
	}

	PreMatchMetric_t::PreMatchMetric_t( const ClusterGroup& p1, const ClusterGroup& p2, const ClusterGroup& p3 ) { // assumes 3-plane match

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

	ChargeVolume ClusterGroupMatchingAlgo::GetIntersectionVolume( const std::vector<larcv::Image2D>& untagged_v, const PreMatchMetric_t& prematch ) {

		bool debug_verbose = false;
		if ( m_debug_targetcombo.size()>0 )
			debug_verbose = true;

		if ( debug_verbose )
  		std::cout << "ClusterGroupMatchingAlgo::GetIntersectionVolume" << std::endl;

		const larcv::ImageMeta& meta = untagged_v.front().meta();

	  std::vector<int> row_interval = GetCommonRowInterval( prematch );
	  if ( debug_verbose ) {
	  	std::cout << "combo[" << prematch.m_index_combo[0] << "," << prematch.m_index_combo[1] << "," << prematch.m_index_combo[2] << "] ";
	    std::cout << " row interval: [" << row_interval[0] << "," << row_interval[1] << "] "
     	  << " ticks: [" << meta.pos_y(row_interval[1]) << "," << meta.pos_y(row_interval[0]) << "]"
     	  << " :: "
     	  << " [" << meta.pos_y(prematch.m_clusters[0]->tick_end) << "," << meta.pos_y(prematch.m_clusters[0]->tick_start) << "] "
     	  << " [" << meta.pos_y(prematch.m_clusters[1]->tick_end) << "," << meta.pos_y(prematch.m_clusters[1]->tick_start) << "] "
     	  << " [" << meta.pos_y(prematch.m_clusters[2]->tick_end) << "," << meta.pos_y(prematch.m_clusters[2]->tick_start) << "] "
     	  << std::endl;
    }
    if ( row_interval[0]>row_interval[1] ) {
    	//std::cout << "  no time overlap" << std::endl;
    	ChargeVolume empty;
    	return empty;
    }

   	int rows = abs( row_interval[1]-row_interval[0] );

   	int slice_size = 4;
   	int ntime_slices = rows/slice_size;
   	if ( rows%slice_size!=0 ) ntime_slices++;

   	SliceList_t area_slices;
   	std::vector<float> tot_plane_charge( untagged_v.size(), 0.0 );

   	for (int islice=0; islice<ntime_slices; islice++) {

   		Slice_t theslice;

   		// define the row slice
   		theslice.row_interval.resize(2);
   		theslice.row_interval[0] = row_interval[0] + slice_size*(islice);
   		theslice.row_interval[1] = theslice.row_interval[0] + slice_size;

   		// define the wire intervals over this slice
   		WireIntervalList_t wire_intervals = GetWireInterval( prematch, theslice.row_interval[0], theslice.row_interval[1], meta );

   		if ( debug_verbose ) {
   		std::cout << " slice[" << islice << "] "
   			<< "tick=[" << meta.pos_y( theslice.row_interval[1] ) << "," << meta.pos_y(theslice.row_interval[0]) << "]"
   			<< " wire intervals: "
   			<< "[" << wire_intervals[0][0] << "," << wire_intervals[0][1] << "] "
   			<< "[" << wire_intervals[1][0] << "," << wire_intervals[1][1] << "] "
   			<< "[" << wire_intervals[2][0] << "," << wire_intervals[2][1] << "] ";
   		}

   		// get the (y,z) point from the center of each wire interval
   		theslice.centroid = GetIntersectionCentroid( wire_intervals );
   		if ( debug_verbose )
    		std::cout << "centroid: (" << theslice.centroid[0] << "," << theslice.centroid[1] << ") ";
   		if ( theslice.centroid[0]<=-10000.0 && theslice.centroid[1]<=-10000.0 ) {
   			if ( debug_verbose )
   			 	std::cout << " [bad slice]" << std::endl;
   			continue;
   		}

   		// get a list of all 12 intersection points made up by the intervals
   		PointList_t yzintersections = GetIntersectionPoints( wire_intervals );
   		//std::cout << "nxsections=" << yzintersections.size() << " ";

   		// we filter out the intersection points along the (y,z) overlap polygon and form the boundary of the polygon
   		PointList_t yzboundary = GetBoundaryPoints( yzintersections, theslice.centroid, wire_intervals );
   		//std::cout << "boundarypts=" << yzboundary.size() << " ";

   		// we make sure the the defined boundary is within the TPC active region (y,z)
   		theslice.inside_tpc_boundary = EnforceTPCBounds( yzboundary );
   		if ( debug_verbose )
   			std::cout << " inside tpc boundary pt=" << theslice.inside_tpc_boundary.size();
   		// we calculate the area of this polygon

   		// we can use the polygon to narrow the wire-ranges
   		theslice.wire_intervals = RecalculateWireIntervalsFromBoundary( theslice.inside_tpc_boundary );
   		if ( debug_verbose ) {
     		std::cout << "overlap intervals: ";
     		for ( auto const& interval : theslice.wire_intervals ) {
   		  	std::cout << "[" << interval[0] << "," << interval[1] << "] ";
   		  }
   		}

   		// within the narrowed slice, we sum pixel charge
   		theslice.plane_charge = SumContainedCharge( prematch, untagged_v, theslice.wire_intervals, theslice.row_interval[0], theslice.row_interval[1] );
   		for (size_t p=0; p<theslice.plane_charge.size(); p++)
   			tot_plane_charge[p] += theslice.plane_charge[p];

   		// store the slice
   		area_slices.emplace_back( std::move(theslice) );

   		if ( debug_verbose )
   		  std::cout << std::endl;
   	}

   	// Build the ChargeVolume Data product
   	// This is poorly coupled. Rethink this.
	  int num_good_slices = 0;
	  for ( auto const& slice : area_slices ) {
	  	if ( slice.inside_tpc_boundary.size()>=3 )
	  		num_good_slices++;
	  }
	  float frac_good_slices = float(num_good_slices)/float(area_slices.size());
		ChargeVolume vol;
		vol.num_good_slices = num_good_slices;
		vol.num_slices = area_slices.size();			
		if ( vol.num_slices>0 ) {
  		vol.frac_good_slices = frac_good_slices;
		}
  	else {
  		vol.frac_good_slices = 0.;
  	}
		vol.m_clustergroups = prematch.m_clusters;
		vol._clustergroup_indices = prematch.m_index_combo;
		vol.m_plane_pixels = vol.GetPixelsInsideVolume( untagged_v );

		std::swap( vol.slices, area_slices );			
		std::swap( vol.plane_charge, tot_plane_charge );

   	// volume is the set of slices
   	return vol;
	}

	std::vector<int> ClusterGroupMatchingAlgo::GetCommonRowInterval( const PreMatchMetric_t& prematch ) {
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

	std::vector< WireInterval > ClusterGroupMatchingAlgo::GetWireInterval( const PreMatchMetric_t& prematch, 
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

	PointList_t ClusterGroupMatchingAlgo::GetIntersectionPoints( const std::vector< WireInterval >& plane_wire_intervals ) {

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

	Point_t ClusterGroupMatchingAlgo::GetIntersectionCentroid( const std::vector< WireInterval >& plane_wire_intervals ) {
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

	PointList_t ClusterGroupMatchingAlgo::GetBoundaryPoints( const PointList_t& crossingpts, 
		const Point_t& centroid, const std::vector<WireInterval>& wireranges  ) {

		// note: if we have only 2 good wire intervals, the 4 intersection points will be boundary points by definition as we don't have 
		// the additional constraint to use
		// we simply copy back if so.
		int numbadplanes = 0;
		for ( size_t p=0; p<wireranges.size(); p++) {
			if ( wireranges[p][0]<0 || wireranges[p][1]<0 ) {
				numbadplanes++;
			}
		}
		if ( numbadplanes>0 )
			return crossingpts; // a simple copy

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

	std::vector<WireInterval> ClusterGroupMatchingAlgo::RecalculateWireIntervalsFromBoundary( const PointList_t& yzboundary ) {
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

	PointList_t ClusterGroupMatchingAlgo::EnforceTPCBounds( const PointList_t& yzboundary ) {
		// In this method, we modify the boundary, if needed, to stay within the wire plane boundary.
		// we look for the first internal yz-boundary point. From there we step until we found a point that goes out of bounds. We find the intersection 
		// point between that line segment and the TPC boundary. 
		// we then go to the next out-of-bounds to in-bounds transition point and do the same.  
		// These boundary points are the inserted into the new yz boundary 
		// We do this until we circl back to the first in-bound point;

		PointList_t inside_boundary_points;

		int idx = 0;
		int first_inbound_idx = -1;
		typedef enum { kTop=0, kDownstream, kBot, kUpstream } TPCBoundary_t; // these are the boundaries we can cross
		std::vector< std::vector<float> > boundaries[4];
		for (int i=0; i<4; i++) {
			boundaries[i].resize(2);
			boundaries[i][0].resize(2);
			boundaries[i][1].resize(2);			
	  }

		// top
		boundaries[kTop][0][1] = 0;
		boundaries[kTop][0][0] = 118.0;
		boundaries[kTop][1][1] = 1037.0;
		boundaries[kTop][1][0] = 118.0;

		// bottom
		boundaries[kBot][0][1] = 0;
		boundaries[kBot][0][0] = -118.0;
		boundaries[kBot][1][1] = 1037.0;
		boundaries[kBot][1][0] = -118.0;

		// kDownstream
		boundaries[kDownstream][0][1] = 1037;
		boundaries[kDownstream][0][0] = -118.0;
		boundaries[kDownstream][1][1] = 1037;
		boundaries[kDownstream][1][0] = 118.0;

		// kUpstream
		boundaries[kUpstream][0][1] = 0.0;
		boundaries[kUpstream][0][0] = -118.0;
		boundaries[kUpstream][1][1] = 0.0;
		boundaries[kUpstream][1][0] = 118.0;


		// find the first in-bound point
		for ( size_t ipt=0; ipt<yzboundary.size(); ipt++ ) {
			auto const& pt = yzboundary.at(ipt);
			if ( pt[1]>0 && pt[1]<1037.0 && pt[0]>-118.0 && pt[0]<118.0 ) {
				first_inbound_idx = ipt;
				break;
			}
		}

		if ( first_inbound_idx<0 ) {
			//std::cout << "no internal boundary points." << std::endl;
			return inside_boundary_points;
		}

		inside_boundary_points.push_back( yzboundary.at(first_inbound_idx) );

		// now move through boundary points finding inbound/outbound transitions (and vice-versa)
		typedef enum { kIn=0, kOut } SearchState_t;
		SearchState_t state = kIn;
		idx = first_inbound_idx+1;
		int last_inbound_idx = first_inbound_idx;
		int last_outbound_idx = 0;
		while ( idx!=first_inbound_idx) {
			if ( idx>=(int)yzboundary.size()) {
				// wrapped back around.
				idx = 0;
				continue;
			}
			const Point_t& pt = yzboundary.at(idx);
			SearchState_t thispt_state = kIn;
			if ( pt[1]<0 || pt[1]>1037.0 || pt[0]<-118.0 || pt[0]>118.0 ) {
				thispt_state = kOut;
			}

			// in-to-out transition
			if ( thispt_state!=state ) {

				// find intersection point
				// we'll use the UBWireTool Line segment tool
				Point_t last_pt(0,0,0,0);
				if ( state==kIn )
					last_pt = yzboundary.at(last_inbound_idx);
				else
					last_pt = yzboundary.at(last_outbound_idx);
				Point_t outpt(0,0,0,0);
				if (thispt_state==kOut)
					outpt = pt;
				else
					outpt = last_pt;

				std::vector< std::vector<float> > lineseg;
				lineseg.push_back( last_pt.as_vec() );
				lineseg.push_back( pt.as_vec() );
				std::vector<float>  intersection(2,-1.0);
				int crosses = 0;

				// easy sectors
				if ( outpt[1]<0 && outpt[0]>-118.0 && outpt[0]<118.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kUpstream], intersection, crosses );
				}
  			else if ( outpt[1]>1037.0 && outpt[0]>-118.0 && outpt[0]<118.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kDownstream], intersection, crosses );
  			}
  			else if ( outpt[0]> 118.0 && outpt[1]>0 && outpt[1]<1037.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kTop], intersection, crosses );
  			}
  			else if ( outpt[0]<-118.0 && outpt[1]>0 && outpt[1]<1037.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kBot], intersection, crosses );
  			}
  			// corner sectors  			
  			else if ( outpt[1]<0 && outpt[0]>118.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kUpstream], intersection, crosses );
  				if ( crosses==0 )
    				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kTop], intersection, crosses );
  			}
  			else if ( outpt[1]<0 && outpt[0]<-118.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kUpstream], intersection, crosses );
  				if ( crosses==0 )
    				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kBot], intersection, crosses );
  			}
  			else if ( outpt[1]>1037.0 && outpt[0]>118.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kDownstream], intersection, crosses );
  				if ( crosses==0 )
    				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kTop], intersection, crosses );  				
  			}
  			else if ( outpt[1]>1037.0 && outpt[0]<-118.0 ) {
  				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kDownstream], intersection, crosses );
  				if ( crosses==0 )
    				larcv::UBWireTool::lineSegmentIntersection2D( lineseg, boundaries[kBot], intersection, crosses );
  			}

  			if (crosses==0) {
  			  std::cout << "in/out transition[" << idx << "]: thispt[" << thispt_state << "]="" << (" << pt[0] << "," << pt[1] << ") "
  				   << "--- (" << intersection[0] << "," << intersection[1] << ") ---"
  				   << " state[" << state << "]=("  << last_pt[0] << "," << last_pt[1] << ") "
  				   << "crosses=" << crosses  << " last_in_idx=" << last_inbound_idx << " last_outbound_idx=" << last_outbound_idx
  				   << std::endl;
  				throw std::runtime_error("in-out boundary crossing missed.");
  			}

  			// switch state
  			// std::cout << "in/out transition[" << idx << "]: (" << pt[0] << "," << pt[1] << ") "
  			//   << "--- (" << intersection[0] << "," << intersection[1] << ") ---"
  			//   << "("  << last_pt[0] << "," << last_pt[1] << ") " 
  			//   << "crosses=" << crosses << std::endl;
  			if ( thispt_state==kOut )
  				last_outbound_idx = idx;
  			else
  				last_inbound_idx = idx;
  			state = thispt_state;
  			// add to vector
  			Point_t boundary_pt( intersection[0], intersection[1], pt.tickstart, pt.tickend );
  			inside_boundary_points.emplace_back( std::move(boundary_pt) );
			}
			else if ( thispt_state==kIn ) {
				// add to vector
				inside_boundary_points.push_back( pt );
				last_inbound_idx = idx;
			}
			else if ( thispt_state==kOut ) {
				// not interested
				last_outbound_idx = idx;
			}
			idx++;
		}

		return inside_boundary_points;
	}

  std::vector<float> ClusterGroupMatchingAlgo::SumContainedCharge( const PreMatchMetric_t& prematch, const std::vector<larcv::Image2D>& untagged_v, 
  	const std::vector<WireInterval>& overlap_intervals, const int row_start, const int row_end ) {

  	std::vector<float> plane_charge(prematch.m_clusters.size(),0.0);

  	for ( size_t p=0; p<prematch.m_clusters.size(); p++ ) {
  		const larcv::ImageMeta& meta = untagged_v.at(p).meta();
			const ClusterGroup& group = *( prematch.m_clusters.at(p) );
			const WireInterval& winterval = overlap_intervals.at(p);

			float plane_cluster_charge = 0.;

			for (auto const& pixels : group.m_clusters_v ) {
				for ( auto const& pix : pixels ) {
					if ( (int)pix.Y()>=row_start && (int)pix.Y()<=row_end
						&& meta.pos_x(pix.X())>=winterval[0] && meta.pos_x(pix.X())<=winterval[1] )
						plane_cluster_charge += untagged_v.at(p).pixel( pix.Y(), pix.X() );
				}
			}

			plane_charge[p] = plane_cluster_charge;
  	}

  	return plane_charge;

  }

#ifndef __CINT__
#ifdef USE_OPENCV
  void ClusterGroupMatchingAlgo::labelCVImageWithMatchedClusters( std::vector<cv::Mat>& cvimgs, const std::vector<larcv::Image2D>& img_v, 
  	const std::vector<ChargeVolume>& vols, const float frac_good_threshold  ) {

    TRandom rand(1);
    for ( auto const& vol : vols ){

      if ( vol.frac_good_slices<frac_good_threshold )
        continue;

      // pick a color at random
      cv::Vec3b color;
      color[0] = (int)(rand.Uniform()*255);
      color[1] = (int)(rand.Uniform()*255);
      color[2] = (int)(rand.Uniform()*255);

      for ( size_t p=0; p<cvimgs.size(); p++ ) {
        cv::Mat& cvimg = cvimgs.at(p);

        // tag the pixesl
        const larlitecv::ClusterGroup& group = *( vol.m_clustergroups.at(p) );
        const larcv::ImageMeta& meta = img_v.at(p).meta();

        for ( auto const& slice : vol.slices ) {

          const larlitecv::WireInterval& winterval     = slice.wire_intervals.at(p);
          const std::vector<int>& rinterval = slice.row_interval;
          for ( auto const& pixels : group.m_clusters_v ) {
            for ( auto const& pix : pixels ) {

              if ( (int)meta.pos_x(pix.X())>=winterval[0] && (int)meta.pos_x(pix.X())<=winterval[1] 
                && (int)pix.Y()>=rinterval[0] && (int)pix.Y()<=rinterval[1] ) {
                cvimg.at<cv::Vec3b>( cv::Point(pix.X(), pix.Y()) ) = color;
              }
            }
          }
        }
      }
    }
  }	    	
#endif
#endif

}