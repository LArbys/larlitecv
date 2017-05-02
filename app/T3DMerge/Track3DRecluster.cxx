#include "Track3DRecluster.h"

#include "GeoAlgo/GeoAABox.h"
#include "GeoAlgo/GeoTrajectory.h"
#include "GeoAlgo/GeoAlgo.h"

#include "ANN/ANNAlgo.h"

namespace larlitecv {

  void Track3DRecluster::addPath( const std::vector< std::vector<float> >& path ) {
    T3DCluster::Builder builder;
    for ( auto const& pt : path ) {
      std::vector<double> v(3);
      for (int i=0; i<3; i++)
	v[i] = pt[i];
      builder.addPoint( v );
    }

    T3DCluster track = builder.build();
    m_tracks.emplace_back( std::move(track) );
  }

  void Track3DRecluster::addPath( const std::vector< std::vector<double> >& path ) {
    T3DCluster::Builder builder;
    builder.setPath( path );
    T3DCluster track = builder.build();
    m_tracks.emplace_back( std::move(track) );
  }

  void Track3DRecluster::recluster() {
    int iterations = 0;
    int num_reclustered = 0;
    do {
      std::vector< T3DCluster > reclustered_tracks;
      num_reclustered = 0;
      for (int itracka=0; itracka<(int)m_tracks.size(); itracka++) {
	for (int itrackb=itracka+1; itrackb<(int)m_tracks.size(); itrackb++) {
	  bool reclustered = ReclusterPair( m_tracks[itracka], m_tracks[itrackb], reclustered_tracks );
	  if (reclustered)
	    num_reclustered++;
	}
      }
      iterations++;
    } while ( num_reclustered==0 || iterations>=10 );
    
  }

  bool Track3DRecluster::ReclusterPair( const T3DCluster& tracka, const T3DCluster& trackb, std::vector<T3DCluster>& tracks_v ) {
					
    // we do a broad test for overlap using bounding box
    // for those that do, we:
    // (1) we use a kdtree to find out if any point on track b is near points on track a
    // (2) likewise for the reverse

    
    const geoalgo::AABox& boxa = tracka.getBBox();
    const geoalgo::AABox& boxb = trackb.getBBox();

    if ( boxa.Contain( boxb.Min() ) || boxa.Contain(boxb.Max() )
	 || boxb.Contain( boxa.Min() ) || boxb.Contain( boxa.Max() ) ) {
      // no overlap. just return.
      tracks_v.push_back( tracka );
      tracks_v.push_back( trackb );
      return false;
    }

    // make fragments. first make tracks based on overlapping regions
    // then for non overlapping regions, build track fragments
    std::vector<T3DCluster> fragments_v;
    
    // we check the overlap of track a with track b
    std::vector<int> point_has_overlap_a;
    std::vector< SegmentOverlap_t > segoverlap_a = getOverlapSegmentsOfBonA( trackb, tracka, point_has_overlap_a );

    if ( segoverlap_a.size()==0 ) {
      // no meaningful overlap
      tracks_v.push_back( tracka );
      tracks_v.push_back( trackb );
      return false;
    }
    
    // we build a track from the overlap
    std::vector<int> useda( tracka.getPath().size(), 0 );
    std::vector<int> usedb( trackb.getPath().size(), 0 );
    
    for ( auto const& overlap : segoverlap_a ) {

      std::vector< std::vector<int> > insert_idx_a( overlap.indices.size() );
      
      // we make a trajectory out of the segment
      geoalgo::Trajectory seg_traj;
      for ( auto const& idx_a : overlap.indices ) {
	geoalgo::Vector pt_a = tracka.getPath()[idx_a];
	seg_traj.push_back( pt_a );
      }

      // we loop over overlaping points in track b
      for (auto const& idx_b : overlap.othertrackmatches ) {
	geoalgo::Point_t pt_b = trackb.getPath()[idx_b];
	int idx_a;
	geoalgo::Point_t traj_pt = m_geoalgo.ClosestPt( seg_traj, pt_b, idx_a );
	std::vector<double> dir_b2a(3);
	float dist_b2a = 0.;
	for (int i=0; i<3; i++) {
	  dir_b2a[i] = tracka.getPath()[idx_a][i] - pt_b[i];
	  dist_b2a += dir_b2a[i]*dir_b2a[i];
	}
	dist_b2a = sqrt( dist_b2a );

	std::vector<double> dir_pta(3,0);
	if ( idx_a<(int)tracka.getPath().size()-1 ) {
	  dir_pta = tracka.getPathDir()[idx_a];
	}
	else {
	  dir_pta = tracka.getPathDir().back();
	}
	
	
	float cosb2a = 0;
	for (int i=0; i<3; i++) {
	  cosb2a += dir_pta[i]*(dir_b2a[i]/dist_b2a);
	}

	int insert_after_idxa = idx_a;
	if ( cosb2a < 0 ) {
	  insert_after_idxa--;
	}
	idx_a = ( idx_a<0 ) ? 0 : idx_a;
	idx_a = ( idx_a>=(int)tracka.getPath().size() ) ? idx_a-1 : idx_a;
	insert_idx_a.at(idx_a).push_back( idx_b );
      }//loop over matches to track-b points in track-b segment

      std::vector< std::vector<double> > newpath;
      for (auto const& idx_a : overlap.indices ) {
	useda[idx_a] = 1;
	newpath.push_back( tracka.getPath()[idx_a] );
	for (auto const& idx_b : insert_idx_a[idx_a] ) {
	  newpath.push_back( trackb.getPath()[idx_b] );
	  usedb[idx_b] = 1;
	}	  
      }
      T3DCluster::Builder builder;
      builder.setPath( newpath ).buildDirList().updateBBox();
      T3DCluster t3d = builder.build();

      fragments_v.emplace_back( std::move(t3d) );
      
    }// loop over segments

    // build track out of unused points
    bool inseg = false;
    T3DCluster::Builder builder;
    for (int idx_a=0; idx_a<(int)useda.size(); idx_a++) {

      if ( !inseg ) {
	if ( useda[idx_a]!=0 ) {
	  inseg = true;
	  builder.clear();
	  builder.addPoint( tracka.getPath()[idx_a] );
	}
      }
      else {
	if ( useda[idx_a]==1 )
	  builder.addPoint( tracka.getPath()[idx_a] );
	else {
	  inseg = false;
	  fragments_v.push_back( builder.build() );
	}
      }
    }

    for (int idx_b=0; idx_b<(int)usedb.size(); idx_b++) {

      if ( !inseg ) {
	if ( useda[idx_b]!=0 ) {
	  inseg = true;
	  builder.clear();
	  builder.addPoint( trackb.getPath()[idx_b] );
	}
      }
      else {
	if ( usedb[idx_b]==1 )
	  builder.addPoint( trackb.getPath()[idx_b] );
	else {
	  inseg = false;
	  fragments_v.push_back( builder.build() );
	}
      }
    }

    // last step, we reconnect fragments on the ends
    // better would be to see if one can fox trot past the ends...
    enum EndPair_t { kStartStart=0, kStartEnd, kEndStart, kEndEnd };
    std::vector<int> used_fragment( fragments_v.size(), 0 );

    for (int frag_a=0; frag_a<(int)fragments_v.size(); frag_a++) {
      if ( used_fragment[frag_a]==1 )
	continue;
      for (int frag_b=frag_a+1; frag_b<(int)fragments_v.size(); frag_b++) {
	if ( used_fragment[frag_b]==1 )
	  continue;
	
	T3DCluster& fa = fragments_v[frag_a];
	T3DCluster& fb = fragments_v[frag_b];
	
	// compare end points and their direction
	float dist2[4]  = {0};
	float cosend[4] = {0};
	
	for (int i=0; i<3; i++) {
	  // distance
	  dist2[kStartStart] += ( fa.getPath().front()[i] - fb.getPath().front()[i] )*( fa.getPath().front()[i] - fb.getPath().front()[i] );
	  dist2[kStartEnd]   += ( fa.getPath().front()[i] - fb.getPath().back()[i] )*(  fa.getPath().front()[i] - fb.getPath().back()[i] );
	  dist2[kEndStart]   += ( fa.getPath().back()[i]  - fb.getPath().front()[i] )*( fa.getPath().back()[i]  - fb.getPath().front()[i] );
	  dist2[kEndEnd]     += ( fa.getPath().back()[i]  - fb.getPath().back()[i] )*(  fa.getPath().back()[i]  - fb.getPath().back()[i] );
	  // cosine
	  cosend[kStartStart] += fa.getPathDir().front()[i]*(-1.0*fb.getPathDir().front()[i]);
	  cosend[kStartEnd]   += fa.getPathDir().front()[i]*fb.getPathDir().back()[i];
	  cosend[kEndStart]   += fa.getPathDir().back()[i]*fb.getPathDir().front()[i];
	  cosend[kEndEnd]     += fa.getPathDir().back()[i]*(-1.0*fb.getPathDir().back()[i]);
	}
	int closest_pair = 0;
	float min_dist2   = dist2[kStartStart];
	for (int i=1; i<4; i++) {
	  if ( dist2[i]<min_dist2 ) {
	    min_dist2 = dist2[i];
	    closest_pair = i;
	  }
	}
	  
	if ( sqrt(min_dist2)<3.0 && cosend[closest_pair]>0 ) {
	  used_fragment[frag_a] = 1;
	  used_fragment[frag_b] = 1;
	  if ( closest_pair==kStartStart ) {
	    fa.reverse();
	    fa.append( fb );
	    tracks_v.emplace_back( std::move(fa) );
	  }
	  else if ( closest_pair==kStartEnd ) {
	    fb.append( fa );
	    tracks_v.emplace_back( std::move(fb) );
	  }
	  else if ( closest_pair==kEndStart ) {
	    fa.append( fb );
	    tracks_v.emplace_back( std::move(fa) );
	  }
	  else if ( closest_pair==kEndEnd ) {
	    fb.reverse();
	    fa.append( fb );
	    tracks_v.emplace_back( std::move(fa) );
	  }
	  
	}
      }
    }

    if ( tracks_v.size()==0 ) {
      // did not recluster
      tracks_v.push_back( tracka );
      tracks_v.push_back( trackb );
      return false;
    }
    
    return true;

  }
  
  std::vector< Track3DRecluster::SegmentOverlap_t > Track3DRecluster::getOverlapSegmentsOfBonA( const T3DCluster& tracka, const T3DCluster& trackb, std::vector<int>& overlap_info ) {
    overlap_info.resize( trackb.getPath().size(), 0 );
    
    // we build a kd-tree for A. then we check matches with b
    ann::ANNAlgo anne_a( tracka.getPath().size(), 3 );
    anne_a.initialize();    
    for ( int idx=0; idx<(int)trackb.getPath().size(); idx++) {
      auto const& pt = tracka.getPath()[idx];
      anne_a.setPoint( idx, pt );
    }
    
    std::vector<SegmentOverlap_t > segments;
    bool inseg = false;
    for ( int idx_b=0; idx_b<(int)trackb.getPath().size(); idx_b++ ) {
      auto const& pt = trackb.getPath()[idx_b];
      std::vector<int> matched_idx_a = anne_a.regionQuery( pt, 5.0, 0.0 );
      bool overlaps = false;
      if ( matched_idx_a.size()>0 ) {
	overlap_info[idx_b] = 1;
	overlaps = true;
      }

      if ( !inseg ) {
	if ( overlaps ) {
	  inseg = true;
	  SegmentOverlap_t newseg;
	  newseg.indices.push_back( idx_b );
	  for (auto const& idx_a : matched_idx_a )
	    newseg.othertrackmatches.insert( idx_a );
	  segments.emplace_back( std::move(newseg) );
	}
      }
      else {
	if ( !overlaps ) {
	  // turn off segment and store
	  inseg = false;
	}
	else {
	  // continue segment
	  segments.back().indices.push_back( idx_b );
	  // log a-indices
	  for (auto const& idx_a : matched_idx_a )
	    segments.back().othertrackmatches.insert( idx_a );
	}
      }
    }

    // clean up
    anne_a.deinitialize();
    ann::ANNAlgo::cleanup();

    return segments;
  }
  
}




