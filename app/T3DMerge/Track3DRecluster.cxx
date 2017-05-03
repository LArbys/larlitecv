#include "Track3DRecluster.h"

#include "GeoAlgo/GeoAABox.h"
#include "GeoAlgo/GeoTrajectory.h"
#include "GeoAlgo/GeoAlgo.h"

#include "ANN/ANNAlgo.h"

namespace larlitecv {

  Track3DRecluster::Track3DRecluster() {
    m_verbosity = 0;
    m_min_segoverlap = 2; // should be dist
    m_ann_dist = 5.0;
    m_max_end_dist = 20.0;
    m_min_end_cos = 0.0;
    m_max_iterations = 10000;
    m_merge = false;
  }
  
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

  std::vector< T3DCluster > Track3DRecluster::recluster() {
    int iterations = 0;
    int num_reclustered = 0;
    bool reclustered = false;    
    do {
      std::vector< T3DCluster > temp_list;
      int old_size = m_tracks.size();
      int used_tracka_idx = -1;
      int used_trackb_idx = -1;
      num_reclustered = 0;
      reclustered = false;      
      for (int itracka=0; itracka<(int)m_tracks.size(); itracka++) {
	for (int itrackb=itracka+1; itrackb<(int)m_tracks.size(); itrackb++) {
	  if ( m_verbosity>1 )
	    std::cout << "Iteration " << iterations << " Recluster (" << itracka << "," << itrackb << ")" << std::endl;
	  std::vector< T3DCluster > reclustered_results;
	  reclustered = ReclusterPair( m_tracks[itracka], m_tracks[itrackb], reclustered_results );
	  if (reclustered) {
	    if ( m_verbosity>1 )
	      std::cout << "  reclustering occured. turned 2 tracks into=" << reclustered_results.size() << std::endl;
	    used_tracka_idx = itracka;
	    used_trackb_idx = itrackb;
	    for ( auto& t : reclustered_results )
	      temp_list.emplace_back( std::move(t) );
	    num_reclustered++;
	    break;
	  }
	}
	if ( reclustered )
	  break;
      }
      for (int itrack=0; itrack<(int)m_tracks.size(); itrack++ ) {
	if ( itrack==used_tracka_idx || itrack==used_trackb_idx )
	  continue;
	temp_list.emplace_back( std::move(m_tracks[itrack]) );
      }
      m_tracks.clear(); // dump the tracks
      for ( auto& track : temp_list ) {
	m_tracks.emplace_back( std::move(track) );
      }
      temp_list.clear();
      if ( m_verbosity>0 )
	std::cout << "RECLUSTER iter=" << iterations << " oldsize=" << old_size << " new size=" << m_tracks.size() << std::endl;
      iterations++;
      if ( m_verbosity>2 )
	std::cin.get();
    } while ( reclustered && iterations<m_max_iterations );
    
    return m_tracks;
  }

  bool Track3DRecluster::ReclusterPair( const T3DCluster& tracka, const T3DCluster& trackb, std::vector<T3DCluster>& tracks_v ) {
					
    // we do a broad test for overlap using bounding box
    // for those that do, we:
    // (1) we use a kdtree to find out if any point on track b is near points on track a
    // (2) likewise for the reverse

    
    const geoalgo::AABox& boxa = tracka.getBBox();
    const geoalgo::AABox& boxb = trackb.getBBox();

    // bool bbox_intersect = false;
    // if ( boxa.Contain( boxb.Min() ) || boxa.Contain(boxb.Max() )
    // 	 || boxb.Contain( boxa.Min() ) || boxb.Contain( boxa.Max() ) ) {
    //   bbox_intersect = true;
    // }

    bool bbox_intersect = tracka.overlaps( trackb );
    
    if (!bbox_intersect) {
      if ( m_verbosity>1) {
	std::cout << "reclusterpair::no bounding box overlap" << std::endl;
	std::cout << "  box a (" << boxa.Min()[0] << "," << boxa.Min()[1] << "," << boxa.Min()[2] << ") "
		  << " to (" << boxa.Max()[0] << "," << boxa.Max()[1] << "," << boxa.Max()[2] << ") " << std::endl;
	std::cout << "  box b (" << boxb.Min()[0] << "," << boxb.Min()[1] << "," << boxb.Min()[2] << ") "
		  << " to (" << boxb.Max()[0] << "," << boxb.Max()[1] << "," << boxb.Max()[2] << ") " << std::endl;
      }
      return false;
    }


    // std::cout << "TRACK A" << std::endl;
    // for ( int i=0; i<(int)tracka.getPath().size(); i++ ) {
    //   auto const& pt_core = tracka.getPath()[i];
    //   std::cout << "  [" << i << "] : (" << tracka.getPath()[i][0]  << "," << tracka.getPath()[i][1]  << "," << tracka.getPath()[i][2]  << ")" << std::endl;      
    // }
    // std::cout << "TRACK B" << std::endl;
    // for ( int i=0; i<(int)trackb.getPath().size(); i++ ) {
    //   auto const& pt_core = trackb.getPath()[i];
    //   std::cout << "  [" << i << "] : (" << trackb.getPath()[i][0]  << "," << trackb.getPath()[i][1]  << "," << trackb.getPath()[i][2]  << ")" << std::endl;      
    // }
    
    // make fragments. first make tracks based on overlapping regions
    // then for non overlapping regions, build track fragments
    std::vector<T3DCluster> init_fragments_v;
    
    // we check the overlap of track a with track b
    std::vector<int> point_has_overlap_a;
    
    std::vector< SegmentOverlap_t > segoverlap_a;
    if ( tracka.getAveStepSize()<trackb.getAveStepSize() )
      segoverlap_a = getOverlapSegmentsOfBonA( tracka, trackb, point_has_overlap_a );
    else
      segoverlap_a = getOverlapSegmentsOfBonA( trackb, tracka, point_has_overlap_a );
    
    // for ( auto const& seg : segoverlap_a ) {
    //   std::cout << "segment" << std::endl;
    //   std::cout << "  segment indices: " << seg.indices.size() << std::endl;
    //   for (const auto& idx : seg.indices ) {
    // 	std::cout << "    idxa=" << idx << ": (" << tracka.getPath()[idx][0]  << "," << tracka.getPath()[idx][1]  << "," << tracka.getPath()[idx][2]  << ")" << std::endl;
    //   }
    //   std::cout << "  matched indices: " << seg.othertrackmatches.size() << std::endl;
    //   for ( const auto& idx : seg.othertrackmatches ) {
    // 	std::cout << "    idxb=" << idx << ": (" << trackb.getPath()[idx][0]  << "," << trackb.getPath()[idx][1]  << "," << trackb.getPath()[idx][2]  << ")" << std::endl;
    //   }
    // }
    

    if ( m_verbosity>1 )
      std::cout << "StepSize A=" << tracka.getAveStepSize() << " StepSize B=" << trackb.getAveStepSize() << std::endl;
    
    if ( segoverlap_a.size()==0 ) {
      // no meaningful overlap
      return false;
    }

    // we build a track from the overlap
    std::vector<int> useda( tracka.getPath().size(), 0 );
    std::vector<int> usedb( trackb.getPath().size(), 0 );
    bool merge_overlap = false;

    for ( auto const& overlap : segoverlap_a ) {
	
      // we orchestrate a combination of the two tracks over the overlap region
      T3DCluster const* core = NULL;
      T3DCluster const* transfer = NULL;
      std::vector<int> indices_core;      
      std::vector<int> indices_transfer;
      std::vector<int>* used_core = NULL;
      std::vector<int>* used_trans = NULL;

      int match_min = *(overlap.othertrackmatches.begin());
      int match_max = match_min;
      for ( auto const& idx : overlap.othertrackmatches ) {
	if ( match_max<idx )
	  match_max = idx; // easier way to get to the max?
	if ( match_min>idx )
	  match_min = idx;
      }
      for ( int idx=match_min; idx<=match_max; idx++ )
	indices_core.push_back(idx);
      indices_transfer = overlap.indices;

      if ( tracka.getAveStepSize()<trackb.getAveStepSize() ) {
	core = &tracka;
	transfer = &trackb;
	used_core = &useda;
	used_trans = &usedb;
      }
      else {
	core = &trackb;
	transfer = &tracka;
	used_core = &usedb;
	used_trans = &useda;	
      }
      
      if ( m_verbosity>1 ) {
	std::cout << "core: " << core << " size=" << core->getPath().size() << std::endl;
	std::cout << "transfer: " << transfer << " size=" << transfer->getPath().size() << std::endl;
	
	std::cout << "core: " << indices_core.size() << std::endl;
	int i=0; 
	for ( auto const& idx : indices_core ) {
	  std::cout << "    [" << i << "] idx=" << idx << ": (" << core->getPath()[idx][0]  << "," << core->getPath()[idx][1]  << "," << core->getPath()[idx][2]  << ")" << std::endl;
	  i++;
	}
	std::cout << "transfer: " << indices_transfer.size() << std::endl;      
	i=0; 
	for ( auto const& idx : indices_transfer ) {
	  std::cout << "    [" << i << "] idx=" << idx << ": (" << transfer->getPath()[idx][0]  << "," << transfer->getPath()[idx][1]  << "," << transfer->getPath()[idx][2]  << ")" << std::endl;
	  i++;
	}
      }


      if ( (int)indices_transfer.size()<m_min_segoverlap ) {
	if (m_verbosity>1) {
	  std::cout << "segment overlap too small" << std::endl;
	  if ( m_verbosity>2 )
	    std::cin.get();
	}
	continue;
      }
      
      std::vector< std::vector<int> > insert_idx_core( core->getPath().size()+1 ); // additional 0-index is the underflow
      
      // we make a trajectory out of the segment
      geoalgo::Trajectory seg_traj;
      //i=0; 
      for ( auto const& pt_core : core->getPath() ) {
	//std::cout << "core pt [" << i << "]: " << pt_core[0] << ", " << pt_core[1] << ", " << pt_core[2] << std::endl;
	seg_traj.push_back( pt_core );
	//i++;
      }
      
      // we loop over overlaping points in track b
      for (auto const& idx_trans : indices_transfer ) {
	geoalgo::Point_t pt_trans = transfer->getPath()[idx_trans];
	int idx_core = -1;
	float min_dist2 = 0;
	//geoalgo::Point_t traj_pt = m_geoalgo.ClosestPt( seg_traj, pt_trans, idx_core ); // this is not what we need, the trajectory returns the wrong segment!
	
	for ( int ipt=0; ipt<(int)core->getPath().size(); ipt++ ) {
	  float dist2 = 0;
	  for (int v=0; v<3; v++) {
	    dist2 += (pt_trans[v]-core->getPath()[ipt][v])*(pt_trans[v]-core->getPath()[ipt][v]);
	  }
	  if ( dist2<min_dist2 || idx_core<0 ) {
	    idx_core = ipt;
	    min_dist2 = dist2;
	  }
	}
	geoalgo::Point_t traj_pt = core->getPath()[idx_core];
	
	std::vector<double> dir_b2a(3);
	float dist_b2a = 0.;
	for (int i=0; i<3; i++) {
	  dir_b2a[i] = pt_trans[i] - traj_pt[i];
	  dist_b2a += dir_b2a[i]*dir_b2a[i];
	}
	dist_b2a = sqrt( dist_b2a );

	std::vector<double> dir_pta(3,0);
	if ( idx_core<(int)core->getPath().size()-1 ) {
	  dir_pta = core->getPathDir()[idx_core];
	}
	else {
	  dir_pta = core->getPathDir().back();
	}
	
	float cosb2a = 0;
	for (int i=0; i<3; i++) {
	  cosb2a += dir_pta[i]*(dir_b2a[i]/dist_b2a);
	}

	int insert_after_idxcore = idx_core; // 1-index offset
	if ( cosb2a < 0 ) {
	  insert_after_idxcore--;
	}

	if ( insert_after_idxcore<0 || insert_after_idxcore>=(int)insert_idx_core.size() )
	  continue;
	//insert_after_idxcore = ( insert_after_idxcore<0 ) ? 0 : insert_after_idxcore;
	//insert_after_idxcore = ( insert_after_idxcore>=(int)insert_idx_core.size() ) ? (int)insert_idx_core.size()-1 : insert_after_idxcore;
	if ( m_verbosity>1 ) {
	  std::cout << "insert idx_trans=" << idx_trans << " after idx_core=" << insert_after_idxcore-1 << " cosb2a=" << cosb2a << " distb2a=" << dist_b2a << std::endl;
	  std::cout << "  trans point: (" << pt_trans[0] << "," << pt_trans[1] << "," << pt_trans[2] << ")" << std::endl;
	  std::cout << "  nearby core point: idx_core=" << idx_core << " (" << traj_pt[0] << "," << traj_pt[1] << "," << traj_pt[2] << ")" << std::endl;
	}
	insert_idx_core.at(insert_after_idxcore).push_back( idx_trans );
      }//loop over matches to track-b points in track-b segment

      std::vector< std::vector<double> > newpath;

      //for (auto const& idx_a : indices_core ) {
      for (int idx_a=0; idx_a<(int)core->getPath().size(); idx_a++ ) {
	(*used_core)[idx_a] = 1;
	newpath.push_back( core->getPath()[idx_a] );
	if ( insert_idx_core[idx_a].size()>0 ) {
	  //std::cout << "inserting " << insert_idx_core[idx_a].size() << " points after idx_a=" << idx_a << std::endl;
	  for (auto const& idx_trans : insert_idx_core[idx_a] ) {
	    //std::cout << "  insert idx_trans=" << idx_trans << " after idxa=" << idx_a << std::endl;
	    if ( m_merge )
	      newpath.push_back( transfer->getPath()[idx_trans] );
	    (*used_trans)[idx_trans] = 1;
	    merge_overlap = true;
	  }
	}
      }
      T3DCluster::Builder builder;
      builder.setPath( newpath );
      T3DCluster t3d = builder.build();

      init_fragments_v.emplace_back( std::move(t3d) );
      
    }// loop over segments
    
    if ( m_verbosity>1 ) {
      std::cout << "TRACK A" << std::endl;
      for ( int i=0; i<(int)tracka.getPath().size(); i++ ) {
	std::cout << "  [" << i << "] : (" << tracka.getPath()[i][0]  << "," << tracka.getPath()[i][1]  << "," << tracka.getPath()[i][2]  << ")"
		  << " used=" << useda[i] << std::endl;      
      }
      std::cout << "TRACK B" << std::endl;
      for ( int i=0; i<(int)trackb.getPath().size(); i++ ) {
	std::cout << "  [" << i << "] : (" << trackb.getPath()[i][0]  << "," << trackb.getPath()[i][1]  << "," << trackb.getPath()[i][2]  << ")"
		  << " used=" << usedb[i] << std::endl;	
      }
    }
    
    if ( !merge_overlap ) {
      std::cout << " no overlapping fragments" << std::endl;
      return false;
    }
    
    // build track out of unused points
    bool inseg = false;
    T3DCluster::Builder builder;
    for (int idx_a=0; idx_a<(int)useda.size(); idx_a++) {

      if ( !inseg ) {
	if ( useda[idx_a]==0 ) {
	  inseg = true;
	  builder.clear();
	  builder.addPoint( tracka.getPath()[idx_a] );
	}
      }
      else {
	if ( useda[idx_a]==0 )
	  builder.addPoint( tracka.getPath()[idx_a] );
	else {
	  inseg = false;
	  init_fragments_v.push_back( builder.build() );
	}
      }
    }
    if ( inseg ) {
      // end it
      init_fragments_v.push_back( builder.build() );
      inseg = false;
    }

    inseg = false;
    for (int idx_b=0; idx_b<(int)usedb.size(); idx_b++) {

      if ( !inseg ) {
	if ( usedb[idx_b]==0 ) {
	  inseg = true;
	  builder.clear();
	  builder.addPoint( trackb.getPath()[idx_b] );
	}
      }
      else {
	if ( usedb[idx_b]==0 )
	  builder.addPoint( trackb.getPath()[idx_b] );
	else {
	  inseg = false;
	  init_fragments_v.push_back( builder.build() );
	}
      }
    }
    if ( inseg ) {
      // end it
      init_fragments_v.push_back( builder.build() );
      inseg = false;
    }

    if ( m_verbosity>1 ) {
      std::cout << "Fragments: " << init_fragments_v.size() << std::endl;
      for (int i=0; i<(int)init_fragments_v.size(); i++) {
	T3DCluster& frag = init_fragments_v[i];
	std::cout << "  Fragment " << i << std::endl;
	for (int ipt=0; ipt<(int)frag.getPath().size(); ipt++) {
	  std::cout << "     [" << ipt << "] (" << frag.getPath()[ipt][0] << "," << frag.getPath()[ipt][1] << "," << frag.getPath()[ipt][2] << ")" << std::endl;
	}
      }
    }
    
    // last step, we reconnect fragments on the ends
    // better would be to see if one can fox trot past the ends...
    enum EndPair_t { kStartStart=0, kStartEnd, kEndStart, kEndEnd };

    std::vector< T3DCluster* > fragment_list;
    for ( auto& frag : init_fragments_v )
      fragment_list.push_back( &frag );

    // we look to append fragments to one another
    // if we do it, the one fragment will have been absorbed and invalid, so we have to restart the loop.
    // we keep checking until no fragments can be merged

    bool performed_merge = true;
    int merger_iterations = 0;
    while ( performed_merge && merger_iterations<20) {
      performed_merge = false;
      
      std::vector< T3DCluster* > temp_fraglist;
      temp_fraglist.clear();
      for (int frag_a=0; frag_a<(int)fragment_list.size(); frag_a++) {
	if ( performed_merge )
	  break;

	T3DCluster& fa = *(fragment_list[frag_a]);
	if ( fa.getPath().size()<2 )
	  continue;
	  
	for (int frag_b=frag_a+1; frag_b<(int)fragment_list.size(); frag_b++) {
	  if ( performed_merge )
	    break;
	  
	  T3DCluster& fb = *(fragment_list[frag_b]);
	  // if ( fb.getPath().size()<2 )
	  //   continue;
	  
	  performed_merge = mergeEnds( fa, fb, m_max_end_dist, m_min_end_cos );
	  
	  if ( performed_merge ) {
	    temp_fraglist.push_back( fragment_list[frag_a] );	    
	    // we need to stop the merger loop, but first we store in the other fragments
	    for ( int frag_c=0; frag_c<(int)fragment_list.size(); frag_c++ ) {
	      if ( fragment_list[frag_c]==fragment_list[frag_a] || fragment_list[frag_c]==fragment_list[frag_b] )
		continue; // skip these
	      temp_fraglist.push_back( fragment_list[frag_c] );
	    }
	    break;
	  }// if merge was performed
	}//end of loop over frag_b
	if ( performed_merge )
	  break;
      }//end of loop over frag_a

      if ( performed_merge ) {
	// clear out the fragmentlist and fill it with the new one
	fragment_list.clear();
	for ( auto& pfrag : temp_fraglist ) {
	  //if ( pfrag->getPath().size()>=2 )
	  fragment_list.push_back( pfrag );
	}
      }
      if ( m_verbosity>1 )
	std::cout << "after merger iterations = " << merger_iterations << " fraglist=" << fragment_list.size() << " templist=" << temp_fraglist.size() << std::endl;
      temp_fraglist.clear();
      merger_iterations++;
    }//end of merger loop

    // move tracks over
    for ( auto& pfinal_fragments : fragment_list ) {
      if ( pfinal_fragments->getPath().size()<2 )
	continue;
      (*pfinal_fragments).makePathDir();
      (*pfinal_fragments).updateBBox();
      tracks_v.push_back( *pfinal_fragments );
    }

    bool recluster_occurred = merge_overlap | performed_merge;
    
    return recluster_occurred;

  }
  
  std::vector< Track3DRecluster::SegmentOverlap_t > Track3DRecluster::getOverlapSegmentsOfBonA( const T3DCluster& tracka, const T3DCluster& trackb, std::vector<int>& overlap_info ) {
    overlap_info.resize( trackb.getPath().size(), 0 );
    
    // we build a kd-tree for A. then we check matches with b
    ann::ANNAlgo anne_a( tracka.getPath().size(), 3 );
    for ( int idx=0; idx<(int)tracka.getPath().size(); idx++) {
      std::vector<double> pt = tracka.getPath()[idx];
      if ( pt.size()!=3 ) {
	throw std::runtime_error("track point does not have 3 dims");
      }
      
      // have to add small jitter
      for (int i=0; i<3; i++) {
	pt[i] += (m_rand.Uniform()-0.5)*1.0e-4;
      }

      anne_a.setPoint( idx, pt );
    }
    anne_a.initialize();
    
    std::vector< SegmentOverlap_t > segments;
    bool inseg = false;

    for ( int idx_b=0; idx_b<(int)trackb.getPath().size(); idx_b++ ) {
      auto const& pt = trackb.getPath()[idx_b];
      std::vector<int> matched_idx_a = anne_a.regionQuery( pt, m_ann_dist*m_ann_dist, 0.0 );

      // std::vector<int> matched_idx_a;
      // for (int idx_a=0; idx_a<(int)tracka.getPath().size(); idx_a++) {
      // 	float distab = 0;
      // 	auto const& pta = tracka.getPath()[idx_a];
      // 	for (int v=0; v<3; v++)
      // 	  distab += (pt[v]-pta[v])*(pt[v]-pta[v]);
      // 	if ( distab<m_ann_dist*m_ann_dist )
      // 	  matched_idx_a.push_back( idx_a );
      // 	std::cout << "   idx_b=" << idx_b << " idx_a=" << idx_a << " dist=" << sqrt(distab)
      // 		  << " ptb=(" << pt[0] << "," << pt[1] << "," << pt[2] << ")"
      // 		  << " pta=(" << pta[0] << "," << pta[1] << "," << pta[2] << ")"
      // 		  << std::endl;
      // }
      
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
	  for (auto const& idx_a : matched_idx_a ) {
	    segments.back().othertrackmatches.insert( idx_a );
	    //std::cout << "    insert idx_a=" << idx_a << " size=" << segments.back().othertrackmatches.size() << std::endl;
	  }
	}
      }

      // std::cout << "check overlap idxb=" << idx_b << " local_matched_idxa=" << matched_idx_a.size() << " inseg=" << inseg;
      // if ( segments.size()>0 )
      // 	std::cout << " last total_matched_idxa=" << segments.back().othertrackmatches.size();
      // std::cout << std::endl;
    }

    // clean up
    anne_a.deinitialize();
    ann::ANNAlgo::cleanup();
    if ( segments.size()>1 ) {
      int iseg = 0;
      int longest = segments.front().othertrackmatches.size();
      std::vector< SegmentOverlap_t > out;
      for (int i=1; i<(int)segments.size(); i++) {
	if ( (int)segments[i].othertrackmatches.size()>longest ) {
	  longest = segments[i].othertrackmatches.size();
	iseg = i;
	}
      }
      out.emplace_back( std::move(segments[iseg]) );
      return out;
    }else {
      return segments;
    }
  }

  bool Track3DRecluster::mergeEnds( T3DCluster& fa, T3DCluster& fb, float max_dist, float min_cos ) {

    enum EndPair_t { kStartStart=0, kStartEnd, kEndStart, kEndEnd };
    
    bool performed_merge = false;

    // compare end points and their direction
    float dist2[4]  = {0};
    float cosend[4] = {0};

    for (int i=0; i<3; i++) {

      // distance
      dist2[kStartStart] += ( fa.getPath().front()[i] - fb.getPath().front()[i] )*( fa.getPath().front()[i] - fb.getPath().front()[i] );
      dist2[kStartEnd]   += ( fa.getPath().front()[i] - fb.getPath().back()[i] )*(  fa.getPath().front()[i] - fb.getPath().back()[i] );
      dist2[kEndStart]   += ( fa.getPath().back()[i]  - fb.getPath().front()[i] )*( fa.getPath().back()[i]  - fb.getPath().front()[i] );
      dist2[kEndEnd]     += ( fa.getPath().back()[i]  - fb.getPath().back()[i] )*(  fa.getPath().back()[i]  - fb.getPath().back()[i] );

      if ( fa.getPath().size()>=2 && fb.getPath().size()>=2 ) {
	// cosine
	cosend[kStartStart] += fa.getPathDir().front()[i]*(-1.0*fb.getPathDir().front()[i]);
	cosend[kStartEnd]   += fa.getPathDir().front()[i]*fb.getPathDir().back()[i];
	cosend[kEndStart]   += fa.getPathDir().back()[i]*fb.getPathDir().front()[i];
	cosend[kEndEnd]     += fa.getPathDir().back()[i]*(-1.0*fb.getPathDir().back()[i]);
      }
      else {
	cosend[kStartStart] = 1.0;
	cosend[kStartEnd]   = 1.0;
	cosend[kEndStart]   = 1.0;
	cosend[kEndEnd]     = 1.0;
      }
    }
    
    int closest_pair = 0;
    float min_dist2   = dist2[kStartStart];
    for (int i=1; i<4; i++) {
      if ( dist2[i]<min_dist2 ) {
	min_dist2 = dist2[i];
	closest_pair = i;
      }
    }

    if ( m_verbosity>1 ) {
      std::cout << "test pair merge " << std::endl;
      std::cout << "  frag A size=" << fa.getPath().size()
		<< " s=(" << fa.getPath().front()[0] << "," << fa.getPath().front()[1] << "," << fa.getPath().front()[2] << ") "
		<< " e=(" << fa.getPath().back()[0] << ", " << fa.getPath().back()[1] << ", " << fa.getPath().back()[2] << ") ";
      if ( fa.getPath().size()>=2 ) {
	std::cout << " sdir=(" << fa.getPathDir().front()[0] << "," << fa.getPathDir().front()[1] << "," << fa.getPathDir().front()[2] << ") "
		  << " edir=(" << fa.getPathDir().back()[0] << "," << fa.getPathDir().back()[1] << "," << fa.getPathDir().back()[2] << ") ";
      }
      std::cout << std::endl;
      
      std::cout << "  frag B size=" << fb.getPath().size()
		<< " s=(" << fb.getPath().front()[0] << "," << fb.getPath().front()[1] << "," << fb.getPath().front()[2] << ") "
		<< " e=(" << fb.getPath().back()[0] << ", " << fb.getPath().back()[1] << ", " << fb.getPath().back()[2] << ") ";
      if ( fb.getPath().size()>=2 ) {
	std::cout << " sdir=(" << fb.getPathDir().front()[0] << "," << fb.getPathDir().front()[1] << "," << fb.getPathDir().front()[2] << ") "
		  << " edir=(" << fb.getPathDir().back()[0] << "," << fb.getPathDir().back()[1] << "," << fb.getPathDir().back()[2] << ") ";
      }
      std::cout << std::endl;
      std::cout << "  closest_pair=" << closest_pair << " mindist=" << sqrt(min_dist2)
		<< " cosend=" << cosend[closest_pair]
		<< std::endl;
    }
    
    if ( sqrt(min_dist2)<max_dist && cosend[closest_pair]>min_cos ) {
      if ( closest_pair==kStartStart ) {
	fa.reverse();
	fa.append(fb);
      }
      else if ( closest_pair==kStartEnd ) {
	fb.append(fa);
	std::swap(fa,fb);
      }
      else if ( closest_pair==kEndStart ) {
	fa.append( fb );
      }
      else if ( closest_pair==kEndEnd ) {
	fb.reverse();
	fa.append( fb );
      }
      performed_merge = true;
    }//end of if merger condition passes

    return performed_merge;
  }
  
}




