#include "T3DPCMerge.h"

namespace larlitecv {

  T3DPCMerge::T3DPCMerge() {
    m_max_endmerge_dist = 30.0;
    m_min_pcacos = 0.8;
    m_max_iterations = 100;
    m_verbose = 0;
  }

  std::vector<T3DCluster> T3DPCMerge::merge( const std::vector<T3DCluster>& tracks ) {

    // techniques we use:
    //  (1) look at shortest distance between primary principle components. if small and cosine between them is below a threshold, then merge
    //  (2) if track is short, end is close to longer track, and pca intersects middle of other track, add it as subtrack
    std::vector<T3DCluster> endpoint_merged = endPointMerge(tracks);

    //std::vector<T3DCluster> submerged = mergeSubTracks( endpoint_merged );

    return endpoint_merged;

  }

  std::vector< T3DCluster > T3DPCMerge::endPointMerge( const std::vector<T3DCluster>& tracks ) {
    bool performed_merge = false;
    int iter=0;

    std::vector<T3DCluster> tracklist;
    // we copy
    for (auto& track : tracks ) {
      tracklist.push_back( track );
    }

    std::vector<T3DCluster*> current_list;
    for (auto& track : tracklist )
      current_list.push_back( &track );

    do {
      performed_merge = false;

      std::vector<T3DCluster*> merge_list;

      for (int i=0; i<(int)current_list.size(); i++) {
      	for (int j=i+1; j<(int)current_list.size(); j++) {

      	  T3DCluster& track_i = *(current_list[i]);
      	  T3DCluster& track_j = *(current_list[j]);

          // --------------------------------------------------------------------
          // Check quality of tracks

          // size: tracks need to be made up of more than one point
          if ( track_i.getPath().size()<2 )
            continue;
          if ( track_j.getPath().size()<2 )
            continue;


      	  // --------------------------------------------------------------------
      	  // END POINT MERGE
      	  double closest_dist = 1.0e6;
      	  std::vector<int> whichends;
	  if ( m_verbose>1 )
	    std::cout << "Testing (" << i << "," << j << ")" << std::endl;
      	  bool valid = shouldWeEndPointMerge( track_i, track_j, closest_dist, whichends );

      	  if ( valid  ) {
      	    performed_merge = true;
      	    if ( whichends[0]==0 && whichends[1]==0 ) {
      	      // start/start
      	      track_i.reverse();
      	      track_i.append( track_j );
      	      track_i.reverse();
      	    }
      	    else if ( whichends[0]==1 && whichends[1]==0 ) {
      	      // end/start
      	      track_i.append( track_j );
      	    }
      	    else if ( whichends[0]==0 && whichends[1]==1 ) {
      	      // start/end
      	      track_j.append( track_i );
      	      std::swap( track_i, track_j ); // dangerous?
      	    }
      	    else if ( whichends[0]==1 && whichends[1]==1 ) {
      	      //end/end
      	      track_j.reverse();
      	      track_i.append( track_j );
      	    }
      	  }

      	  // --------------------------------------------------------------------
      	  // GATHER MERGED TRACKS
      	  if ( performed_merge ) {
      	    // fill the merge list to go to the next round
      	    merge_list.push_back( &track_i );
      	    for ( auto& ptrack : current_list ) {
      	      if ( ptrack==&track_i || ptrack==&track_j )
      		continue;
      	      merge_list.push_back( ptrack );
      	    }
      	    break;
      	  }
      	}

      	if ( performed_merge )
      	  break;
      }

      // update track list
      if ( performed_merge ) {
      	current_list.clear();
      	current_list = merge_list;
      }

      iter++;
    } while ( performed_merge && iter<m_max_iterations );


    // copy current track list objects to output vector
    std::vector<T3DCluster> out;
    for ( auto& ptrack : current_list ) {
      out.push_back( *ptrack );
    }
    return out;
  }

  bool T3DPCMerge::shouldWeEndPointMerge( const T3DCluster& ta, const T3DCluster& tb, double& closest_dist, std::vector<int>& whichends ) {
    
    if ( m_verbose>1 ) {
      std::cout << "----------------------" << std::endl;
      std::cout << " should we merge" << std::endl;
      std::cout << "  TA size: " << ta.getPath().size() << std::endl;
      std::cout << "  TB size: " << tb.getPath().size() << std::endl;
      std::cout << "  TA: (" << ta.getPath().front()[0] << "," << ta.getPath().front()[1] << "," << ta.getPath().front()[2] << ") -> "
    		<< "(" << ta.getPath().back()[0] << "," << ta.getPath().back()[1] << "," << ta.getPath().back()[2] << ")"
    		<< std::endl;
      std::cout << "  TB: (" << tb.getPath().front()[0] << "," << tb.getPath().front()[1] << "," << tb.getPath().front()[2] << ") -> "
    		<< "(" << tb.getPath().back()[0] << "," << tb.getPath().back()[1] << "," << tb.getPath().back()[2] << ")"
    		<< std::endl;

      // std::cout << "  TA PCA0: (" << linea.Pt1()[0] << "," << linea.Pt1()[1] << "," << linea.Pt1()[2] << ") -> "
      // 		<< "(" << linea.Pt2()[0] << "," << linea.Pt2()[1] << "," << linea.Pt2()[2] << ")"
      // 		<< std::endl;
      // std::cout << "  TB PCA0: (" << lineb.Pt1()[0] << "," << lineb.Pt1()[1] << "," << lineb.Pt1()[2] << ") -> "
      // 		<< "(" << lineb.Pt2()[0] << "," << lineb.Pt2()[1] << "," << lineb.Pt2()[2] << ")"
      // 		<< std::endl;

      std::cout << "  TA pc-0 dir: (" << ta.getPCADir(0)[0] << "," << ta.getPCADir(0)[1] << "," << ta.getPCADir(0)[2] << ")" << std::endl;
      std::cout << "  TB pc-0 dir: (" << tb.getPCADir(0)[0] << "," << tb.getPCADir(0)[1] << "," << tb.getPCADir(0)[2] << ")" << std::endl;
      std::cout << "  TA mean: (" << ta.getMean()[0] << "," << ta.getMean()[1] << "," << ta.getMean()[2] << ")" << std::endl;
      std::cout << "  TB mean: (" << tb.getMean()[0] << "," << tb.getMean()[1] << "," << tb.getMean()[2] << ")" << std::endl;

      // std::cout << "  Mid point between close end: (" << midpt[0] << "," << midpt[1] << "," << midpt[2] << ") " << std::endl;
      // std::cout << "  closest point distance=" << closest_dist << " cm" << std::endl;
      std::cout << "--[TB]-------------" << std::endl;
      for (int i=0; i<(int)tb.getPath().size(); i++) {
	const auto& tbpt = tb.getPath()[i];
	std::cout << "  #" << i << ": (" << tbpt[0] << "," << tbpt[1] << "," << tbpt[2] << ")" << std::endl;
      }
      
    }
    //return false;// bail!
    
    // Make PCA line
    std::vector<double> staa(3);
    std::vector<double> stab(3);
    std::vector<double> enda(3);
    std::vector<double> endb(3);
    for (int i=0; i<3; i++) {
      staa[i] = ta.getMean()[i]+ta.getPCABounds(0)[0]*ta.getPCADir(0)[i];
      enda[i] = ta.getMean()[i]+ta.getPCABounds(0)[1]*ta.getPCADir(0)[i];

      stab[i] = tb.getMean()[i]+tb.getPCABounds(0)[0]*tb.getPCADir(0)[i];
      endb[i] = tb.getMean()[i]+tb.getPCABounds(0)[1]*tb.getPCADir(0)[i];
    }

    // return 'false' and do not merge the endpoints if one of the lines have the same endpoints, meaning that they are just a single point.
    // use 0.002 cm for endpoints that differ slightly from rounding.
    if ( ( fabs(staa[0] - enda[0]) < .002 && fabs(staa[1] - enda[1]) < 0.002 && fabs(staa[2] - enda[2]) < 0.002 ) || ( fabs(stab[0] - endb[0]) < 0.002 && fabs(stab[1] - endb[1]) < 0.002 && fabs(stab[2] - endb[2]) < 0.002 ) ) {
      return false;
    }
    
    ::larlite::geoalgo::Line_t linea( staa, enda );
    ::larlite::geoalgo::Line_t lineb( stab, endb );

    // Distance between end points
    double enddist[4] = {0};
    for (int i=0; i<3; i++) {
      enddist[0] += (ta.getPath().front()[i]-tb.getPath().front()[i])*(ta.getPath().front()[i]-tb.getPath().front()[i]); // start to start
      enddist[1] += (ta.getPath().front()[i]-tb.getPath().back()[i])*(ta.getPath().front()[i]-tb.getPath().back()[i]); // start to end
      enddist[2] += (ta.getPath().back()[i]-tb.getPath().front()[i])*(ta.getPath().back()[i]-tb.getPath().front()[i]); // end to start
      enddist[3] += (ta.getPath().back()[i]-tb.getPath().back()[i])*(ta.getPath().back()[i]-tb.getPath().back()[i]); // end to end
    }
    enum EndType_t { kStartStart=0, kStartEnd, kEndStart, kEndEnd, kNumEndTypes };
    EndType_t endtype = kStartStart;
    double endmin = enddist[0];
    for (int i=kStartEnd;i<kNumEndTypes; i++) {
      if ( enddist[i]<endmin ) {
      	endtype = (EndType_t)i;
      	endmin = enddist[i];
      }
    }
    endmin = sqrt(endmin);
    if ( m_verbose>1 ) {
      std::cout << "End point match Type: " << endtype << std::endl;
      std::cout << "Distance between closest ends: " << endmin << " cm" << std::endl;
    }

    // pass back, which ends to use
    whichends.resize(2);

    // set which ends. also calculate the mid point between the two ends
    ::larlite::geoalgo::Point_t midpt(3);
    for (int i=0; i<3; i++) {
      switch (endtype) {
      case kStartStart:
      	midpt[i] = 0.5*(ta.getPath().front()[i]+tb.getPath().front()[i]);
      	whichends[0] = 0;
      	whichends[1] = 0;
      	break;
      case kStartEnd:
        midpt[i] = 0.5*(ta.getPath().front()[i]+tb.getPath().back()[i]);
        whichends[0] = 0;
        whichends[1] = 1;
        break;
      case kEndStart:
      	midpt[i] = 0.5*(ta.getPath().back()[i]+tb.getPath().front()[i]);
      	whichends[0] = 1;
      	whichends[1] = 0;
      	break;
      case kEndEnd:
      	midpt[i] = 0.5*(ta.getPath().back()[i]+tb.getPath().back()[i]);
      	whichends[0] = 1;
      	whichends[1] = 1;
      	break;
      case kNumEndTypes:
      default:
      	throw std::runtime_error("bad end type");
      	break;
      }//end of switch
      // cruft: in case we do this using ends of the first-PC
      // case 0:
      // 	midpt[i] = 0.5*(linea.Pt1()[i]+lineb.Pt1()[i]);
      // 	l2pt[i] = 0.5*(linea.Pt1()[i]+lineb.Pt1()[i]);
      // 	break;
      // case 1:
      // 	midpt[i] = 0.5*(linea.Pt1()[i]+lineb.Pt2()[i]);
      // 	l2pt[i] = 0.5*(linea.Pt1()[i]+lineb.Pt2()[i]);
      // 	break;
      // case 2:
      // 	midpt[i] = 0.5*(linea.Pt2()[i]+lineb.Pt1()[i]);
      // 	l2pt[i] = 0.5*(linea.Pt2()[i]+lineb.Pt1()[i]);
      // 	break;
      // case 3:
      // 	midpt[i] = 0.5*(linea.Pt2()[i]+lineb.Pt2()[i]);
      // 	l2pt[i] = 0.5*(linea.Pt2()[i]+lineb.Pt2()[i]);
      // 	break;
      // default:
      // 	throw std::runtime_error("bad end type");
      // 	break;
      // }//end of switch
    }
    closest_dist = 0.;

    if ( m_verbose>1 ) {
      // std::cout << "TA: (" << ta.getPath().front()[0] << "," << ta.getPath().front()[1] << "," << ta.getPath().front()[2] << ") -> "
      // 		<< "(" << ta.getPath().back()[0] << "," << ta.getPath().back()[1] << "," << ta.getPath().back()[2] << ")"
      // 		<< std::endl;
      // std::cout << "TB: (" << tb.getPath().front()[0] << "," << tb.getPath().front()[1] << "," << tb.getPath().front()[2] << ") -> "
      // 		<< "(" << tb.getPath().back()[0] << "," << tb.getPath().back()[1] << "," << tb.getPath().back()[2] << ")"
      // 		<< std::endl;

      std::cout << "TA PCA0: (" << linea.Pt1()[0] << "," << linea.Pt1()[1] << "," << linea.Pt1()[2] << ") -> "
    		<< "(" << linea.Pt2()[0] << "," << linea.Pt2()[1] << "," << linea.Pt2()[2] << ")"
    		<< std::endl;
      std::cout << "TB PCA0: (" << lineb.Pt1()[0] << "," << lineb.Pt1()[1] << "," << lineb.Pt1()[2] << ") -> "
    		<< "(" << lineb.Pt2()[0] << "," << lineb.Pt2()[1] << "," << lineb.Pt2()[2] << ")"
    		<< std::endl;

      // std::cout << "TA pc-0 dir: (" << ta.getPCADir(0)[0] << "," << ta.getPCADir(0)[1] << "," << ta.getPCADir(0)[2] << ")" << std::endl;
      // std::cout << "TB pc-0 dir: (" << tb.getPCADir(0)[0] << "," << tb.getPCADir(0)[1] << "," << tb.getPCADir(0)[2] << ")" << std::endl;
      // std::cout << "TA mean: (" << ta.getMean()[0] << "," << ta.getMean()[1] << "," << ta.getMean()[2] << ")" << std::endl;
      // std::cout << "TB mean: (" << tb.getMean()[0] << "," << tb.getMean()[1] << "," << tb.getMean()[2] << ")" << std::endl;

      std::cout << "Mid point between close end: (" << midpt[0] << "," << midpt[1] << "," << midpt[2] << ") " << std::endl;
      std::cout << "closest point distance=" << closest_dist << " cm" << std::endl;
    }

    // Calculate cosine between PCA directions
    float pcacos = 0.;
    for (int i=0; i<3; i++)
      pcacos += ta.getPCADir(0)[i]*tb.getPCADir(0)[i];

    // correct sign based on track orientations
    switch (endtype) {
    case kStartStart:
    case kEndEnd:
      pcacos *= -1;
      break;
    case kStartEnd:
    case kEndStart:
    case kNumEndTypes:
      break;
    }

    if ( m_verbose>1 )
      std::cout << "PCA cosine: " << pcacos << std::endl;

    // midpt and l2pt need to extend beyond
    double pca_l1 = ta.getPCvalue(midpt,0);
    double pca_l2 = tb.getPCvalue(midpt,0);
    bool outside_segment = false;
    if ( m_verbose>1 ) {
      std::cout << "PtA PC-0 value: " << pca_l1 << " bounds: [" << ta.getPCABounds(0)[0] << "," << ta.getPCABounds(0)[1] << "]" << std::endl;
      std::cout << "PtB PC-0 value: " << pca_l2 << " bounds: [" << tb.getPCABounds(0)[0] << "," << tb.getPCABounds(0)[1] << "]" << std::endl;
    }

    // dist to midpoint
    float middist_a = sqrt(m_geoalgo.SqDist( midpt, linea )); // not reasonable
    float middist_b = sqrt(m_geoalgo.SqDist( midpt, lineb ));

    if ( m_verbose>1 ) {
      std::cout << "TA: dist to midpoint=" << middist_a << " cm" << std::endl;
      std::cout << "TB: dist to midpoint=" << middist_b << " cm" << std::endl;
    }

    if ( (pca_l1<ta.getPCABounds(0)[0] || pca_l1>ta.getPCABounds(0)[1] )
      && (pca_l2<tb.getPCABounds(0)[0] || pca_l2>tb.getPCABounds(0)[1] ) ) {
      outside_segment = true;
    }
    if ( m_verbose>1 )
      std::cout << "Is midpoint outside segment: " << outside_segment << std::endl;


    if ( m_verbose>1 )
      std::cout << "End Matching Code: A=" << whichends[0] << " B=" << whichends[1] << std::endl;

    // should we match
    bool match = false;
    if ( endmin<m_max_endmerge_dist && pcacos>m_min_pcacos && outside_segment )
      match = true;

    if ( m_verbose>1 )
      std::cout << "MATCH=" << match << std::endl;

    return match;
  }

  void T3DPCMerge::mergeTracks( T3DCluster& track_i, T3DCluster& track_j, const std::vector<int>& whichends ) {
    if ( whichends[0]==0 && whichends[1]==0 ) {
      // start/start
      track_i.reverse();
      track_i.append( track_j );
      track_i.reverse();
    }
    else if ( whichends[0]==1 && whichends[1]==0 ) {
      // end/start
      track_i.append( track_j );
    }
    else if ( whichends[0]==0 && whichends[1]==1 ) {
      // start/end
      track_j.append( track_i );
      std::swap( track_i, track_j ); // dangerous?
    }
    else if ( whichends[0]==1 && whichends[1]==1 ) {
      //end/end
      track_j.reverse();
      track_i.append( track_j );
    }
  }

}
