
#include "FoxTrotTrackerAlgo.h"

namespace larlitecv {

  FoxTrotTrackerAlgo::FoxTrotTrackerAlgo( const FoxTrotTrackerAlgoConfig& cfg )
    : m_config(cfg) {
    m_user_lead = NULL;
  }

  FoxTrotTrackerAlgo::~FoxTrotTrackerAlgo() {
  }

  FoxTrack FoxTrotTrackerAlgo::followTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
					    const BoundarySpacePoint& start ) {
    // function where the caller doesn't specify a starting direction
    std::vector<float> start_dir(1,0);
    return followTrack( img_v, badch_v, tagged_v, start, start_dir );
  }

  FoxTrack FoxTrotTrackerAlgo::followTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
					    const BoundarySpacePoint& start, const std::vector<float>& start_dir ) {
    FoxTrack track;

    // need first step dir. define based on end
    std::vector<float> firstdir(3,0);

    if ( start_dir.size()!=3 ) {
      // if user doesn't provide a starting direction, we infer one based on the boundary crossing
      switch ( start.type() ) {
      case kTop:
        firstdir[1] = -1.0;
        break;
      case kBottom:
        firstdir[1] = 1.0;
        break;
      case kUpstream:
        firstdir[2] = 1.0;
        break;
      case kDownstream:
        firstdir[2] = -1.0;
        break;
      case kAnode:
        firstdir[0] = 1.0;
        break;
      case kCathode:
        firstdir[0] = -1.0;
        break;
      case kImageEnd:
        if (start.pos()[0]>0)
          firstdir[0] = -1.0;
        else
          firstdir[0] = 1.0;
        break;
      default:
        throw std::runtime_error("unrecognize boundary type");
        break;
      }//end of switch
    }
    else {
      firstdir = start_dir;
    }

    // Define First Step
    FoxStep firststep( start.pos(), firstdir );
    track.push_back( firststep );

    if ( m_config.verbosity>0 )
      std::cout << "start of fox track: (" << firststep.pos()[0] << "," << firststep.pos()[1] << "," <<  firststep.pos()[2] << ")" << std::endl;
    float min_dcos = -1.0;
    do {
      // if we've moved beyond the first step, we enforce a forward-going track
      if ( track.size()>0 )
        min_dcos = m_config.min_cosine;
      m_default_lead.setMinCos( min_dcos );
      FoxStep next = getNextStep( track.back(), img_v, badch_v, track );
      if ( m_config.verbosity>0 && next.isgood() )
        std::cout << "next fox step: (" << next.pos()[0] << "," << next.pos()[1] << "," << next.pos()[2] << ") "
        	  << " dir=(" << next.dir()[0] << "," << next.dir()[1] << "," << next.dir()[2] << ") "
        	  << " good=" << next.isgood() << std::endl;
      else if ( m_config.verbosity>0 && !next.isgood() )
        std::cout << "next fox step is bad." << std::endl;
      track.emplace_back( std::move(next) );
    } while ( track.back().isgood() && track.size()<m_config.max_steps );

    // pop off the end which is bad
    if ( !track.back().isgood() )
      track.pop_back();

    return track;
  }


  FoxStep FoxTrotTrackerAlgo::getNextStep( const FoxStep& current, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const FoxTrack& path) {

    // we take a step by using the radialsegment search algo.
    // we find a 3d segment that is on the straightest path for us
    // we know it doesn't always work out, so we shrink the step a few times

    float current_radius = m_config.step_size;

    if ( m_config.verbosity>1 ) {
      std::cout << "fox step ---------------" << std::endl;
      std::cout << "  current pos: (" << current.pos()[0] << "," << current.pos()[1] << "," << current.pos()[2] << ")" << std::endl;
    }
    for (int iattempt=0; iattempt<m_config.num_step_attempts; iattempt++) {

      // Use the radial algo to get candidate 3D segments
      std::vector< Segment3D_t > candidate_segs = m_radialalgo.find3Dsegments( img_v, badch_v, current.pos(), current_radius, m_config.pixel_thresholds,
									       m_config.min_hit_width, m_config.hit_neighborhood, m_config.segment_frac_w_charge, m_config.verbosity );
      if ( m_config.verbosity>1 ) {
        std::cout << "  attempt=" << iattempt << " numsegs=" << candidate_segs.size() << std::endl;
      }

      bool goodseg = false;
      int ibest = -1;
      if ( m_user_lead!=0 ) {
        goodseg = m_user_lead->chooseBestSegment( current, candidate_segs, path, m_config, ibest );
      }
      else {
        goodseg = m_default_lead.chooseBestSegment( current, candidate_segs, path, m_config, ibest );
      }

      if ( !goodseg || ibest<0 ) {
        current_radius *= m_config.radius_reduction_factor;
      }
      else {
        Segment3D_t& seg = candidate_segs[ibest];
        if ( m_config.verbosity>1 ) {
          std::cout << "    chosen best segment idx=" << ibest
        	    << " (" << seg.start[0] << "," << seg.start[1] << "," << seg.start[2] << ") -> "
        	    << " (" << seg.end[0] << "," << seg.end[1] << "," << seg.end[2] << ")"
        	    << std::endl;
        }
        std::vector<float> canddir(3,0);
        std::vector<float> candpos(3,0);
        float dist = 0.;
        for (int v=0; v<3; v++) {
          candpos[v] = seg.end[v];
          canddir[v] = seg.end[v]-seg.start[v];
          dist += canddir[v]*canddir[v];
        }
        dist = sqrt(dist);
        for (int v=0; v<3; v++)
          canddir[v] /= dist;

        return FoxStep( candpos, canddir );
      }

    }//end of attempt loop

    /// no step :(. Return empty step.
    FoxStep empty;
    //std::cout << "return empty" << std::endl;
    return empty;

  }

  bool FoxTrotLeadStraight::_chooseBestSegment_( const FoxStep& current, std::vector<Segment3D_t>& candidate_segs, const FoxTrack& past,
						 const FoxTrotTrackerAlgoConfig& config, int& best_seg_idx ) {

    best_seg_idx = -1;
    float bestdcos = 2.;
    std::vector<float> bestdir(3,0);
    std::vector<float> bestpos(3,0);
    if ( config.verbosity>1 ) {
      std::cout << "    FoxTrotLeadStraight::chooseBestSegment" << std::endl;
    }
    for ( size_t iseg=0; iseg<candidate_segs.size(); iseg++) {
      Segment3D_t& seg = candidate_segs[iseg];

      // // first, which is the near and far points?
      // float dist_start=0;
      // float dist_end=0;
      // for (int v=0; v<1; v++) {
      // 	dist_start += ( current.pos()[v]-seg.start[v] )*( current.pos()[v]-seg.start[v] );
      // 	dist_end   += ( current.pos()[v]-seg.end[v] )*( current.pos()[v]-seg.end[v] );
      // }
      // dist_start = sqrt(dist_start);
      // dist_end   = sqrt(dist_end);

      // if ( config.verbosity>1 ) {
      // 	std::cout << "    iseg " << iseg << ": start=(" << seg.start[0] << "," << seg.start[1] << "," << seg.start[2] << ")"
      // 		  << " end=(" << seg.end[0] << "," << seg.end[1] << "," << seg.end[2] << ")"
      // 		  << std::endl;
      // }

      // if ( dist_start > dist_end ) {
      // 	// flip
      // 	if ( config.verbosity>1 )
      // 	  std::cout << "    flip seg start/end. startdist=" << dist_start << " enddist=" << dist_end << std::endl;
      // 	float tmp[3] = { seg.start[0], seg.start[1], seg.start[2] };
      // 	for (int v=0; v<3; v++) {
      // 	  seg.start[v] = seg.end[v];
      // 	  seg.end[v] = tmp[v];
      // 	}
      // }

      std::vector<float> canddir(3,0);
      float dist = 0.;
      for (int v=0; v<3; v++) {
        canddir[v] = seg.end[v]-seg.start[v];
        dist += canddir[v]*canddir[v];
      }
      dist = sqrt(dist);
      for (int v=0; v<3; v++)
        canddir[v] /= dist;
      float cosseg = 0;
      for (int v=0; v<3; v++) {
        cosseg += canddir[v]*current.dir()[v];
      }
      float dcos = cosseg;
      if ( config.verbosity>1 )
      std::cout << "    seg " << iseg << " dcos=" << dcos
      	  << " (" << seg.start[0] << "," << seg.start[1] << "," << seg.start[2] << ") -> "
      	  << " (" << seg.end[0] << "," << seg.end[1] << "," << seg.end[2] << ")"
      	  << " segdir=(" << canddir[0] << "," << canddir[1] << "," << canddir[2] << ")"
      	  << " current=(" << current.dir()[0] << "," << current.dir()[1] << "," << current.dir()[2] << ")"
      	  << std::endl;
      if ( dcos >min_dcos ) {
        if ( dcos>bestdcos || best_seg_idx<0 ) {
          best_seg_idx = iseg;
          bestdcos = dcos;
          bestdir = canddir;
          for (int v=0; v<3; v++)
            bestpos[v] = seg.end[v];
        }
      }
    }

    return true;
  }

}
