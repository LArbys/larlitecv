
#include "FoxTrotTrackerAlgo.h"

namespace larlitecv {

  FoxTrotTrackerAlgo::FoxTrotTrackerAlgo( const FoxTrotTrackerAlgoConfig& cfg )
    : m_config(cfg) {
  }

  FoxTrack FoxTrotTrackerAlgo::followTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
					    const BoundarySpacePoint& start ) {
    FoxTrack track;

    // need first step dir. define based on end
    std::vector<float> firstdir(3,0);
    switch ( start.type() ) {
    case kTop:
      firstdir[1] = -1.0;
      break;
    case kBottom:
      firstdir[1] = -1.0;
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

    FoxStep firststep( start.pos(), firstdir );
    track.push_back( firststep );
    std::cout << "start of fox track: (" << firststep.pos()[0] << "," << firststep.pos()[1] << "," <<  firststep.pos()[2] << ")" << std::endl;
    float min_dcos = -1.0;
    do {
      if ( track.size()>0 )
	min_dcos = 0.5;
      FoxStep next = getNextStep( track.back(), img_v, badch_v, min_dcos );
      if ( next.isgood() )
	std::cout << "next fox step: (" << next.pos()[0] << "," << next.pos()[1] << "," << next.pos()[2] << ") good=" << next.isgood() << std::endl;
      track.emplace_back( std::move(next) );
    } while ( track.back().isgood() && track.size()<100 );

    // pop off the end which is bad
    if ( !track.back().isgood() )
      track.pop_back();

    return track;
  }


  FoxStep FoxTrotTrackerAlgo::getNextStep( const FoxStep& current, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, float min_dcos ) {
    // we take a step by using the radialsegment search algo.
    // we find a 3d segment that is on the straightest path for us
    // we know it doesn't always work out, so we shrink the step a few times

    float current_radius = m_config.step_size;

    for (int iattempt=0; iattempt<m_config.num_step_attempts; iattempt++) {

      std::vector< Segment3D_t > candidate_segs = m_radialalgo.find3Dsegments( img_v, badch_v, current.pos(), current_radius, m_config.pixel_thresholds,
									       m_config.min_hit_width, m_config.segment_frac_w_charge );
      std::cout << "attempt=" << iattempt << " numsegs=" << candidate_segs.size() << std::endl;
      //if ( candidate_segs.size()>2 )
      //continue;
      int ibest = -1;
      float bestdcos = 2.;
      std::vector<float> bestdir(3,0);
      std::vector<float> bestpos(3,0);
      for ( size_t iseg=0; iseg<candidate_segs.size(); iseg++) {
        Segment3D_t& seg = candidate_segs[iseg];

        // first, which is the near and far points?
        float dist_start=0;
        float dist_end=0;
        for (int v=0; v<3; v++) {
          dist_start += ( current.pos()[v]-seg.start[v] )*( current.pos()[v]-seg.start[v] );
          dist_end   += ( current.pos()[v]-seg.end[v] )*( current.pos()[v]-seg.end[v] );
        }
        if ( dist_start > dist_end ) {
          // flip
          float tmp[3] = { seg.start[0], seg.start[1], seg.start[2] };
          for (int v=0; v<3; v++) {
            seg.start[v] = seg.end[v];
            seg.end[v] = tmp[v];
          }
        }

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
	std::cout << "  seg " << iseg << " dcos=" << dcos << std::endl;
        if ( dcos >min_dcos ) {
	  if ( dcos>bestdcos || ibest<0 ) {
	    ibest = iseg;
	    bestdcos = dcos;
	    bestdir = canddir;
	    for (int v=0; v<3; v++)
	      bestpos[v] = seg.end[v];
	  }
	}
      }

      if ( ibest<0 ) {
        current_radius *= m_config.radius_reduction_factor;
      }
      else {
        return FoxStep( bestpos, bestdir );
      }

    }//end of attempt loop
    
    /// no step :(. Return empty step.
    FoxStep empty;
    std::cout << "return empty" << std::endl;
    return empty;

  }

}
