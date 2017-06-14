#include "FoxTrotLead.h"

#include <cmath>

namespace larlitecv {

  bool FoxTrotLead::chooseBestSegment( const FoxStep& current, std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past,
				       const FoxTrotTrackerAlgoConfig& config, const std::vector<larcv::Image2D>& stepped_v, int& best_seg_idx ) {
    // set initial best seg
    best_seg_idx = -1; // set to sentinal value

    // first condition segments. Start should be from current position
    for ( size_t iseg=0; iseg<seg3d_v.size(); iseg++) {
      Segment3D_t& seg = seg3d_v[iseg];

      // first, which is the near and far points?
      float dist_start=0;
      float dist_end=0;
      if ( seg.start[0]==seg.end[0] ) {
        for (int v=0; v<3; v++) {
          dist_start += ( current.pos()[v]-seg.start[v] )*( current.pos()[v]-seg.start[v] );
          dist_end   += ( current.pos()[v]-seg.end[v] )*( current.pos()[v]-seg.end[v] );
        }
        dist_start = sqrt(dist_start);
        dist_end   = sqrt(dist_end);
      }
      else {
        for (int v=0; v<1; v++) {
          dist_start += ( current.pos()[v]-seg.start[v] )*( current.pos()[v]-seg.start[v] );
          dist_end   += ( current.pos()[v]-seg.end[v] )*( current.pos()[v]-seg.end[v] );
        }
        dist_start = sqrt(dist_start);
        dist_end   = sqrt(dist_end);
      }

      if ( config.verbosity>1 ) {
        std::cout << "    iseg " << iseg << ": "
                  << "q(" << seg.plane_frac_w_charge[0] << "," << seg.plane_frac_w_charge[1] << "," << seg.plane_frac_w_charge[2] << ")"
                  << " start=(" << seg.start[0] << "," << seg.start[1] << "," << seg.start[2] << ")"
		              << " end=(" << seg.end[0] << "," << seg.end[1] << "," << seg.end[2] << ")"
		              << std::endl;
      }

      if ( dist_start > dist_end ) {
        // flip
        if ( config.verbosity>1 )
          std::cout << "    flip seg start/end. startdist=" << dist_start << " enddist=" << dist_end << std::endl;
        float tmp[3] = { seg.start[0], seg.start[1], seg.start[2] };
        for (int v=0; v<3; v++) {
          seg.start[v] = seg.end[v];
          seg.end[v] = tmp[v];
        }
      }
    }

    return _chooseBestSegment_( current, seg3d_v, past, config, stepped_v, best_seg_idx );
  }


}
