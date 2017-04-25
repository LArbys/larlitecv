#include "RadialEndpointFilter.h"

// std lib
#include <cmath>

// larlitecv
#include "ChargeSegmentAlgos/Segment3DAlgoTypes.h"
#include "ChargeSegmentAlgos/RadialSegmentSearch.h"

namespace larlitecv {

  bool RadialEndpointFilter::isWithinStraightSegment( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
						      const float segment_radius, const float segment_min_width, const float segment_frac_w_charge,
						      const std::vector<float>& pixel_thresholds, const float acceptance_angle, int& num_seg3d ) {

    larlitecv::RadialSegmentSearch segalgo;
    std::vector< Segment3D_t > seg3d_v = segalgo.find3Dsegments( img_v, badch_v, pos3d, segment_radius, pixel_thresholds, segment_min_width, segment_frac_w_charge );
    num_seg3d = (int)seg3d_v.size();

    if ( seg3d_v.size()<=1 ) {
      //std::cout << " less than or only 1 3d segment: " << seg3d_v.size() << std::endl;
      return false;
    }

    bool found_straight_angle = false; // use this to mark if found and to break out of loop
    for (int i=0; i<(int)seg3d_v.size(); i++) {
      if ( found_straight_angle ) break;
      for (int j=i+1; j<(int)seg3d_v.size(); j++) {
	if ( found_straight_angle ) break;

	
	Segment3D_t& seg_i = seg3d_v[i];
	Segment3D_t& seg_j = seg3d_v[j];

	// we need to know which end of the segment (start or end) corresponds to the pos3d
	float dist_start_i = 0;
	float dist_start_j = 0;
	for (int v=0; v<1; v++) {
	  dist_start_i += (seg_i.start[v]-pos3d[v])*(seg_i.start[v]-pos3d[v]);
	  dist_start_j += (seg_j.start[v]-pos3d[v])*(seg_j.start[v]-pos3d[v]);
	}
	dist_start_i = sqrt(dist_start_i);
	dist_start_j = sqrt(dist_start_j);	
	//std::cout << "dist start: i=" << dist_start_i << " j=" << dist_start_j << std::endl;
	
	if ( dist_start_i>0.4 ) {
	  //reverse (i)
	  float tmp[3] = { seg_i.start[0], seg_i.start[1], seg_i.start[2] };
	  for (int v=0; v<3; v++) {
	    seg_i.start[v] = seg_i.end[v];
	    seg_i.end[v]   = tmp[v];
	  }
	}
	if ( dist_start_j>0.4 ) {
	  //reverse (i)
	  float tmp[3] = { seg_j.start[0], seg_j.start[1], seg_j.start[2] };
	  for (int v=0; v<3; v++) {
	    seg_j.start[v] = seg_j.end[v];
	    seg_j.end[v]   = tmp[v];
	  }
	}

	// get norm vectors for direction
	float dir_i[3] = {0};
	float dir_j[3] = {0};
	float dist_i = 0;
	float dist_j = 0;
	for (int v=0; v<3; v++) {
	  //dir_i[v] = seg_i.end[v]-pos3d[v];
	  //dir_j[v] = seg_j.end[v]-pos3d[v];
	  dir_i[v] = seg_i.end[v]-seg_i.start[v];
	  dir_j[v] = seg_j.end[v]-seg_j.start[v];
	  dist_i += dir_i[v]*dir_i[v];
	  dist_j += dir_j[v]*dir_j[v];
	}
	dist_i = sqrt(dist_i);
	dist_j = sqrt(dist_j);
	for (int v=0; v<3; v++) {
	  dir_i[v] /= dist_i;
	  dir_j[v] /= dist_j;
	}

	float cosdir = 0.;
	for (int v=0; v<3; v++)
	  cosdir += dir_i[v]*dir_j[v];	
	//std::cout << "acceptance angle(" << i << "," << j << "): "
	//	  << " cos=" << cosdir
	//	  << " " << fabs(cosdir-(-1.0))*180.0/3.14159 << " deg. out of " << seg3d_v.size() << " 3d segments." << std::endl;
	if ( std::fabs( cosdir-(-1.0) )<acceptance_angle ) {
	  //std::cout << "found straight line instance: " << fabs(cosdir-(-1.0))*180.0/3.14159 << " deg. out of " << seg3d_v.size() << " 3d segments." << std::endl;
	  found_straight_angle = true;
	}
	
      }
    }
    
    return found_straight_angle;
  }

  bool RadialEndpointFilter::isWithinStraightSegment( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
						      const RadialEndpointFilterConfig& config, int& num_seg3d ) {
    return isWithinStraightSegment( pos3d, img_v, badch_v, config.segment_radius, config.segment_min_width, config.segment_frac_w_charge, config.pixel_thresholds, config.acceptance_angle, num_seg3d );
  }
  

  
}
