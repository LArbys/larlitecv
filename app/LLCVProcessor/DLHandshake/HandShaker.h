#ifndef __HANDSHAKER_H__
#define __HANDSHAKER_H__

#include <iostream>

#include "Base/GeoTypes.h"

#include "DataFormat/pfpart.h"
#include "DataFormat/vertex.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/gtruth.h"
#include "DataFormat/simch.h"

#include "DataFormat/EventPGraph.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventPixel2D.h"
#include "HandShakeUtils.h"

namespace llcv {

  class HandShaker{
    
  public:
    
    HandShaker(){ this->reset(); }
    ~HandShaker(){}

    void set_larlite_pointers(larlite::event_pfpart*  ev_pfpart,
			      larlite::event_vertex*  ev_vertex,
			      larlite::event_shower*  ev_shower,
			      larlite::event_track*   ev_track,
			      larlite::event_cluster* ev_cluster,
			      larlite::event_hit*     ev_hit,
			      larlite::event_ass*     ev_ass);

    void construct_contour(const larcv::EventPGraph&  ev_pgraph,
			   const larcv::EventPixel2D& ev_pixel2d,
			   const larlite::event_hit*  ev_hit);
        

    void construct_image(const larcv::EventPGraph&  ev_pgraph,
			 const larcv::EventPixel2D& ev_pixel2d,
			 const larlite::event_hit*  ev_hit);


    void reset();

    void pixel_distance_threshold(double dist)
    { _dist_thresh = -1. * dist; }

    bool ready() const;

    void set_pfparticle_types(const std::vector<int>& ptype_v) { _ptype_v = ptype_v; }

  protected:

    double _dist_thresh;
    
    std::vector<int> _ptype_v;

    larlite::AssSet_t _ass_pfpart_to_vertex;  // many to 1
    larlite::AssSet_t _ass_vertex_to_track;   // 1 to many
    larlite::AssSet_t _ass_vertex_to_shower;   // 1 to many
    larlite::AssSet_t _ass_pfpart_to_track;   // 1 to 1
    larlite::AssSet_t _ass_pfpart_to_shower;  // 1 to 1
    larlite::AssSet_t _ass_pfpart_to_cluster;  // 1 to many
    larlite::AssSet_t _ass_track_to_cluster;  // 1 to many
    larlite::AssSet_t _ass_shower_to_cluster; // 1 to many
    larlite::AssSet_t _ass_track_to_hit;  // 1 to many
    larlite::AssSet_t _ass_shower_to_hit; // 1 to many
    larlite::AssSet_t _ass_cluster_to_hit;    // 1 to many
    larlite::event_pfpart*  _ev_pfpart;
    larlite::event_vertex*  _ev_vertex;
    larlite::event_shower*  _ev_shower;
    larlite::event_track*   _ev_track;
    larlite::event_cluster* _ev_cluster;
    larlite::event_hit*     _ev_hit;
    larlite::event_ass*     _ev_ass;

  public:

    void copy_here_to_there(larlite::event_base* ev_in,larlite::event_base* ev_out) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

