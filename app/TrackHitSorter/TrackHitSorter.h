#ifndef __TRACKHITSORTER_H__
#define __TRACKHITSORTER_H__

#include <vector>

// larlite
#include "larlite/core/DataFormat/vertex.h"
#include "larlite/core/DataFormat/track.h"
#include "larlite/core/DataFormat/hit.h"

// Geo2D
#include "Geo2D/Core/Geo2D.h"
#include "Geo2D/Core/LineSegment.h"


namespace larlitecv {

  class HitOrder {
  public:
  HitOrder() : phit(NULL), s(0), r(0) {};
    HitOrder ( const larlite::hit* phit_, float s_, float r_ ) : phit(phit_), s(s_), r(r_) {};
    ~HitOrder() {};
    
    const larlite::hit* phit; // pointer to hit
    float s; // path position metric
    float r; // distance from track seg to hit

    bool operator< ( const HitOrder& rh ) const {
      if ( s < rh.s ) return true;
      return false;
    };
    
  };
  
  class TrackHitSorter {
    
  public:

  TrackHitSorter() : _ismc(true) {};
    ~TrackHitSorter(){};
    
    void buildSortedHitList( const larlite::vertex& vtx, const larlite::track& track, const std::vector<larlite::hit>& hit_v,
			     const float max_radius, std::vector<int>& hitmask_v );
    void getPathBinneddEdx( const float binstep, const float binwidth, std::vector< std::vector<float> >& dedx_per_plane );
    void dump() const;
    float q2MeV( const float q, const std::vector<float>& xyz );
    const std::vector< std::vector<float> >& getBinCentersXYZ( int generating_plane ) { return bincenters_xyz[generating_plane]; };

    // track segments. 3d and 2d projected
    std::vector< std::vector<float> > path3d[3]; // per plane. corresponding 3d point at a
    std::vector< float > dist3d[3]; // per plane. corresponding 3d point at a
    std::vector< geo2d::LineSegment<float> > seg_v[3]; // segment per plane
    std::vector< float > segdist_v[3]; // distance to the segment

    std::vector< std::vector<float> > bincenters_xyz[3]; // per plane. 3d position of bin centers over witch we calculated de/dx
    std::vector<HitOrder> pathordered[3]; // per plane. ordered by path length
    std::vector<HitOrder> distordered[3]; // per plane. ordered by distance from vertex

    void clear();

    void SetIsMC(bool ismc) { _ismc = ismc; }

  private:
    bool _ismc;
    
  };


}


#endif
