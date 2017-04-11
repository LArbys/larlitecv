#ifndef __BM_TRACK_CLUSTER_3D__
#define __BM_TRACK_CLUSTER_3D__

/*
  This class aggregates BoundarySpacePoint end points and a 3D path between them.
  What's the point of this object?

 */

#include <vector>

// larlite
#include "DataFormat/track.h"

// LArCV
#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2DCluster.h"

//  larlitecv/app/Thrumu
#include "BoundaryMuonTaggerTypes.h"
#include "BoundarySpacePoint.h"


namespace larlitecv {

  class BMTrackCluster3D {

  public:

    BMTrackCluster3D(); // default constructor. empty object
    BMTrackCluster3D( const BoundarySpacePoint& startpt, const BoundarySpacePoint& endpt, const std::vector< std::vector<double> >& path3d );
    virtual ~BMTrackCluster3D();

    bool isempty() const;
    int tick_start() const;
    int tick_end() const;

    bool operator< (const BMTrackCluster3D& rhs ) const {
      if ( tick_start()<rhs.tick_start() ) return true;
      if ( tick_start()==rhs.tick_start() && tick_end()<rhs.tick_end() ) return true;
      return false;
    };

    bool operator== (const BMTrackCluster3D& rhs ) const  {
      if ( tick_start()==rhs.tick_start() && tick_end()==rhs.tick_end() )
        return true;
      return false;
    };

    void markImageWithTrack( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
			     const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
			     std::vector<larcv::Image2D>& markedimgs, const float stepsize, const float markvalue );

    const std::vector<larcv::Pixel2DCluster>& getTrackPixelsFromImages( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
									const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
									const bool recalc, const float stepsize );


    larlite::track makeTrack() const;

    BoundarySpacePoint start_endpts;
    BoundarySpacePoint end_endpts;
    std::vector< std::vector<double> > path3d;
    std::vector< larcv::Pixel2DCluster > plane_pixels; // only filled if asked for

  };


}

#endif
