#ifndef __STOPMU_TRACKER__
#define __STOPMU_TRACKER__

/*! \brief Tracker for Stop Muons with attempt at 3D/2D coordination
 *         
 *
 *  Tries to step through 3D space and coordinate with hits found on 2D plane views.
 *  Once 3D direction no longer consistent, uses 2D plane views to attempt to reorient direction and continue tracking.
 *  Written specificially with stopping muons in mind.
 *  Assumes you can deliver a reasonably good track end.
 */

#include <vector>
#include <array>
#include <exception>
#include <algorithm>

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2D.h"

// larcv/app
#include "dbscan/DBSCANAlgo.h"

#include "StopMuAlgoTypes.h"

namespace larlitecv {

  class Step3D {
    // step 3d with connection to closest 2D plane pixel with charge
  public:
    Step3D() {
      next = NULL;
      prev = NULL;
    };
    virtual ~Step3D() {};
    Point3D_t pos; //< 3d position
    std::vector<Point2D_t> closesthits; //< closest 2d plane pixels (can be different times)
    std::vector<Point2D_t> planepositions; //< position in 2D plane view
    Step3D* next; //< next step
    Step3D* prev; //< prev step
  };
  
  class Hit2D : public std::array<int,2> {
    // represents location in image with pixel of interest
    // stores a distance (for 2D distance from start location)
  public:
  
    Hit2D() : distance(0.0) {};
    virtual ~Hit2D() {};
    double distance;
    inline bool operator<( const Hit2D& rhs ) const {
      if ( distance<rhs.distance ) return true;
      return false;
    };
  };
  
  class Hit2DList : public std::vector<Hit2D> {
  public:
    Hit2DList(size_t init_size=0) {
      if ( init_size!=0 )
	resize(init_size);
    };
    virtual ~Hit2DList() {};

    void sort() {
      std::sort(begin(),end()); // will use hit2d operator<
    };
    void add( const Hit2D& ahit ) {
      push_back( ahit );
      sort();
    };
    void emplace( Hit2D&& ahit ) {
      emplace_back( ahit );
      sort();
    };
  };

  class StopMuTracker {

    StopMuTracker() {};

  public:
      
    StopMuTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v );
    virtual ~StopMuTracker() {};

    void trackStopMu(const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
		     const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d );
    
    std::vector<larcv::Image2D> skel_v; // skeleton images
    std::vector< dbscan::dbPoints > imghits;
    std::vector< dbscan::dbscanOutput > clusters;


    Hit2DList* hitlists[3];
    int current_hit[3];
    Step3D* trackstart;
    Step3D* current;

  };

}


#endif
