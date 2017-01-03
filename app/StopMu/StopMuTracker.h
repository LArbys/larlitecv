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
#include <utility>

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
    Step3D( const std::vector<float>& _pos, const std::vector<float>& _dir, 
	    const std::vector<Point2D_t>& _closesthits, const std::vector<Point2D_t>& _planepositions, 
	    Step3D& _prev ) : pos(_pos), dir(_dir), closesthits(_closesthits), planepositions(_planepositions)
    {
      next = NULL;
      prev = &_prev;
      _prev.next = this;
    };
    virtual ~Step3D() {};
    
    std::vector<float> pos; //< 3d position
    std::vector<float> dir;
    std::vector<Point2D_t> closesthits; //< closest 2d plane pixels (can be different times)
    std::vector<Point2D_t> planepositions; //< position in 2D plane view
    Step3D* next; //< next step
    Step3D* prev; //< prev step

    Step3D& GetNext() const { 
      if ( next==NULL )
	throw std::runtime_error("Step3D::Next() - next step is NULL");
      return *next; 
    };
    Step3D& GetPrev() const { 
      if ( prev==NULL )
	throw std::runtime_error("Step3D::Prev() - prev step is NULL");
      return *prev; 
    };
    bool isEnd() const { 
      if ( next==NULL )
	return true;
      return false;
    };
    bool isStart() const {
      if ( prev==NULL )
	return true;
      return false;
    };

    void removeNext() {
      if ( next )
	delete next;
      next = NULL;
    };

  };
  
  class Hit2D : public std::array<int,2> {
    // represents location in image with pixel of interest
    // stores a distance (for 2D distance from start location)
  public:
  
    Hit2D() 
      : distance(0.0), marked(0) {};
    virtual ~Hit2D() {};
    double distance;
    int marked;
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
      _verbose_ = false;
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
    // closest hit in the list
    void closestHits( std::vector<double>& test_pos, const larcv::ImageMeta& meta, const float cm_per_tick, const float cm_per_wire, 
		      std::vector< std::pair<int,double> >& hitlist, const int max_nhits=-1, const int ignore_marked=0 ) const;
    void markUpTo( const size_t index ) {
      for (size_t i=0; i<size(); i++) {
	if ( i<index ) at(i).marked = 1;
	else break;
      }
    };

    bool _verbose_;
  };

  class FitDataHack {
    // an admitted hack to pass data to the minimizer.
    // if this ever works, need a better way.

    FitDataHack() 
      : data_hitlist(NULL), data_anchor(NULL), data_prev_dir(NULL), data_meta(NULL), data_images(NULL), dist_weight(1.0), bend_weight(1.0), step_size_cm(1.0), charge_weight(1.0)
      {};
    virtual ~FitDataHack() {};

    static FitDataHack* _global_instance;

  public:

    static FitDataHack* getMe() { 
      if ( _global_instance==NULL )
	_global_instance = new FitDataHack;
      return _global_instance; 
    };
    
    const std::vector<Hit2DList>* data_hitlist;
    const std::vector<float>* data_anchor;
    const std::vector<float>* data_prev_dir;
    const larcv::ImageMeta*   data_meta;
    const std::vector<larcv::Image2D>* data_images;
    float dist_weight;
    float bend_weight;
    float step_size_cm;
    float charge_weight;

  };

  class StopMuTracker {

    StopMuTracker() { 
      m_verbosity = 0;
    };

  public:
      
    StopMuTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, int verbosity=0 );
    virtual ~StopMuTracker() {};

    void trackStopMu(const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
		     const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d, Step3D& start_step );
    void makeProposedPos( const std::vector<float>& currentpos, const std::vector<float>& currentdir, std::vector<float>& proposedpos, const float stepsize );    
    static void imagePositions( const std::vector<float>& currentpos, int& tick, std::vector<int>& wid );
    void getClosestHitsInPlane( const int clusterid, const std::vector<int>& test_pos,
				const ::dbscan::dbPoints& src_data, const ::dbscan::dbscanOutput& cluster_info, const larcv::ImageMeta& meta,
				std::vector< std::pair<int,double> >& hitlist );
    void getClosestHitsInList( const std::vector<int>& test_pos, const Hit2DList& src_list, const larcv::ImageMeta& meta, 
			       std::vector< std::pair<int,double> >& hitlist );
    std::vector<larcv::Image2D> fillSortedHit2Dlist( const larcv::ImageMeta& meta, 
						     const std::vector< std::vector<int> >& start2d,  const std::vector< std::vector<float> >& start_dir2d,
						     std::vector<Hit2DList>& hitlists, std::vector<int>& clusterid );
    bool findNewDirection( Step3D& current_step, const larcv::ImageMeta& meta,
			   const std::vector< std::vector< std::pair<int,double> > >& closest_pixels, const std::vector<Hit2DList>& sorted_hits,
			   Step3D& proposed_step );


    void stopMuString(const std::vector<larcv::Image2D>& img_v, const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
		      const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d, Step3D& start_step );
    void setVerbosity( int v ) { m_verbosity = v; };
    
    std::vector<larcv::Image2D> skel_v; // skeleton images
    std::vector< dbscan::dbPoints > m_imghits;
    std::vector< dbscan::dbscanOutput > m_clusters;
    
    
    int current_hit[3];
    Step3D* trackstart;
    Step3D* current;
    int m_verbosity;

  protected:
    float _norm( std::vector<float>& vec );
    void _wire2pixel( const int tick, const std::vector<int>& wid, const larcv::ImageMeta& meta, std::vector<int>& pixel_col, int& pixel_row );
  };

}


#endif
