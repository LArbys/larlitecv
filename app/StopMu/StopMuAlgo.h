#ifndef __STOP_MU_ALGO__
#define __STOP_MU_ALGO__

#include <vector>
#include <array>

#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2D.h"

#include "StopMuAlgoTypes.h"

namespace larlitecv {

  class HitNeighborhood {
  public:
  HitNeighborhood( int s, int e, int mc, float mv ) 
    : start(s), end(e), maxcol(mc), maxval(mv) {};
    virtual ~HitNeighborhood() {};
    int start;
    int end;
    int maxcol;
    float maxval;
  };
  
  class StopMuAlgo {

  public:

    /*     // for convenience */
    /*     typedef std::array<float,2>  Vec2D_t; */
    /*     typedef std::array<float,3>  Vec3D_t; */
    /*     typedef std::vector< Vec2D_t > PlaneVec2D_t; */
    /*     typedef std::list< PlaneVec2D_t > Vec2DList_t; */
    /*     typedef std::list< Vec3D_t > Vec3DList_t; */
    
    StopMuAlgo();
    virtual ~StopMuAlgo();

    void runTimeStepTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Pixel2D>& start,
			     std::vector< std::vector<larcv::Pixel2D> >& pixellist,  std::list< std::array<float,3> >& spacepoints ); //< Run this
    
    // time step tracker functions
    void getStartDirection( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
			    const std::vector<larcv::Pixel2D>& start, int rneighbor, int cneighbor, float fThreshold,
			    std::array<float,3>& start_spacepoint, 
			    std::vector< std::array<float,2> >& start_dir2d, std::array<float,3>& start_dir3d ); //< for a given starting point get start direction
    void getStartDirectionV( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
			     const std::vector<larcv::Pixel2D>& start, int rneighbor, int cneighbor, float fThreshold,
			     std::vector<float>& start_spacepoint, 
			     std::vector< std::vector<float> >& start_dir2d, std::vector<float>& start_dir3d ); //< for a given starting point get start direction
			    
    void defineNeighborHoodHitSimple( const std::vector<larcv::Image2D >& img_v, const std::vector<larcv::Pixel2D>& start_pix_v, float threshold,
				      std::vector< HitNeighborhood >& hits ); //< for a given pixel, define the hit around it
    bool getProposedStep( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Pixel2D>& current,
			  float threshold, int neighborhood, std::vector<larcv::Pixel2D>& next_pix_v,
			  std::vector< std::array<float,2> >& local_dir2d );
    bool calcNextPos( const std::vector<larcv::Image2D>& img_v, const std::vector< std::array<float,2> >& dir2d, const std::vector< larcv::Pixel2D >& hit_v,
		      std::vector<larcv::Pixel2D>& next_hit_v );
    bool findNeighborhoodMax( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Pixel2D>& pix_v, float threshold, int neighborhood,
			      std::vector<larcv::Pixel2D>& max_pix_v );
    bool getProposedStep( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Pixel2D>& current,
			  float threshold, int neighborhood, std::vector< std::array<float,2> >& local_dir2d );
    bool amGoingBackwards( const Vec3D_t& nextsp, const Vec3DList_t& history_dir3d, const Vec3DList_t& spacepoints, const Vec3D_t& lastdir3d, 
			   const float fbackward_cos_threshold, Vec3D_t& local_dir3d );
    bool findNearestPix( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Pixel2D>& pix_v, const PlaneVec2D_t& dir2d,
			 const float threshold, const int neighborhood,
			 std::vector<larcv::Pixel2D>& max_pix_v );
    
    void aveVec3DList( const Vec3DList_t& history_dir3d, Vec3D_t& ave_dir3d );

    void setVerbose( int v ) { verbose = v; };
  protected:

    int verbose;

  };
}

#endif
