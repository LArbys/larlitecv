#ifndef __STOP_MU_START__
#define __STOP_MU_START__

#include <vector>
#include <list>
#include <array>
#include <queue>

#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Pixel2D.h"

#include "StopMuAlgoTypes.h"

namespace larlitecv {
    
  class StopMuStart {

  public:
    
    StopMuStart();
    virtual ~StopMuStart();
    
    // time step tracker functions
    void getStartDirection( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
			    const std::vector<larcv::Pixel2D>& start, const int rneighbor, const int cneighbor, const float fThreshold,
			    std::array<float,3>& start_spacepoint, 
			    std::vector< std::array<float,2> >& start_dir2d, 
			    std::array<float,3>& start_dir3d ); //< for a given starting point get start direction
    void getStartDirectionV( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
			     const std::vector<larcv::Pixel2D>& start, const int rneighbor, const int cneighbor, const float fThreshold,
			     std::vector<float>& start_spacepoint, 
			     std::vector< std::vector<float> >& start_dir2d, std::vector<float>& start_dir3d ); //< for a given starting point get start direction

    void setVerbose( int v ) { verbose = v; };
  protected:

    int verbose;

  };

}

#endif
