#ifndef __STOPMU_FILTERSPACEPOINTS__
#define __STOPMU_FILTERSPACEPOINTS__

#include <vector>

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2D.h"
#include "DataFormat/EventPixel2D.h"

namespace larlitecv {

  class StopMuFilterSpacePointsConfig {
  public:
    StopMuFilterSpacePointsConfig() {};
    virtual ~StopMuFilterSpacePointsConfig() {};

    float duplicate_radius_cm;
    float pixel_threshold;
    int track_width_pixels;
    int row_tag_neighborhood;
    int col_tag_neighborhood;
  };

  class StopMuFilterSpacePoints {

    StopMuFilterSpacePoints() {}; // empty constructor
    
  public:
    StopMuFilterSpacePoints( const StopMuFilterSpacePointsConfig& config  );
    virtual ~StopMuFilterSpacePoints() {};


    const StopMuFilterSpacePointsConfig* m_config; //< configuration parameters

    std::vector<larcv::Pixel2D> filterSpacePoints( const std::vector< larcv::EventPixel2D* >& spacepoints_list, 
						   const std::vector<larcv::Image2D>& thrumu_pixels,
						   const std::vector<larcv::Image2D>& badch_imgs);

    void removeThroughGoingEndPoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
				      const std::vector<larcv::Image2D>& thrumu_pixels,
				      std::vector< std::vector<const larcv::Pixel2D*> >& passlist );
    void removeThroughGoingEndPointsFromPixVectors( std::vector< std::vector<const larcv::Pixel2D*> >& spacepoints_list,
						    const std::vector<larcv::Image2D>& thrumu_pixels,
						    std::vector< std::vector<const larcv::Pixel2D*> >& passlist );
    bool isEndPtNearThruMuTag( const std::vector<const larcv::Pixel2D*>& pix_v, const std::vector<larcv::Image2D>& thrumu_pixels );
    
    
    void removeDuplicateEndPoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
				   const std::vector< larcv::Image2D >& img_v,
				   std::vector< std::vector<const larcv::Pixel2D*> >& outputlist );
    
    float dwall( std::vector<float>& pos3d );



  };

}

#endif
