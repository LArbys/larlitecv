#ifndef __FLASH_MUON_TAGGER_ALGO__
#define __FLASH_MUON_TAGGER_ALGO__

#include <vector>

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "ANN/ANNAlgo.h"
#include "dbscan/DBSCANAlgo.h"

// larlitecv header
#include "BoundaryEndPt.h"

namespace larlitecv {
  
  class FlashMuonTaggerConfig {
  public:
    FlashMuonTaggerConfig() {};
    virtual ~FlashMuonTaggerConfig() {};
    
    // values per plane
    std::vector<float>  pixel_value_threshold;
    std::vector<int>    clustering_time_neighborhood;
    std::vector<int>    clustering_wire_neighborhood;
    std::vector<int>    clustering_minpoints;
    std::vector<double> clustering_radius;
    std::vector<int>    endpoint_time_neighborhood;
    int                 verbosity;
    float               trigger_tick;
    float               usec_per_tick;
    float               drift_distance;
    float               drift_velocity;
  };
  
  class FlashMuonTaggerAlgo {

  public:
    
    typedef enum { kAnode, kCathode } SearchMode_t;

    FlashMuonTaggerAlgo( SearchMode_t mode_) {
      fSearchMode = mode_;
      fconfigured = false;
    };
    virtual ~FlashMuonTaggerAlgo() {};
    
    void configure( const FlashMuonTaggerConfig& config_ ) { fConfig=config_; fconfigured = true; };
    bool findTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const larcv::Image2D& image,
			std::vector< BoundaryEndPt >& trackendpts, larcv::Image2D& markedimg );
    
    
  protected:
    
    SearchMode_t fSearchMode;
    bool fconfigured;
    FlashMuonTaggerConfig fConfig;
    
    // subroutines
    bool getClusters( const larcv::Image2D& tpc_img, int query_row, int query_col, dbscan::dbscanOutput& cluster_info, dbscan::dbPoints& windowpts, int& containing_cluster );
    
  };
  
  
}

#endif
