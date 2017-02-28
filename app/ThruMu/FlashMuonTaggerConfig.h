#ifndef __FLASHMUONTAGGERCONFIG__
#define __FLASHMUONTAGGERCONFIG__

#include <vector>

// larcv
#include "Base/PSet.h"

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
    int                 search_row_radius; ///< neighborhood to gather charge for clustering
    float               flash_zrange_extension;
    int                 flash_pixelcluster_minsize; //< min. number of pixels required to make a flash-tagged cluster
    float               max_triarea;
    float               max_triarea_tight;
    float               cathode_drift_tick_correction;

    void setdefaults();
    
  };
  
  FlashMuonTaggerConfig MakeFlashMuonTaggerConfigFromPSet( const larcv::PSet& flashtagger_pset );

}

#endif
