/* Function to compare a flash of light with a track using the flash matching infrastructure. 

   This is meant to introduce flash matching into the tagger's operations. */
#ifndef __GENERAL_FLASHMATCH_ALGO_H__
#define __GENERAL_FLASHMATCH_ALGO_H__

#include <vector>
#include <cmath>

// larlite 
#include "OpT0Finder/Algorithms/QLLMatch.h"

// larcv
#include "Base/PSet.h"

// larlitecv
#include "ContainedROI/TaggerFlashMatchAlgo.h"
#include "ContainedROI/TaggerFlashMatchAlgoConfig.h"
#include "TaggerCROI/TaggerCROITypes.h"
#include "TaggerTypes/BMTrackCluster3D.h"

namespace larlitecv {

  // define the class
  class GeneralFlashMatchAlgo {
  public:
    GeneralFlashMatchAlgo() : m_config() {};
    virtual ~GeneralFlashMatchAlgo();

    ThruMuFlashMatchAlgo m_config;

    // This might need more information based on the needs of converting a flash to an opflash.
    std::vector< int > non_anodecathode_flash_idx(const std::vector< larlite::event_opflash* >& opflash_v, std::vector < int > anode_flash_idx_v, std::vector < int > cathode_flash_idx_v);
    larlite::track     generate_tracks_between_passes(std::vector< BMTrackCluster3D > trackcluster3d_v);
    
    
    

  };


} 

#endif

