#ifndef __GENERAL_FLASH_MATCH_ALGO_CONFIG_H__
#define __GENERAL_FLASH_MATCH_ALGO_CONFIG_H__

#include <vector>

// LArCV
#include "Base/PSet.h"

// LArLite
#include "FhiclLite/PSet.h"

namespace larlitecv {

  class GeneralFlashMatchAlgoConfig {
  public:
    GeneralFlashMatchAlgoConfig();
    virtual ~GeneralFlashMatchAlgoConfig() {};
    
    int verbosity;
    float qcluster_stepsize;
    float MeV_per_cm;
    float fudge_factor;
    float pmtflash_thresh;
    bool use_gaus2d;
    fcllite::PSet m_flashmatch_config;
    std::vector < float > gain_correction;
    static const std::string m_flashman_default;
    
    static GeneralFlashMatchAlgoConfig MakeConfigFromPSet( const larcv::PSet& pset );
  };

}

#endif
