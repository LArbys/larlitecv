#ifndef __GENERAL_FLASH_MATCH_ALGO_CONFIG_H__
#define __GENERAL_FLASH_MATCH_ALGO_CONFIG_H__

#include <vector>

// LArCV
#include "Base/PSet.h"

// LArLite
#include "FhiclLite/PSet.h"

namespace larlitecv {

  class TaggerFlashMatchAlgoConfig {
  public:
    TaggerFlashMatchAlgoConfig();
    virtual ~TaggerFlashMatchAlgoConfig() {};
    
    int verbosity;
    float qcluster_stepsize;
    float MeV_per_cm;
    float fudge_factor;
    float pmtflash_thresh;
    bool use_gaus2d;
    fcllite::PSet m_flashmatch_config;
    vector < float > gain_correction;
    static const std::string m_flashman_default;
    
    static TaggerFlashMatchAlgoConfig MakeTaggerFlashMatchAlgoConfigFromPSet( const larcv::PSet& pset );
  };

}

#endif
