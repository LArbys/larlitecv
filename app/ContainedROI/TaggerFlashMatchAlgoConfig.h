#ifndef __TAGGER_FLASH_MATCH_ALGO_CONFIG_H__
#define __TAGGER_FLASH_MATCH_ALGO_CONFIG_H__

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
    float fudge_factor_cosmic;
    float us_per_tick;
    float pmtflash_thresh;
    float flashmatch_chi2_cut;
    float flashpe_thresh;
    float totpe_sigma_cut;
    float bbox_pad;
    bool use_gaus2d;
    std::vector<int> beam_tick_range;
    std::vector<float> FVCutX;
    std::vector<float> FVCutY;
    std::vector<float> FVCutZ;
    std::vector<float> gain_correction;
    fcllite::PSet m_flashmatch_config;
    static const std::string m_flashman_default;
    
    static TaggerFlashMatchAlgoConfig MakeTaggerFlashMatchAlgoConfigFromPSet( const larcv::PSet& pset );
  };

}

#endif
