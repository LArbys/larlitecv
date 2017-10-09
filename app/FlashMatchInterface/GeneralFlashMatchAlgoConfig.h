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

    float chi2_anode_cathode_cut;
    float chi2_yz_flash_cut;
    int   verbosity;
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
    bool  use_gaus2d;
    std::vector<int> beam_tick_range;
    std::vector<float> FVCutX;
    std::vector<float> FVCutY;
    std::vector<float> FVCutZ;
    fcllite::PSet m_flashmatch_config;
    std::vector < float > gain_correction;
    static const std::string m_flashman_default;
    
    static GeneralFlashMatchAlgoConfig MakeConfigFromPSet( const larcv::PSet& pset );
  };

}

#endif
