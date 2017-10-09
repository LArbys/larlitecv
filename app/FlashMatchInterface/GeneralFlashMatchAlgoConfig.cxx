#include "GeneralFlashMatchAlgoConfig.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larlitecv {

  const std::string GeneralFlashMatchAlgoConfig::m_flashman_default = "FlashMatchManager: {AllowReuseFlash: true CustomAlgo: [0,0] FlashFilterAlgo: \"\" HypothesisAlgo: \"PhotonLibHypothesis\" MatchAlgo: \"QLLMatch\" ProhibitAlgo: \"\" StoreFullResult: true TPCFilterAlgo: \"TimeCompatMatch\" Verbosity: 02 }";

  GeneralFlashMatchAlgoConfig::GeneralFlashMatchAlgoConfig():
    m_flashmatch_config(fcllite::PSet("FlashMatchManager",m_flashman_default)) {
    chi2_anode_cathode_cut = 100.0;
    chi2_yz_flash_cut      = 20.0; // Arbitrarily chosen for now.
    verbosity = 0;
    qcluster_stepsize = 0.3;
    MeV_per_cm = 2.3;
    fudge_factor = 66666.0; // This has been increased from 33333.0 (increased by a factor of 2).
    fudge_factor_cosmic = 16666.5;
    beam_tick_range.resize(2);
    beam_tick_range[0] = 150;
    beam_tick_range[1] = 400;
    us_per_tick = 0.015625;
    flashpe_thresh = 10.0;
    bbox_pad = 20.0;
    flashmatch_chi2_cut = 10.0;
    FVCutX.resize(2);
    FVCutY.resize(2);
    FVCutZ.resize(2);
    FVCutX[0] = 0.0;
    FVCutX[1] = 275.0;
    FVCutY[0] = -118.0;
    FVCutY[1] = 118.0;
    FVCutZ[0] = 0.0;
    FVCutZ[1] = 1037.0;
    totpe_sigma_cut = 4.0;
    gain_correction.resize(32,1.0);
    pmtflash_thresh = 3;
    use_gaus2d = false;
  }

  GeneralFlashMatchAlgoConfig GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( const larcv::PSet& pset ) {

    std::cout << "Making the 'GeneralFlashMatchAlgoConfig' file into a PSet." << std::endl;

    GeneralFlashMatchAlgoConfig cfg;
    cfg.verbosity              = pset.get<int>("Verbosity");
    cfg.fudge_factor           = pset.get<float>("FudgeFactor");
    cfg.fudge_factor_cosmic = pset.get<float>("CosmicDiscFudgeFactor");
    cfg.pmtflash_thresh   = pset.get<float>("PMTFlashThreshold");
    cfg.flashpe_thresh    = pset.get<float>("FlashPEThreshold");
    cfg.beam_tick_range   = pset.get< std::vector<int> >("BeamTickRange");
    cfg.flashmatch_chi2_cut = pset.get<float>("FlashMatchChi2Cut");
    cfg.totpe_sigma_cut   = pset.get<float>("TotalPESigmaCut");
    cfg.use_gaus2d        = pset.get<bool>("UseGaus2D");
    //std::cout << "converting pset from larcv::PSet to fcllite::PSet" << std::endl;
    fcllite::PSet vox( "GenFlashMatchAlgoManager", pset.data_string() );
    cfg.m_flashmatch_config    = vox.get<fcllite::PSet>("GeneralFlashMatchAlgo");
    //std::cout << cfg.m_flashmatch_config.dump() << std::endl;
    return cfg;
  }
}

