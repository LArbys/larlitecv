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
    gain_correction.resize(32,1.0);
    pmtflash_thresh = 3;
    use_gaus2d = false;
  }

  GeneralFlashMatchAlgoConfig GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( const larcv::PSet& pset ) {

    std::cout << "Making the 'GeneralFlashMatchAlgoConfig' file into a PSet." << std::endl;

    GeneralFlashMatchAlgoConfig cfg;
    cfg.chi2_anode_cathode_cut = pset.get<float>("Chi2_Anode_Cathode_Cut");
    cfg.chi2_yz_flash_cut      = pset.get<float>("Chi2_YZ_Flash_Cut");
    cfg.verbosity              = pset.get<int>("Verbosity");
    cfg.qcluster_stepsize      = pset.get<float>("QClusterStepSize");
    cfg.MeV_per_cm             = pset.get<float>("MeV_per_cm");
    cfg.fudge_factor           = pset.get<float>("FudgeFactor");
    cfg.pmtflash_thresh        = pset.get<float>("PMTFlashThreshold");
    cfg.use_gaus2d             = pset.get<bool>("UseGaus2D");
    fcllite::PSet vox( "Manager", pset.data_string() );
    cfg.m_flashmatch_config    = vox.get<fcllite::PSet>("GeneralFlashMatchAlgo");
    return cfg;
  }
}

