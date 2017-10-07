#include "TaggerFlashMatchAlgoConfig.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larlitecv {

  TaggerFlashMatchAlgoConfig::TaggerFlashMatchAlgoConfig() {
    verbosity = 0;
    qcluster_stepsize = 0.3;
    MeV_per_cm = 2.3;
    fudge_factor = 33333.0;
    fudge_factor_cosmic = 16666.5;
    beam_tick_range.resize(2);
    beam_tick_range[0] = 150;
    beam_tick_range[1] = 400;
    us_per_tick = 0.015625;
    pmtflash_thresh = 3;
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
    gain_correction.resize(32,1.0);
    totpe_sigma_cut = 4.0;
    use_gaus2d = false;
  }

  TaggerFlashMatchAlgoConfig TaggerFlashMatchAlgoConfig::FromPSet( const larcv::PSet& pset ) {
    TaggerFlashMatchAlgoConfig cfg;
    cfg.verbosity         = pset.get<int>("Verbosity");
    cfg.qcluster_stepsize = pset.get<float>("QClusterStepSize");
    cfg.MeV_per_cm        = pset.get<float>("MeV_per_cm");
    cfg.fudge_factor      = pset.get<float>("FudgeFactor");
    cfg.fudge_factor_cosmic = pset.get<float>("CosmicDiscFudgeFactor");
    cfg.pmtflash_thresh   = pset.get<float>("PMTFlashThreshold");
    cfg.flashpe_thresh    = pset.get<float>("FlashPEThreshold");
    cfg.beam_tick_range   = pset.get< std::vector<int> >("BeamTickRange");
    cfg.FVCutX            = pset.get< std::vector<float> >("FVCutX");
    cfg.FVCutY            = pset.get< std::vector<float> >("FVCutY");
    cfg.FVCutZ            = pset.get< std::vector<float> >("FVCutZ");
    cfg.flashmatch_chi2_cut = pset.get<float>("FlashMatchChi2Cut");
    cfg.totpe_sigma_cut   = pset.get<float>("TotalPESigmaCut");
    cfg.use_gaus2d        = pset.get<bool>("UseGaus2D");
    cfg.bbox_pad          = pset.get<float>("BBoxPadcm");
    cfg.gain_correction   = pset.get< std::vector<float> >("GainCorrection");
    cfg.genflashmatch_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( pset.get< larcv::PSet >("GeneralFlashMatchAlgo") );
    return cfg;
  }
}
