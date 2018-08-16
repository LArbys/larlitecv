#include "TaggerFlashMatchAlgoConfig.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larlitecv {

  TaggerFlashMatchAlgoConfig::TaggerFlashMatchAlgoConfig() {
    verbosity = 0;
    use_version = 1;
    qcluster_stepsize = 0.3;
    MeV_per_cm = 2.3;
    fudge_factor = 33333.0;
    fudge_factor_cosmic = 16666.5;
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
    extend_cosmics = true;
  }

  TaggerFlashMatchAlgoConfig TaggerFlashMatchAlgoConfig::FromPSet( const larcv::PSet& pset ) {
    TaggerFlashMatchAlgoConfig cfg;
    cfg.verbosity         = pset.get<int>("Verbosity");
    cfg.FVCutX            = pset.get< std::vector<float> >("FVCutX");
    cfg.FVCutY            = pset.get< std::vector<float> >("FVCutY");
    cfg.FVCutZ            = pset.get< std::vector<float> >("FVCutZ");
    cfg.use_version       = pset.get< int >("Version");
    cfg.extend_cosmics    = pset.get< bool >("ExtendCosmicTracks");
    cfg.use_fixed_croi    = pset.get< bool >("UseFixedCROI");           //< CROI are made using flash position. No selection based on flash-matching.
    cfg.split_fixed_ycroi = pset.get< bool >("SplitFixedYCROI");        //< CROI are made using flash position. No selection based on flash-matching.
    if ( cfg.use_version<1 || cfg.use_version>2 )
      throw std::runtime_error( "TaggerFlashMatchAlgoConfig::FromPSet : Version must be 1 or 2" );
    //cfg.gain_correction   = pset.get< std::vector<float> >("GainCorrection");
    cfg.genflashmatch_cfg = larlitecv::GeneralFlashMatchAlgoConfig::MakeConfigFromPSet( pset.get< larcv::PSet >("GeneralFlashMatchAlgo") );
    return cfg;
  }
}
