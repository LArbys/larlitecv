#include "TaggerFlashMatchAlgoConfig.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larlitecv {

	TaggerFlashMatchAlgoConfig::TaggerFlashMatchAlgoConfig() 
	 : flashmatch_config("TaggerFlashMatchAlgo") {
		verbosity = 0;
		qcluster_stepsize = 0.3;
		MeV_per_cm = 2.3;
		fudge_factor = 1.0;
		beam_tick_range.resize(2);
		beam_tick_range[0] = 150;
		beam_tick_range[1] = 400;
		us_per_tick = 0.015625;
		pmtflash_thresh = 3;
		flashmatch_chi2_cut = 10.0;
		FVCutX.resize(2);
		FVCutY.resize(2);
		FVCutZ.resize(2);				
		FVCutX[0] = 5.0;
		FVCutX[1] = 270.0;
		FVCutY[0] = -113.0;
		FVCutY[1] = 113.0;
		FVCutZ[0] = 5.0;
		FVCutZ[1] = 1032.0;
		gain_correction.resize(32,1.0);
	}

	TaggerFlashMatchAlgoConfig TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( const larcv::PSet& pset ) {
		TaggerFlashMatchAlgoConfig cfg;
		cfg.verbosity         = pset.get<int>("Verbosity");
		cfg.qcluster_stepsize = pset.get<float>("QClusterStepSize");
		cfg.MeV_per_cm        = pset.get<float>("MeV_per_cm");
		cfg.fudge_factor      = pset.get<float>("FudgeFactor");
		cfg.pmtflash_thresh   = pset.get<float>("PMTFlashThreshold");
		cfg.beam_tick_range   = pset.get< std::vector<int> >("BeamTickRange");
		cfg.FVCutX            = pset.get< std::vector<float> >("FVCutX");
		cfg.FVCutY            = pset.get< std::vector<float> >("FVCutY");
		cfg.FVCutZ            = pset.get< std::vector<float> >("FVCutZ");
		cfg.flashmatch_chi2_cut = pset.get<float>("FlashMatchChi2Cut");
		cfg.gain_correction   = pset.get< std::vector<float> >("GainCorrection");
		larcv::PSet flashmatchman_cfg = pset.get<larcv::PSet>("FlashMatchManager");
		cfg.flashmatch_config = fcllite::PSet( "temp", flashmatchman_cfg.data_string() ).get<fcllite::PSet>("FlashMatchManager");
		return cfg;
	}

}