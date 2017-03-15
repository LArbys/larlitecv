#include "TaggerFlashMatchAlgoConfig.h"

namespace larlitecv {

	TaggerFlashMatchAlgoConfig::TaggerFlashMatchAlgoConfig() {
		qcluster_stepsize = 0.3;
		MeV_per_cm = 2.3;
		fudge_factor = 1.0;
		beam_tick_range.resize(2);
		beam_tick_range[0] = 150;
		beam_tick_range[1] = 400;
		FVCutX.resize(2);
		FVCutX[0] = 5.0;
		FVCutX[1] = 270.0;
		FVCutCurveYtop.
	}

	TaggerFlashMatchAlgoConfig TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( const larcv::PSet& pset ) {
		TaggerFlashMatchAlgoConfig cfg;
		cfg.qcluster_stepsize = pset.get<int>("Verbosity");
		cfg.MeV_per_cm        = pset.get<float>("MeV_per_cm");
		cfg.fudge_factor      = pset.get<float>("FudgeFactor");
		cfg.beam_tick_range   = pset.get< std::vector<int> >("BeamTickRange");
		cfg.flashmatch_config = larlite::PSet( "temp", pset.get<std::string>("FlashMatchManager").data_string() ).get<larlite::PSet>("FlashMatchManager");
		return cfg;
	}

}