#include "FoxTrotTrackerAlgoConfig.h"

namespace larlitecv {

  FoxTrotTrackerAlgoConfig::FoxTrotTrackerAlgoConfig() {
    // defaults
    step_size = 3.0;
    num_step_attempts = 2;
    pixel_thresholds.resize(3,10.0);
    min_hit_width = 1;
    segment_frac_w_charge = 0.5;
    radius_reduction_factor = 0.5;
    min_cosine = 0.0;
    max_steps = 10000;
    verbosity = 0;
  }

  FoxTrotTrackerAlgoConfig FoxTrotTrackerAlgoConfig::makeFromPSet( const larcv::PSet& pset ) {
    FoxTrotTrackerAlgoConfig cfg;
    cfg.step_size = pset.get<float>("StepSizecm");
    cfg.num_step_attempts = pset.get<int>("NumStepAttempts");
    cfg.pixel_thresholds = pset.get<std::vector<float> >("PixelThresholds");
    cfg.min_hit_width = pset.get<int>("SegmentMinHitWidth");
    cfg.segment_frac_w_charge = pset.get<float>("SegmentFractionWithCharge");
    cfg.radius_reduction_factor = pset.get<float>("StepReductionFactor");
    cfg.min_cosine = pset.get<float>("MinCosine");
    cfg.max_steps = pset.get<int>("MaxSteps");
    cfg.verbosity = pset.get<int>("Verbosity");
    return cfg;
  }

}

