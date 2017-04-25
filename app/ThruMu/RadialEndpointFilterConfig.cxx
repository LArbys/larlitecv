#include "RadialEndpointFilterConfig.h"

namespace larlitecv {

  RadialEndpointFilterConfig::RadialEndpointFilterConfig() {
    // default constructor with some sensible defaults
    segment_radius = 10.0; // cm
    segment_min_width = 1; // pixels
    segment_frac_w_charge = 0.8; // fraction
    acceptance_angle = (10.0/180.0)*3.14159; // 5 degrees in radians
    pixel_thresholds.resize(3,10.0);
    min_segments = 0;
    max_segments = 2;
  }

  RadialEndpointFilterConfig::~RadialEndpointFilterConfig() {
    pixel_thresholds.clear();
  }

  RadialEndpointFilterConfig RadialEndpointFilterConfig::makeFromPSet( const larcv::PSet& pset ) {
    RadialEndpointFilterConfig cfg;
    cfg.segment_radius        = pset.get<float>("SegmentRadius_cm");
    cfg.segment_min_width     = pset.get<int>("SegmentMinWidth");
    cfg.segment_frac_w_charge = pset.get<float>("SegmentFractionWithCharge");
    cfg.acceptance_angle      = pset.get<float>("AcceptanceAngle");
    cfg.pixel_thresholds      = pset.get<std::vector<float> >("PixelThresholds");
    cfg.min_segments          = pset.get<int>("MinNumSegments");
    cfg.max_segments          = pset.get<int>("MaxNumSegments");
    return cfg;
  }

}
