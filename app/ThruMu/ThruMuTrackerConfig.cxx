#include "ThruMuTrackerConfig.h"

// stdlib
#include <sstream>

namespace larlitecv {

  ThruMuTrackerConfig::ThruMuTrackerConfig() {
    setDefaults();
  }

  void ThruMuTrackerConfig::setDefaults() {

    verbosity = 0;
    num_passes = 1;
    tag_neighborhood.resize(3,5);
    pixel_threshold.resize(3,10);
    compression_mode = 2;
    downsampling_factor = 4;
    thrumu_flashmatch = true;

    // setup example pass config
    ThruMuPassConfig passcfg;
    passcfg.run_linear_tagger = true;
    passcfg.run_astar_tagger  = true;
    passcfg.run_radial_filter = false;
    passcfg.run_foxtrot_extender = false;
    passcfg.min_point_separation = 0.0;

    // pars for controlling track acceptance and astar running
    passcfg.linear3d_min_tracksize = 15;
    passcfg.linear3d_min_goodfraction =  0.9;
    passcfg.linear3d_min_majoritychargefraction = 0.8;
    passcfg.astar3d_min_goodfrac = 0.2;
    passcfg.astar3d_min_majfrac  = 0.2;

    // astar configuration
    passcfg.astar3d_cfg.astar_threshold.resize(3);
    passcfg.astar3d_cfg.astar_threshold[0] =  50.0;
    passcfg.astar3d_cfg.astar_threshold[1] =  50.0;
    passcfg.astar3d_cfg.astar_threshold[2] = 100.0;
    passcfg.astar3d_cfg.astar_neighborhood.resize(3,5);
    passcfg.astar3d_cfg.astar_start_padding = 3;
    passcfg.astar3d_cfg.astar_end_padding = 3;
    passcfg.astar3d_cfg.lattice_padding = 10;
    passcfg.astar3d_cfg.accept_badch_nodes = true;
    passcfg.astar3d_cfg.min_nplanes_w_hitpixel = 3;
    passcfg.astar3d_cfg.min_nplanes_w_charge = 3;
    passcfg.astar3d_cfg.restrict_path = true;
    passcfg.astar3d_cfg.verbosity = 0;
    passcfg.astar3d_cfg.path_restriction_radius = 30.0;

    // linear 3d tagger
    passcfg.linear3d_cfg.trigger_tpc_tick = 3200.0;
    passcfg.linear3d_cfg.min_ADC_value = 10.0;
    passcfg.linear3d_cfg.step_size = 1.5;
    passcfg.linear3d_cfg.neighborhood_square = 5;
    passcfg.linear3d_cfg.neighborhood_posttick = 5;

    // foxtrot tracker
    passcfg.foxextend_cfg.step_size = 4.0;
    passcfg.foxextend_cfg.num_step_attempts = 3;
    passcfg.foxextend_cfg.pixel_thresholds.resize(3,10.0);
    passcfg.foxextend_cfg.min_hit_width = 1;
    passcfg.foxextend_cfg.hit_neighborhood = 2;
    passcfg.foxextend_cfg.segment_frac_w_charge = 0.5;
    passcfg.foxextend_cfg.radius_reduction_factor = 0.5;
    passcfg.foxextend_cfg.min_cosine = 0.5;
    passcfg.foxextend_cfg.max_steps = 100;
    passcfg.foxextend_cfg.verbosity = 0;

    pass_configs.emplace_back( std::move(passcfg) );

  }

  ThruMuTrackerConfig ThruMuTrackerConfig::MakeFromPSet( const larcv::PSet& pset ) {
    ThruMuTrackerConfig cfg;
    cfg.verbosity  = pset.get<int>("Verbosity");
    cfg.num_passes = pset.get<int>("NumPasses");
    cfg.tag_neighborhood = pset.get< std::vector<int> >("TaggingNeighborhood");
    cfg.pixel_threshold  = pset.get< std::vector<float> >("PixelThresholds");
    cfg.downsampling_factor = pset.get<int>("DownsampleFactor");
    cfg.compression_mode = pset.get<int>("CompressionMode");
    cfg.thrumu_flashmatch = pset.get<bool>("ThruMuFlashMatch");
    cfg.pass_configs.clear();
    for ( int ipass=0; ipass<cfg.num_passes; ipass++ ) {
      std::stringstream ss;
      ss << "ThruMuPassConfig" << ipass;
      larcv::PSet pass_pset = pset.get< larcv::PSet >( ss.str() );
      ThruMuTrackerConfig::ThruMuPassConfig passcfg;
      passcfg.run_linear_tagger                   = pass_pset.get<bool>("RunLinearTagger");
      passcfg.run_astar_tagger                    = pass_pset.get<bool>("RunAStarTagger");
      passcfg.run_radial_filter                   = pass_pset.get<bool>("RunRadialFilter");
      passcfg.run_foxtrot_extender                = pass_pset.get<bool>("RunFoxTrotExtender");
      passcfg.min_point_separation                = pass_pset.get<float>("MinPointSeparationcm");
      passcfg.linear3d_min_tracksize              = pass_pset.get<int>("Linear3DMinTrackSize");
      passcfg.linear3d_min_goodfraction           = pass_pset.get<float>("Linear3DMinGoodFraction");
      passcfg.linear3d_min_majoritychargefraction = pass_pset.get<float>("Linear3DMinMajorityChargeFraction");
      passcfg.astar3d_min_goodfrac                = pass_pset.get<float>("AStar3DMinGoodFraction");
      passcfg.astar3d_min_majfrac                 = pass_pset.get<float>("AStar3DMinMajorityChargeFraction");
      passcfg.astar3d_cfg   = larlitecv::AStar3DAlgoConfig::MakeFromPSet( pass_pset.get<larcv::PSet>( "AStarConfig" ) );
      passcfg.linear3d_cfg  = larlitecv::Linear3DChargeTaggerConfig::makeFromPSet( pass_pset.get<larcv::PSet>("Linear3DConfig" ) );
      passcfg.radial_cfg    = larlitecv::RadialEndpointFilterConfig::makeFromPSet( pass_pset.get<larcv::PSet>("RadialFilterConfig") );
      if ( passcfg.run_foxtrot_extender ) {
        passcfg.foxextend_cfg = larlitecv::FoxTrotTrackerAlgoConfig::makeFromPSet( pass_pset.get<larcv::PSet>("FoxTrotConfig") );
      }
      // std::cout << ss.str() << ": runradial=" << passcfg.run_radial_filter
                // << " runlinear=" << passcfg.run_linear_tagger
                // << " runastar=" << passcfg.run_astar_tagger
                // << " runfoxextender=" << passcfg.run_foxtrot_extender
                // << std::endl;
      cfg.pass_configs.emplace_back( std::move(passcfg) );
    }
    return cfg;
  }

}
