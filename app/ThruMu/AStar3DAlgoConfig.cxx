#include "AStar3DAlgoConfig.h"

namespace larlitecv {

  AStar3DAlgoConfig AStar3DAlgoConfig::MakeFromPSet( const larcv::PSet& pset ) {
    AStar3DAlgoConfig cfg;

    cfg.astar_threshold        = pset.get< std::vector<float> >( "PixelThresholds" );
    cfg.astar_neighborhood     = pset.get< std::vector<int> >( "NeighborhoodSize" );
    cfg.astar_start_padding    = pset.get< int >( "StartPadding" );
    cfg.astar_end_padding      = pset.get< int >( "EndPadding" );
    cfg.lattice_padding        = pset.get< int >( "LatticePadding" );
    cfg.accept_badch_nodes     = pset.get< bool >( "AcceptBadChannelNodes" );
    cfg.min_nplanes_w_hitpixel = pset.get< int >( "MinNumPlanesWithHitPixel" );
    cfg.restrict_path          = pset.get< bool >( "RestrictPath", false );
    cfg.verbosity              = pset.get< int >( "Verbosity" );
    if ( cfg.restrict_path ) {
      cfg.path_restriction_radius = pset.get<float>("PathRestrictionRadius");
    }
    else
      cfg.path_restriction_radius = pset.get<float>("PathRestrictionRadius",0.0);

    return cfg;
  }

  
}
