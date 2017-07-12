#include "BoundaryMuonTaggerAlgoConfig.h"

namespace larlitecv {

  void BoundaryMuonTaggerAlgoConfig::setdefaults() {

    save_endpt_images = false;
    hitsearch_uses_badchs = false;
    type_modifier.resize(4,1.0);
    type_modifier[2] = 0.5;
    type_modifier[3] = 0.5;

    neighborhoods.resize(3,10);
    thresholds.resize(3,10.0);
    boundary_cluster_minpixels.resize(3,10);
    boundary_cluster_radius.resize(3,10);
    fKernelRadius = 2;
    verbosity = 0;    
  }

  void BoundaryMuonTaggerAlgoConfig::print() {
    std::cout << "==================================================" << std::endl;
    std::cout << "BoundaryMuonTaggerAlgoConfig" << std::endl;

    std::cout << " Pixel Thresholds: "; for ( size_t i=0; i<thresholds.size(); i++ ) std::cout << thresholds[i] << " "; std::cout << std::endl;
    std::cout << " Search Neighborhood: "; for ( size_t i=0; i<neighborhoods.size(); i++ ) std::cout << neighborhoods[i] << " "; std::cout << std::endl;
    std::cout << " Hit Search Uses BadChs: " << hitsearch_uses_badchs << std::endl;
    std::cout << " Boundary Cluster Min Pixels: "; for ( size_t i=0; i<boundary_cluster_minpixels.size(); i++ ) std::cout << boundary_cluster_minpixels[i] << " "; std::cout << std::endl;
    std::cout << " Boundary Cluster Radius: "; for ( size_t i=0; i<boundary_cluster_radius.size(); i++ ) std::cout << boundary_cluster_radius[i] << " "; std::cout << std::endl;
    std::cout << " Type Modifer: "; for ( size_t i=0; i<type_modifier.size(); i++ ) std::cout << type_modifier[i] << " "; std::cout << std::endl;
    std::cout << " Kernel Radius: " << fKernelRadius << std::endl;
    std::cout << " Verbosity: " << verbosity << std::endl;
    std::cout << " Save Endpoint Images: " << save_endpt_images << std::endl;
    std::cout << "==================================================" << std::endl;
  }

  BoundaryMuonTaggerAlgoConfig MakeBoundaryMuonTaggerAlgoConfigFromPSet( const larcv::PSet& pset ) {

    BoundaryMuonTaggerAlgoConfig config;
    config.neighborhoods              = pset.get< std::vector<int> >("Neighborhoods");
    config.thresholds                 = pset.get< std::vector<float> >( "Thresholds" );
    config.boundary_cluster_minpixels = pset.get< std::vector<int> >( "BoundaryClusterMinPixels" );
    config.boundary_cluster_radius    = pset.get< std::vector<float> >( "BoundaryClusterRadius" );
    config.save_endpt_images          = pset.get<bool>("SaveMatchImages",false);
    config.hitsearch_uses_badchs      = pset.get<bool>("UseBadChannels",true);
    config.fKernelRadius              = pset.get<int>("KernelRadius",2);
    config.verbosity                  = pset.get<int>("Verbosity",0);

    return config;
  }

}
