#include "TaggerCROIAlgo.h"

// larlite

// larcv

// larlitecv
#include "ThruMu/ConfigBoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerConfig.h"

#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/EndPointFilter.h"

namespace larlitecv {

  ThruMuPayload TaggerCROIAlgo::runThruMu( const InputPayload& input ) {

    ThruMuPayload output;

    // configure different stages

    // side tagger
    larlitecv::BoundaryMuonTaggerAlgo sidetagger;
    sidetagger.configure( m_config.sidetagger_cfg );
    sidetagger.printConfiguration();

    // flash tagger
    larlitecv::FlashMuonTaggerAlgo anode_flash_tagger(   larlitecv::FlashMuonTaggerAlgo::kAnode );
    larlitecv::FlashMuonTaggerAlgo cathode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kCathode );
    larlitecv::FlashMuonTaggerAlgo imgends_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kOutOfImage );

    anode_flash_tagger.configure(   m_config.flashtagger_cfg );
    cathode_flash_tagger.configure( m_config.flashtagger_cfg );
    imgends_flash_tagger.configure( m_config.flashtagger_cfg );

    // end point filter
    larlitecv::EndPointFilter endptfilter;

	  // RUN THE THRUMU ALGOS

    // run side tagger
    sidetagger.searchforboundarypixels3D( input.img_v, input.badch_v, output.side_spacepoint_v, output.boundarypixel_image_v, output.realspacehit_image_v );
    int nsides[4] = {0};
    for ( auto const& sp : output.side_spacepoint_v ) {
    	nsides[ sp.at(0).type ]++;
    }
    std::cout << " Side Tagger End Points: " << output.side_spacepoint_v.size() << std::endl;
    std::cout << "  Top: "        << nsides[0] << std::endl;
    std::cout << "  Bottom: "     << nsides[1] << std::endl;
    std::cout << "  Upstream: "   << nsides[2] << std::endl;
    std::cout << "  Downstream: " << nsides[3] << std::endl;

    // run flash tagger
    anode_flash_tagger.flashMatchTrackEnds(   input.opflashes_v, input.img_v, input.badch_v, output.anode_spacepoint_v );
    cathode_flash_tagger.flashMatchTrackEnds( input.opflashes_v, input.img_v, input.badch_v, output.cathode_spacepoint_v );
    imgends_flash_tagger.flashMatchTrackEnds( input.opflashes_v, input.img_v, input.badch_v, output.imgends_spacepoint_v );

    // we collect pointers to all the end points
    std::vector< const larlitecv::BoundarySpacePoint* > all_endpoints;

    // gather endpoints from space points
    for (int isp=0; isp<(int)output.side_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.side_spacepoint_v.at( isp ));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.anode_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.anode_spacepoint_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.cathode_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.cathode_spacepoint_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.imgends_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.imgends_spacepoint_v.at(isp));
      all_endpoints.push_back( pts );
    }

    // filter spacepoints
    std::vector<int> endpoint_passes( all_endpoints.size(), 1 );
    endptfilter.removeBoundaryAndFlashDuplicates( all_endpoints, input.img_v, input.gapch_v, endpoint_passes );

    // remove the filtered end points
    std::vector< const larlitecv::BoundarySpacePoint* > filtered_endpoints;
    for ( size_t idx=0; idx<endpoint_passes.size(); idx++ ) {
      if ( endpoint_passes.at(idx)==1 ) {
        filtered_endpoints.push_back( all_endpoints.at(idx) );
      }
    }

    // make track clusters
    std::vector<int> used_filtered_endpoints( filtered_endpoints.size(), 0 );
    sidetagger.makeTrackClusters3D( input.img_v, input.gapch_v, filtered_endpoints, output.trackcluster3d_v, output.tagged_v, used_filtered_endpoints );

    // collect unused endpoints
    for ( size_t isp=0; isp<filtered_endpoints.size(); isp++ ) {
    	if ( used_filtered_endpoints.at(isp)==1 ) 
    		output.used_spacepoint_v.push_back( *(filtered_endpoints.at(isp)) );
    	else
    		output.unused_spacepoint_v.push_back( *(filtered_endpoints.at(isp)) );
    }

    // return stage output
    return output;
  }
}