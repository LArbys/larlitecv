#include "TaggerCROIAlgo.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "DataFormat/Pixel2D.h"

// larlitecv
#include "ThruMu/ConfigBoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerConfig.h"

#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/EndPointFilter.h"
#include "StopMu/StopMuFilterSpacePoints.h"
#include "StopMu/StopMuCluster.h"

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
    std::cout << "   Top: "        << nsides[0] << std::endl;
    std::cout << "   Bottom: "     << nsides[1] << std::endl;
    std::cout << "   Upstream: "   << nsides[2] << std::endl;
    std::cout << "   Downstream: " << nsides[3] << std::endl;

    // run flash tagger
    anode_flash_tagger.flashMatchTrackEnds(   input.opflashes_v, input.img_v, input.badch_v, output.anode_spacepoint_v );
    cathode_flash_tagger.flashMatchTrackEnds( input.opflashes_v, input.img_v, input.badch_v, output.cathode_spacepoint_v );
    imgends_flash_tagger.findImageTrackEnds( input.img_v, input.badch_v, output.imgends_spacepoint_v );
    int totalflashes = (int)output.anode_spacepoint_v.size() + (int)output.cathode_spacepoint_v.size() + (int)output.imgends_spacepoint_v.size();
    std::cout << " Flash Tagger End Points: " << totalflashes << std::endl;
    std::cout << "  Anode: "      << output.anode_spacepoint_v.size() << std::endl;
    std::cout << "  Cathode: "    << output.cathode_spacepoint_v.size() << std::endl;
    std::cout << "  Image Ends: " << output.imgends_spacepoint_v.size() << std::endl;

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

  StopMuPayload TaggerCROIAlgo::runStopMu( const InputPayload& input, const ThruMuPayload& thrumu ) {
  	StopMuPayload output;

  	float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;

  	// Algos
  	larlitecv::StopMuFilterSpacePoints stopmu_filterpts( m_config.stopmu_filterpts_cfg );
  	larlitecv::StopMuCluster           stopmu_cluster( m_config.stopmu_cluster_cfg );

  	// We strip the unused pixel locations into vector of pixels.
  	std::vector< std::vector< const larcv::Pixel2D* > > unused_spacepoint_v;

  	const larcv::ImageMeta& meta = input.img_v.front().meta();

  	for ( auto const& pt : thrumu.unused_spacepoint_v ) {

  		Double_t xyz[3];
  		for (int i=0; i<3; i++)
  			xyz[i] = pt.pos()[i];

  		float tick = xyz[0]/cm_per_tick + 3200.0;
  		if ( tick<=meta.min_y() ) tick = meta.min_y()+1;
  		if ( tick>=meta.max_y() ) tick = meta.max_y()-1;
  		int row = meta.row( tick );
  		std::vector<float> wid( input.img_v.size(),0.0);
  		for (size_t p=0; p<input.img_v.size(); p++) {
        float wire = larutil::Geometry::GetME()->WireCoordinate(xyz,p);
        if ( wire<0 ) wire = 0;
        if ( wire>=3456 ) wire = 3455;
        int col = input.img_v.at(p).meta().col( wire );
        larcv::Pixel2D pix( col, row );
  		  output.stopmu_pixel_endpt_v.Emplace( p, std::move(pix) );
  		}

  	}

  	std::vector< larcv::EventPixel2D* > unused_spacepoints;
  	unused_spacepoints.push_back( &(output.stopmu_pixel_endpt_v) );

  	output.stopmu_candidate_endpt_v = stopmu_filterpts.filterSpacePoints( unused_spacepoints, thrumu.tagged_v, input.badch_v );
  	std::cout << "  Number of candidate stop-mu start points: " << output.stopmu_candidate_endpt_v.size() << std::endl;

  	output.stopmu_trackcluster_v = stopmu_cluster.findStopMuTracks( input.img_v, input.gapch_v, thrumu.tagged_v, output.stopmu_candidate_endpt_v );
  	std::cout << "  Number of candidate StopMu tracks: " << output.stopmu_trackcluster_v.size() << std::endl;

  	return output;

  }

  CROIPayload TaggerCROIAlgo::runCROISelection( const InputPayload& input, const ThruMuPayload& thrumu, const StopMuPayload& stopmu ) {
  	
  }
}