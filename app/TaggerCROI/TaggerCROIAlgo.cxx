#include "TaggerCROIAlgo.h"

#include <sstream>
#include <ctime>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "DataFormat/Pixel2D.h"
#include "DataFormat/ROI.h"

// larlitecv
#include "ThruMu/FlashMuonTaggerAlgoConfig.h"

#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/EndPointFilter.h"
#include "ThruMu/ThruMuTracker.h"
#include "StopMu/StopMuFilterSpacePoints.h"
#include "StopMu/StopMuCluster.h"
#include "StopMu/StopMuFoxTrot.h"
#include "UntaggedClustering/ClusterGroupAlgo.h"
#include "UntaggedClustering/ClusterGroupMatchingAlgo.h"
#include "ContainedROI/TaggerFlashMatchAlgo.h"

namespace larlitecv {

  TaggerCROIAlgo::TaggerCROIAlgo( const TaggerCROIAlgoConfig& config )
   : m_config(config) {
    m_time_tracker.resize( kNumStages, 0.0 );
  }

  ThruMuPayload TaggerCROIAlgo::runThruMu( const InputPayload& input ) {

    ThruMuPayload output;

    // configure different stages of the Thrumu Tagger
    std::clock_t timer;

    // (1) side tagger
    timer = std::clock();
    larlitecv::BoundaryMuonTaggerAlgo sidetagger;
    sidetagger.configure( m_config.sidetagger_cfg );
    sidetagger.printConfiguration();

    // (2) flash tagger
    larlitecv::FlashMuonTaggerAlgo anode_flash_tagger(   larlitecv::FlashMuonTaggerAlgo::kAnode );
    larlitecv::FlashMuonTaggerAlgo cathode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kCathode );
    larlitecv::FlashMuonTaggerAlgo imgends_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kOutOfImage );

    anode_flash_tagger.configure(   m_config.flashtagger_cfg );
    cathode_flash_tagger.configure( m_config.flashtagger_cfg );
    imgends_flash_tagger.configure( m_config.flashtagger_cfg );

    // (3) end point filter
    larlitecv::EndPointFilter endptfilter;

    // (4) thrumu tracker
    larlitecv::ThruMuTracker thrumu_tracker( m_config.thrumu_tracker_cfg );

    m_time_tracker[kThruMuConfig] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;

    // RUN THE THRUMU ALGOS

    // run side tagger
    timer = std::clock();
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
    m_time_tracker[kThruMuBMT] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // run flash tagger
    timer = std::clock();
    anode_flash_tagger.flashMatchTrackEnds(   input.opflashes_v, input.img_v, input.badch_v, output.anode_spacepoint_v );
    cathode_flash_tagger.flashMatchTrackEnds( input.opflashes_v, input.img_v, input.badch_v, output.cathode_spacepoint_v );
    imgends_flash_tagger.findImageTrackEnds( input.img_v, input.badch_v, output.imgends_spacepoint_v );
    int totalflashes = (int)output.anode_spacepoint_v.size() + (int)output.cathode_spacepoint_v.size() + (int)output.imgends_spacepoint_v.size();
    std::cout << " Flash Tagger End Points: " << totalflashes << std::endl;
    std::cout << "  Anode: "      << output.anode_spacepoint_v.size() << std::endl;
    std::cout << "  Cathode: "    << output.cathode_spacepoint_v.size() << std::endl;
    std::cout << "  Image Ends: " << output.imgends_spacepoint_v.size() << std::endl;
    m_time_tracker[kThruMuFlash] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

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
    if ( m_config.run_thrumu_tracker ) {
      timer = std::clock();
      thrumu_tracker.makeTrackClusters3D( input.img_v, input.gapch_v, filtered_endpoints, output.trackcluster3d_v, output.tagged_v, used_filtered_endpoints );
      m_time_tracker[kThruMuTracker]  +=  (std::clock()-timer)/(double)CLOCKS_PER_SEC;
    }

    // collect unused endpoints
    for ( size_t isp=0; isp<filtered_endpoints.size(); isp++ ) {
      if ( used_filtered_endpoints.at(isp)==1 )
        output.used_spacepoint_v.push_back( *(filtered_endpoints.at(isp)) );
      else
        output.unused_spacepoint_v.push_back( *(filtered_endpoints.at(isp)) );
    }

    // copy track and pixels into separate containers.
    for ( auto const& bmtrack : output.trackcluster3d_v ) {
      output.track_v.push_back( bmtrack.makeTrack() );
      std::vector< larcv::Pixel2DCluster > cluster_v;
      for ( auto const& track2d : bmtrack.plane_pixels )
        cluster_v.push_back( track2d );
      output.pixelcluster_v.emplace_back( std::move(cluster_v) );
    }

    std::cout << "End of ThruMu Algo" << std::endl;

    // return stage output
    return output;
  }

  StopMuPayload TaggerCROIAlgo::runStopMu( const InputPayload& input, const ThruMuPayload& thrumu ) {
    StopMuPayload output;

    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;

    // Algos
    larlitecv::StopMuFilterSpacePoints stopmu_filterpts( m_config.stopmu_filterpts_cfg );
    larlitecv::StopMuCluster           stopmu_cluster( m_config.stopmu_cluster_cfg );
    larlitecv::StopMuFoxTrot           stopmu_foxtrot( m_config.stopmu_foxtrot_cfg );
    
    if ( m_config.stopmu_cluster_cfg.save_pass_images || m_config.stopmu_cluster_cfg.dump_tagged_images  ) {
      std::stringstream ss;
      ss << "smc_r" << input.run << "_s" << input.subrun << "_e" << input.event << "_i" << input.entry;
      stopmu_cluster.setOpenCVImageStemName( ss.str() );
    }

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

    std::clock_t timer = std::clock();
    //output.stopmu_trackcluster_v = stopmu_cluster.findStopMuTracks( input.img_v, input.gapch_v, thrumu.tagged_v, output.stopmu_candidate_endpt_v );    
    output.stopmu_trackcluster_v = stopmu_foxtrot.findStopMuTracks( input.img_v, input.gapch_v, thrumu.tagged_v, thrumu.unused_spacepoint_v );
    m_time_tracker[kStopMuTracker] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;
    std::cout << "  Number of candidate StopMu tracks: " << output.stopmu_trackcluster_v.size() << std::endl;

    // make stopmu-tagged image
    for (size_t p=0; p<input.img_v.size(); p++) {
      larcv::Image2D stopmu_img( input.img_v.at(p).meta() );
      stopmu_img.paint(0);
      output.stopmu_v.emplace_back( std::move(stopmu_img) );
    }
    for ( size_t itrack=0; itrack<output.stopmu_trackcluster_v.size(); itrack++ ) {
      larlitecv::BMTrackCluster3D& track3d = output.stopmu_trackcluster_v.at(itrack);
      std::vector<larcv::Pixel2DCluster> trackpix_v = track3d.getTrackPixelsFromImages( input.img_v, input.gapch_v,
											m_config.stopmu_foxtrot_cfg.foxtrotalgo_cfg.pixel_thresholds,
											m_config.thrumu_tracker_cfg.tag_neighborhood,
											true, 0.3 );
      for (size_t p=0; p<trackpix_v.size(); p++) {
        const larcv::Pixel2DCluster& trackpixs = trackpix_v[p];
        for ( auto const& pix : trackpixs ) {
          output.stopmu_v.at(p).set_pixel( pix.Y(), pix.X(), 255 );
        }
      }
    }

    // copy track and pixels into separate containers.
    for ( auto const& bmtrack : output.stopmu_trackcluster_v ) {
      output.track_v.push_back( bmtrack.makeTrack() );
      std::vector< larcv::Pixel2DCluster > cluster_v;
      for ( auto const& track2d : bmtrack.plane_pixels )
        cluster_v.push_back( track2d );
      output.pixelcluster_v.emplace_back( std::move(cluster_v) );
    }

    std::cout << "End of StopMu Algo" << std::endl;

    return output;
  }

  CROIPayload TaggerCROIAlgo::runCROISelection( const InputPayload& input, const ThruMuPayload& thrumu, const StopMuPayload& stopmu ) {

    std::cout << "[TaggerCROIAlgo::runCROISelection]" << std::endl;

    CROIPayload output;

    ClusterGroupAlgo         clusteralgo(   m_config.untagged_cluster_cfg );
    ClusterGroupMatchingAlgo matchingalgo;
    TaggerFlashMatchAlgo     selectionalgo( m_config.croi_selection_cfg );

    // Make Thrumu/StopMu tagged and subtracted images
    //std::vector< larcv::Image2D > tagged_v;
    //std::vector< larcv::Image2D > subimg_v;
    for ( size_t p=0; p<input.img_v.size(); p++) {
      larcv::Image2D tagged( input.img_v.at(p).meta() );
      tagged.paint(0.0);
      larcv::Image2D sub( input.img_v.at(p) );

      for ( size_t r=0; r<tagged.meta().rows(); r++ ) {
        for ( size_t c=0; c<tagged.meta().cols(); c++ ) {
          // tagged image
          if ( thrumu.tagged_v.at(p).pixel(r,c)>0 || stopmu.stopmu_v.at(p).pixel(r,c)>0 )
	         tagged.set_pixel(r,c,255);

          // subtraction image: below threshold and tagged pixels get zeroed (for clustering)
          if ( sub.pixel(r,c)<10.0 || thrumu.tagged_v.at(p).pixel(r,c)>0 || stopmu.stopmu_v.at(p).pixel(r,c)>0 )
            sub.set_pixel(r,c,0.0);
        }
      }
      output.tagged_v.emplace_back( std::move(tagged) );
      output.subimg_v.emplace_back( std::move(sub) );
    }

    // ----------------------
    //  RUN ALGOS

    output.plane_groups_v = clusteralgo.MakeClusterGroups( input.img_v, input.gapch_v, output.tagged_v );

    output.vols_v = matchingalgo.MatchClusterGroups( output.subimg_v, output.plane_groups_v );

    // ------------------------------------------------------------------
    // WE COLLECT OUR CLUSTER DATA, forming TaggerFlashMatchData objects
    // Collect Pixel2D clusters for each plane. Define LArLite Track
    //std::vector< larlitecv::TaggerFlashMatchData > flashdata_v;

    // ThruMu
    for ( int itrack=0; itrack<(int)thrumu.trackcluster3d_v.size(); itrack++ ) {
      //const BMTrackCluster3D& track = thrumu.trackcluster3d_v.at(itrack);
      //std::vector< larcv::Pixel2DCluster > pixels;
      //for ( size_t p=0; p<input.img_v.size(); p++)
      //  pixels.push_back( track.plane_paths.at(p).pixelpath );
      //larlitecv::TaggerFlashMatchData thrumu_track( larlitecv::TaggerFlashMatchData::kThruMu, pixels, track.makeTrack() );
      larlitecv::TaggerFlashMatchData thrumu_track( larlitecv::TaggerFlashMatchData::kThruMu, thrumu.pixelcluster_v.at(itrack), thrumu.track_v.at(itrack) );
      output.flashdata_v.emplace_back( std::move(thrumu_track) );
    }

    // StopMu
    for ( int itrack=0; itrack<(int)stopmu.stopmu_trackcluster_v.size(); itrack++ ) {
      //const BMTrackCluster3D& track = stopmu.stopmu_trackcluster_v.at(itrack);
      //std::vector< larcv::Pixel2DCluster > pixels;
      //for ( size_t p=0; p<input.img_v.size(); p++)
      //  pixels.push_back( track.plane_paths.at(p).pixelpath );
      //larlitecv::TaggerFlashMatchData stopmu_track( larlitecv::TaggerFlashMatchData::kStopMu, pixels, track.makeTrack() );
      larlitecv::TaggerFlashMatchData stopmu_track( larlitecv::TaggerFlashMatchData::kStopMu, stopmu.pixelcluster_v.at(itrack), stopmu.track_v.at(itrack) );
      output.flashdata_v.emplace_back( std::move(stopmu_track) );
    }

    // Find Contained Clusters
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    for ( auto const& vol : output.vols_v ) {

      // there should be a selection here...
      // select by (1) good fraction (2) charge even-ness
      if ( vol.frac_good_slices<0.8 )
        continue;

      std::cout << "VOL: clgroup[" << vol._clustergroup_indices[0] << "," << vol._clustergroup_indices[1] << "," << vol._clustergroup_indices[2] << "] "
                << " numslices=" << vol.num_slices
                << " goodslices=" << vol.num_good_slices
                << " fracgood=" << vol.frac_good_slices
                << " planecharge=[" << vol.plane_charge[0] << "," << vol.plane_charge[1] << "," << vol.plane_charge[2] << "]"
                << std::endl;

      // we need to make a larlite::track object for this. we use the centroid of the slices
      std::vector< std::vector<float> > xyz_v;
      std::vector< float > slice_charge;
      for ( size_t islice=0; islice<vol.slices.size(); islice++ ) {

        const larlitecv::Slice_t& slice = vol.slices.at(islice);

        if ( slice.inside_tpc_boundary.size()==0 )
          continue;

        // for each slice volume, we are going to do the easy thing first and represent charge at centroid
        // of volume. not great, I know.
        float centroid[2] = {0.0};
        for ( auto const pt : slice.inside_tpc_boundary ) {
          for (int i=0; i<2; i++)
            centroid[i] += pt[i];
        }
        for (int i=0; i<2; i++) {
          centroid[i] /= float(slice.inside_tpc_boundary.size());
        }

        float tick = input.img_v.front().meta().pos_y( 0.5*(slice.row_interval[0]+slice.row_interval[1]) );
        float x = (tick-3200.0)*cm_per_tick;
        std::vector< float > xyz(3,0);
        xyz[0] = x;
        xyz[1] = centroid[0];
        xyz[2] = centroid[1];

        xyz_v.emplace_back( std::move(xyz) );

        // slice charge
        float ave_charge = 0.;
        for ( int i=0; i<3; i++ ) {
          ave_charge += slice.plane_charge[i];
        }
        ave_charge /= 3.0;
        slice_charge.push_back( ave_charge );
      }

      larlite::track contained_track;
      for ( int ipt=0; ipt<(int)xyz_v.size()-1; ipt++ ) {
        const std::vector<float>& xyz = xyz_v.at(ipt);
        TVector3 pos( xyz[0], xyz[1], xyz[2] );
        const std::vector<float>& xyz_next = xyz_v.at(ipt+1);
        float dir[3] = {0};
        float norm = 0.;
        for (int i=0; i<3; i++) {
          dir[i] = xyz_next[i] - xyz[i];
          norm += dir[i]*dir[i];
        }
        norm = sqrt(norm);
        for (int i=0; i<3; i++)
          dir[i] /= norm;
        TVector3 dirv( dir[0], dir[1], dir[2] );
        contained_track.add_vertex( pos );
        contained_track.add_direction( dirv );
        contained_track.add_momentum( slice_charge.at(ipt) );

        if ( ipt==(int)xyz_v.size()-2) {
          TVector3 nextpos( xyz_next[0], xyz_next[1], xyz_next[2] );
          contained_track.add_vertex( nextpos );
          contained_track.add_direction( dirv );
          contained_track.add_momentum( slice_charge.at(ipt+1) );
        }
      }
      if ( contained_track.NumberTrajectoryPoints()==0 )
        continue;
      larlitecv::TaggerFlashMatchData contained_cluster( larlitecv::TaggerFlashMatchData::kUntagged, vol.m_plane_pixels, contained_track );
      output.flashdata_v.emplace_back( std::move(contained_cluster) );
    }// end of vol loop

    // --------------------------------------------------------------------//
    // SELECT ROIs

    output.flashdata_selected_v.resize( output.flashdata_v.size(), 0 );
    std::vector<larcv::ROI> selected_rois = selectionalgo.FindFlashMatchedContainedROIs( output.flashdata_v, input.opflashes_v, output.flashdata_selected_v );

    output.croi_v.clear();
    for ( size_t itrack=0; itrack<output.flashdata_v.size(); itrack++ ){
      const larlite::track& track3d = output.flashdata_v.at(itrack).m_track3d;
      if ( output.flashdata_selected_v.at(itrack)==0 || track3d.NumberTrajectoryPoints()==0)
        continue;

      larcv::ROI croi = output.flashdata_v.at(itrack).MakeROI( input.img_v, m_config.croi_selection_cfg.bbox_pad , true );

      std::cout << "[Selected CROI]" << std::endl;
      for ( size_t p=0; p<3; p++ ) {
        std::cout << "  " << croi.BB(p).dump() << std::endl;
      }

      output.croi_v.emplace_back( std::move(croi) );
    }

    // ------------------------------------------------------------------------//
    // Make Combined Tagged Image

    if ( m_config.croi_write_cfg.get<bool>("WriteCombinedTaggedImage") ) {


      for ( size_t p=0; p<input.img_v.size(); p++ ) {
        larcv::Image2D combined( input.img_v.at(p).meta() );
        combined.paint(0.0);
        output.combined_v.emplace_back( std::move(combined) );
      }

      for ( size_t itrack=0; itrack<output.flashdata_v.size(); itrack++ ) {
        const larlitecv::TaggerFlashMatchData& flashdata = output.flashdata_v.at(itrack);
        if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 ) {
          continue;
        }
        const std::vector<larcv::Pixel2DCluster>& pixels = flashdata.m_pixels;
        int tagval = 10.0*((int)flashdata.m_type + 1); // 0= nothing, 10 = thrumu, 20=stopmu, 30=untagged/contained, 40=selected CROI
        if ( output.flashdata_selected_v.at(itrack)==1 ) {
          tagval = 40.0;
        }
        for ( size_t p=0; p<output.combined_v.size(); p++ ) {
          const larcv::ImageMeta& meta = input.img_v.at(p).meta();
          for ( auto const& pix : pixels.at(p) ) {
            if ( (int)pix.X()<0 ||  pix.X()>=meta.cols() || (int)pix.Y()<0 || pix.Y()>=meta.rows() ) continue;
            int pixval = output.combined_v.at(p).pixel( pix.Y(), pix.X() );
            if ( tagval > pixval )
              output.combined_v.at(p).set_pixel( pix.Y(), pix.X(), tagval );
          }
        }
      }
    }

    if ( m_config.croi_write_cfg.get<bool>("WriteTrackOpFlashes")) {
      // collect flash info
      output.track_opflash_v.clear();
      for ( auto const& ophypo : selectionalgo.getOpFlashHypotheses() )
        output.track_opflash_v.push_back( ophypo );
    }

    return output;
  }

  void TaggerCROIAlgo::printTimeTracker( int num_events ) {
    const std::string stage_names[] = { "ThruMuConfig", "ThruMuBMT", "ThruMuFlash", "ThruMuTracker", "StopMuTracker", "Untagged", "CROI" };
    float tot_time = 0.;
    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "TaggerCROIAlgoConfig::printTimeTracker" << std::endl;
    std::cout << "Number of Events: " << num_events << std::endl;
    for (int i=0; i<kNumStages; i++) {
      std::cout << stage_names[i] << " : " << m_time_tracker[i] << " secs";
      if ( num_events>0 )
        std::cout << "  " << m_time_tracker[i]/float(num_events) << " secs/event";
      std::cout << std::endl;
      tot_time += m_time_tracker[i];
    }
    std::cout << "Total: " << tot_time << " secs";
    if ( num_events>0 )
      std::cout << "  " << tot_time/float(num_events) << " secs/event";
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
  }
}
