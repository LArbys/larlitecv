#include <iostream>
#include <string>
#include <vector>

// larcv
#include "DataFormat/EventROI.h"
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/chstatus.h"
#include "LArUtil/LArProperties.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "UntaggedClustering/ClusterGroupAlgo.h"
#include "UntaggedClustering/ClusterGroupMatchingTypes.h"
#include "UntaggedClustering/ClusterGroupMatchingAlgo.h"
#include "ContainedROI/TaggerFlashMatchTypes.h"
#include "ContainedROI/TaggerFlashMatchAlgoConfig.h"
#include "ContainedROI/TaggerFlashMatchAlgo.h"

#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif


int main( int nargs, char** argv ) {

  std::cout << "[[ ROI SELECTION ]]" << std::endl;

  if ( nargs!=2) {
    std::cout << "./run_roi_selection [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet cfg      = cfg_root.get<larcv::PSet>( "ContainedROI" );

  larlitecv::DataCoordinator dataco_source;
  dataco_source.set_filelist( cfg.get<std::string>( "SourceLArCVFilelist" ), "larcv" );   // original image/(segment image if we need it)
  dataco_source.set_filelist( cfg.get<std::string>( "SourceLArLiteFilelist"), "larlite"); // source larlite file

  // Deprecate?
  //larlitecv::DataCoordinator dataco_thrumu;
  //dataco_thrumu.set_filelist( , "larcv" ); // thrumu-tagger info, larcv
  //dataco_thrumu.set_filelist( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.set_filelist( cfg.get<std::string>( "StopMuLArCVFilelist"),   "larcv" );   //stopmu-tagger output
  dataco_stopmu.set_filelist( cfg.get<std::string>( "StopMuLArLiteFilelist"), "larlite" ); //stopmu-tagger output

  larlitecv::DataCoordinator dataco_output;

  dataco_source.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  //dataco_thrumu.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_stopmu.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_output.configure( "containedroi.cfg", "StorageManagerOutput", "IOManagerOutput", "ContainedROI" );

  dataco_source.initialize();
  //dataco_thrumu.initialize();
  dataco_stopmu.initialize();
  dataco_output.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  //std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;

  // configuration parameters
  larcv::PSet cluster_pset   = cfg.get<larcv::PSet>("ContainedGroupAlgo");
  larcv::PSet selection_pset = cfg.get<larcv::PSet>("TaggerFlashMatchAlgo");
  bool isMC = cfg.get<bool>("IsMC");
  std::string tpc_image_producername    = cfg.get<std::string>("StopMuInputLArCVImages");
  std::string gapchs_image_producername = cfg.get<std::string>("StopMuBadChLArCVImages");
  std::string stopmutagged_image_producername = cfg.get<std::string>("StopMuTaggedImages");
  std::string thrumutagged_image_producername = cfg.get<std::string>("ThruMuTaggedImages");
  std::vector<std::string> flashproducers = cfg.get< std::vector<std::string> >("OpFlashProducers");

  // ----------------------------------------------------------------

  // Cluster Group
  larlitecv::ClusterGroupAlgoConfig cluster_cfg = larlitecv::ClusterGroupAlgoConfig::MakeClusterGroupAlgoConfigFromPSet( cluster_pset );  
  larlitecv::ClusterGroupAlgo clusteralgo( cluster_cfg );

  // Cluster Group Match
  larlitecv::ClusterGroupMatchingAlgo matchingalgo;

  // Selection Algo
  larlitecv::TaggerFlashMatchAlgoConfig selection_cfg = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( selection_pset );
  larlitecv::TaggerFlashMatchAlgo selectionalgo( selection_cfg );

  // Event Loop
  int nentries        = dataco_stopmu.get_nentries("larcv");
  int user_nentries   = cfg.get<int>("NumEntries",-1);
  int user_startentry = cfg.get<int>("StartEntry",-1);
  int start_entry     = 0;
  if ( user_startentry>=0 ) {
    start_entry = user_startentry;
  }
  int end_entry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+start_entry<nentries )
      end_entry = user_nentries+start_entry;
  }

  for ( int ientry=start_entry; ientry<end_entry; ientry++ ) {

    // load same event in all three data sets
    dataco_stopmu.goto_entry(ientry,"larcv");
    int run,subrun,event;
    dataco_stopmu.get_id(run,subrun,event);
    std::cout << "entry " << ientry << std::endl;
    std::cout << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    //dataco_thrumu.goto_event(run,subrun,event,"larcv");
    dataco_source.goto_event(run,subrun,event,"larcv");
    dataco_output.set_id(run, subrun, event);

    // algo.m_run = run;
    // algo.m_subrun = subrun;
    // algo.m_event = event;
    // algo.m_entry = ientry;

    // ------------------------------------------------------------------------------------------//
    // STOPMU DATA PRODUCTS: images, pixels, and tracks

    // images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D, tpc_image_producername);
    larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D, thrumutagged_image_producername);
    larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D, stopmutagged_image_producername); 
    larcv::EventImage2D* ev_gapchs = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D, gapchs_image_producername);

    // retrieve the thrumu and stopmu clusters as well
    larcv::EventPixel2D*  ev_thrumu_pixels = (larcv::EventPixel2D*) dataco_stopmu.get_larcv_data(larcv::kProductPixel2D,   "thrumu2d" );
    larcv::EventPixel2D*  ev_stopmu_pixels = (larcv::EventPixel2D*) dataco_stopmu.get_larcv_data(larcv::kProductPixel2D,   "stopmupixels");

    // tracks
    larlite::event_track* ev_thrumu_tracks = (larlite::event_track*)dataco_stopmu.get_larlite_data( larlite::data::kTrack, "thrumu3d" );
    larlite::event_track* ev_stopmu_tracks = (larlite::event_track*)dataco_stopmu.get_larlite_data( larlite::data::kTrack, "stopmu3d");

    // MC information: ROI and segment image
    larcv::EventROI* rois       = nullptr;
    larcv::EventImage2D* ev_seg = nullptr;
    if ( isMC ) {
      // source
      rois   = (larcv::EventROI*)    dataco_source.get_larcv_data(larcv::kProductROI, "tpc" );
      ev_seg = (larcv::EventImage2D*)dataco_source.get_larcv_data(larcv::kProductImage2D,"segment");
    }

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
    const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();
    const std::vector<larcv::Image2D>& gapchs_v = ev_gapchs->Image2DArray();

    // ------------------------------------------------------------------------------------------//
    // LARLITE OPFLASH DATA

    std::vector< larlite::event_opflash* > opflash_containers;
    for ( auto &flashproducer : flashproducers ) {
      std::cout << "search for flash hits from " << flashproducer << ". ";
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco_source.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << " has " << opdata->size() << " flashes" << std::endl;
      opflash_containers.push_back( opdata );
    }

    // ------------------------------------------------------------------------------------------//
    // COMBINE STOPMU/THRUMU TAGGER INFO

    std::vector< larcv::Image2D > tagged_v;
    std::vector< larcv::Image2D > subimg_v;    
    for ( size_t p=0; p<imgs_v.size(); p++) {
      larcv::Image2D tagged( imgs_v.at(p).meta() );
      tagged.paint(0.0);
      larcv::Image2D sub( imgs_v.at(p) );

      for ( size_t r=0; r<tagged.meta().rows(); r++ ) {
        for ( size_t c=0; c<tagged.meta().cols(); c++ ) {
          // tagged image
          if ( thrumu_v.at(p).pixel(r,c)>0 || stopmu_v.at(p).pixel(r,c)>0 )
            tagged.set_pixel(r,c,255);

          // subtraction image: below threshold and tagged pixels get zeroed (for clustering)
          if ( sub.pixel(r,c)<10.0 || thrumu_v.at(p).pixel(r,c)>0 || stopmu_v.at(p).pixel(r,c)>0 )
            sub.set_pixel(r,c,0.0);
        }
      }
      tagged_v.emplace_back( std::move(tagged) );
      subimg_v.emplace_back( std::move(sub) );
    }

    // ------------------------------------------------------------------------------------------//
    // WE MAKE THE UNTAGGED/CONTAINED (3D) CLUSTERS and MATCH OVER PLANES

    std::vector< larlitecv::PlaneClusterGroups > plane_groups_v = clusteralgo.MakeClusterGroups( imgs_v, gapchs_v, tagged_v );

    std::vector<larlitecv::ChargeVolume> vols_v = matchingalgo.MatchClusterGroups( subimg_v, plane_groups_v );

    // ------------------------------------------------------------------------------------------//
    // OPTIONAL: DUMP JPG

    // image the output of clustergroup matcher
    std::vector<cv::Mat> cv_matched;
    for ( size_t p=0; p<imgs_v.size(); p++) {
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( imgs_v.at(p), 10, 60 );
      cv_matched.emplace_back( std::move(cvimg) );
    }
    matchingalgo.labelCVImageWithMatchedClusters( cv_matched, imgs_v, vols_v, 0.75 );
    for ( size_t p=0; p<cv_matched.size(); p++) {
      cv::Mat& cvimg = cv_matched.at(p);
      std::stringstream ss;
      ss << "clustermatch_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".jpg";
      imwrite( ss.str(), cvimg );
    }   

    // ------------------------------------------------------------------------------------------//
    // WE COLLECT OUR CLUSTER DATA, forming TaggerFlashMatchData objects

    // ThruMu
    if ( ev_thrumu_pixels->Pixel2DClusterArray(0).size()!=ev_thrumu_tracks->size() )
      throw std::runtime_error("size of ThruMu pixelclusters and tracks are not the same!");

    std::vector< larlitecv::TaggerFlashMatchData > flashdata_v;
    for ( int itrack=0; itrack<(int)ev_thrumu_tracks->size(); itrack++ ) {
      std::vector< larcv::Pixel2DCluster > pixels;
      for ( size_t p=0; p<imgs_v.size(); p++)
        pixels.push_back( ev_thrumu_pixels->Pixel2DClusterArray(p).at(itrack) );
      larlitecv::TaggerFlashMatchData thrumu_track( larlitecv::TaggerFlashMatchData::kThruMu, pixels, ev_thrumu_tracks->at(itrack) );
      flashdata_v.emplace_back( std::move(thrumu_track) );
    }

    // StopMu
    if ( ev_stopmu_pixels->Pixel2DClusterArray(0).size()!=ev_stopmu_tracks->size() )
      throw std::runtime_error("size of StopMu pixelclusters and tracks are not the same!");

    for ( int itrack=0; itrack<(int)ev_stopmu_tracks->size(); itrack++ ) {
      std::vector< larcv::Pixel2DCluster > pixels;
      for ( size_t p=0; p<imgs_v.size(); p++)
        pixels.push_back( ev_stopmu_pixels->Pixel2DClusterArray(p).at(itrack) );
      larlitecv::TaggerFlashMatchData stopmu_track( larlitecv::TaggerFlashMatchData::kStopMu, pixels, ev_stopmu_tracks->at(itrack) );
      flashdata_v.emplace_back( std::move(stopmu_track) );
    }

    // Find Contained Clusters
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;        
    for ( auto const& vol : vols_v ) {

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

        float tick = imgs_v.front().meta().pos_y( 0.5*(slice.row_interval[0]+slice.row_interval[1]) );
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
      larlitecv::TaggerFlashMatchData contained_cluster( larlitecv::TaggerFlashMatchData::kUntagged, vol.m_plane_pixels, contained_track );
      flashdata_v.emplace_back( std::move(contained_cluster) );
    }// end of vol loop

    // --------------------------------------------------------------------//
    // SELECT ROIs

    std::vector<int> flashdata_selected_v( flashdata_v.size(), 0 );
    std::vector<larcv::ROI> selected_rois = selectionalgo.FindFlashMatchedContainedROIs( flashdata_v, opflash_containers, flashdata_selected_v );

    std::vector<larcv::ROI> croi_v;
    for ( size_t itrack=0; itrack<flashdata_v.size(); itrack++ ){
      const larlite::track& track3d = flashdata_v.at(itrack).m_track3d;
      if ( flashdata_selected_v.at(itrack)==0 || track3d.NumberTrajectoryPoints()==0)
        continue;
      larcv::ROI croi = flashdata_v.at(itrack).MakeROI( imgs_v, true );

      std::cout << "[Selected CROI]" << std::endl;
      for ( size_t p=0; p<3; p++ ) {
        std::cout << "  " << croi.BB(p).dump() << std::endl;
      }
      croi_v.emplace_back( std::move(croi) );
    }


    // --------------------------------------------------------------------//
    // SAVE DATA TO OUTPUT FILE
    // We must transfer the data explicitly

    // Save New Untagged Clusters
    larlite::event_track* evout_tracks_thrumu   = (larlite::event_track*)dataco_output.get_larlite_data( larlite::data::kTrack, "thrumu3d" );
    larlite::event_track* evout_tracks_stopmu   = (larlite::event_track*)dataco_output.get_larlite_data( larlite::data::kTrack, "stopmu3d" );
    larlite::event_track* evout_tracks_untagged = (larlite::event_track*)dataco_output.get_larlite_data( larlite::data::kTrack, "untagged3d" );
    larlite::event_track* evout_tracks_selected = (larlite::event_track*)dataco_output.get_larlite_data( larlite::data::kTrack, "croi3d" );

    for ( auto const& track : *ev_thrumu_tracks ) {
      evout_tracks_thrumu->push_back( track );
    }
    for ( auto const& track : *ev_stopmu_tracks ) {
      evout_tracks_stopmu->push_back( track );
    }
    for ( size_t itrack=0; itrack<flashdata_v.size(); itrack++) {
      auto const& flashdata = flashdata_v.at(itrack);
      if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 )
        continue;
      if ( flashdata.m_type==larlitecv::TaggerFlashMatchData::kUntagged ) {
        evout_tracks_untagged->push_back( flashdata.m_track3d );
      }
      if ( flashdata_selected_v[itrack] ) {
        evout_tracks_selected->push_back( flashdata.m_track3d );
      }
    }


    // ROIs, StopMu Tracks and clusters, ThruMu Tracks and Clusters, Truth ROI, Bad channels
    larcv::EventROI* out_ev_roi = (larcv::EventROI*)dataco_output.get_larcv_data( larcv::kProductROI, "croi" );
    out_ev_roi->Set( croi_v );

    // save the image and bad channels
    larcv::EventImage2D* out_ev_images = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "tpc" );
    larcv::EventImage2D* out_ev_gapchs = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "gapchs" );
    for ( auto& img : imgs_v )
      out_ev_images->Append( img );
    for ( auto& img : gapchs_v )
      out_ev_gapchs->Append( img );

    // Save Pixels
    larcv::EventPixel2D* evout_pixels_thrumu   = (larcv::EventPixel2D*)dataco_output.get_larcv_data( larcv::kProductPixel2D, "thrumupixels" );
    larcv::EventPixel2D* evout_pixels_stopmu   = (larcv::EventPixel2D*)dataco_output.get_larcv_data( larcv::kProductPixel2D, "stopmupixels" );
    larcv::EventPixel2D* evout_pixels_untagged = (larcv::EventPixel2D*)dataco_output.get_larcv_data( larcv::kProductPixel2D, "untaggedpixels" );
    larcv::EventPixel2D* evout_pixels_croi     = (larcv::EventPixel2D*)dataco_output.get_larcv_data( larcv::kProductPixel2D, "croipixels" );

    for ( size_t itrack=0; itrack<flashdata_v.size(); itrack++ ) {
      auto const& trackdata = flashdata_v.at(itrack);

      if ( trackdata.m_track3d.NumberTrajectoryPoints()==0 )
        continue;

      if ( flashdata_selected_v.at(itrack)==0 ) {
        if ( trackdata.m_type==larlitecv::TaggerFlashMatchData::kThruMu ) {
          for ( size_t p=0; p<3; p++ )
            evout_pixels_thrumu->Append( (larcv::PlaneID_t)p, trackdata.m_pixels.at(p) );
        }
        else if ( trackdata.m_type==larlitecv::TaggerFlashMatchData::kStopMu ) {
          for ( size_t p=0; p<3; p++ )
            evout_pixels_stopmu->Append( (larcv::PlaneID_t)p, trackdata.m_pixels.at(p) );      
        }
        else if ( trackdata.m_type==larlitecv::TaggerFlashMatchData::kUntagged ) {
          for ( size_t p=0; p<3; p++ )
            evout_pixels_untagged->Append( (larcv::PlaneID_t)p, trackdata.m_pixels.at(p) );
        }
      }
      else {
        for ( size_t p=0; p<3; p++ ) {
          evout_pixels_croi->Append( (larcv::PlaneID_t)p, trackdata.m_pixels.at(p) );
        }
      }
    }

    // Produce Tagged Image
    std::vector< larcv::Image2D > combined_tagged_v;
    for ( size_t p=0; p<imgs_v.size(); p++) {
      larcv::Image2D combined_tagged( imgs_v.at(p).meta() );
      combined_tagged.paint(0);

      for ( size_t itrack=0; itrack<flashdata_v.size(); itrack++ ) {
        auto const& trackdata = flashdata_v.at(itrack);
        auto const& trackpixels = trackdata.m_pixels.at(p);
        int tag_val = 0;
        if ( flashdata_selected_v.at(itrack)==1 )
          tag_val = 1;
        else
          tag_val = 2;

        for ( auto const& pix : trackpixels ) {
          combined_tagged.set_pixel( pix.Y(), pix.X(), tag_val );
        }
      }

      combined_tagged_v.emplace_back( std::move(combined_tagged) );
    }
    larcv::EventImage2D* out_ev_combined_tagged = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "combinedtagged" );
    out_ev_combined_tagged->Emplace( std::move(combined_tagged_v) );

    // copy over truth quantities
    if ( isMC ) {
      larcv::EventImage2D* out_ev_segmnt = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "segment");
      larcv::EventROI*     out_ev_roi    = (larcv::EventROI*)    dataco_output.get_larcv_data( larcv::kProductROI, "tpc");
      for ( auto& img : ev_seg->Image2DArray() )
         out_ev_segmnt->Append( img );
      out_ev_roi->Set( rois->ROIArray() );
    }

    // copy over opdigits
    larlite::event_opdetwaveform* ev_opdigit     = (larlite::event_opdetwaveform*) dataco_source.get_larlite_data( larlite::data::kOpDetWaveform, "saturation" );
    larlite::event_opdetwaveform* out_ev_opdigit = (larlite::event_opdetwaveform*) dataco_output.get_larlite_data( larlite::data::kOpDetWaveform, "saturation" );
    for ( auto const& opdigit : *ev_opdigit )
      out_ev_opdigit->push_back( opdigit );

    // copy over trigger
    larlite::trigger* ev_trigger    = (larlite::trigger*) dataco_source.get_larlite_data( larlite::data::kTrigger, "triggersim" );
    larlite::trigger* evout_trigger = (larlite::trigger*) dataco_output.get_larlite_data( larlite::data::kTrigger, "triggersim" );


    dataco_output.save_entry();

  }//end of event loop

  dataco_output.finalize();

  return 0;
}
