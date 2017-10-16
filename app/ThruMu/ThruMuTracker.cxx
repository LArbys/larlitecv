#include "ThruMuTracker.h"

#include "RadialEndpointFilter.h"
#include "AStarNodes2BMTrackCluster3D.h"
#include "ThruMuFoxExtender.h"
#include "PushBoundarySpacePoint.h"

// larcv 
#include "UBWireTool/UBWireTool.h"
#include "Reco3D/AStar3DAlgo.h"


namespace larlitecv {

  ThruMuTracker::ThruMuTracker( const ThruMuTrackerConfig& config )
    : m_config(config)
  {}

  // Add additional arguments: 
  //                           
  void ThruMuTracker::makeTrackClusters3D( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector<larcv::Image2D>& img_v,  
					   const std::vector<larcv::Image2D>& badchimg_v, const std::vector< const BoundarySpacePoint* >& spacepts,
					   std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
					   std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices, const std::vector< larlite::event_opflash* >& opflashsets ) {

    // This method takes in the list of boundaryspacepoints and pairs them up in order to try to find through-going muons
    // input:
    //   flash_match_config': configuration for 'GeneralFlashMatchAlgo' instance within the 'flashMatchTracks' function.    
    //   img_v: vector images, one for each plane. assumed to be in U,V,Y order
    //   badchimg_v: vector images, one for each plane. assumed to be in U,V,Y order    
    //   spacepts: all possible boundary crossing points
    // output:
    //   trackclusters: thrumu tracks. collection of space points and pixel2dclusters
    //   tagged_v: pixels tagged as thrumu.  the value indicates the pass when the muon was tagged.
    //   used_endpoints_indices: indicates which endpoints in spacepts was used. 0=not used. 1=used.

    //std::cout << "At the top of 'makeTrackClusters3D' function in 'ThruMuTracker'." << std::endl;

    // Declare a set, 'impossible_match_endpoints', meant to save the indices of endpoints that together constituted an impossible match.
    // This removes from consideration endpoints that were found to fail the flashmatching stage of reconstruction.                                                                                   
    std::vector< std::vector< BoundarySpacePoint > > impossible_match_endpoint_v;
    impossible_match_endpoint_v.clear();

    // Declare a new vector, 'track_endpoint_indices', based off the indices in 'spacepts', for the endpoints that are used to make tracks.
    std::vector < std::vector< BoundarySpacePoint >  > track_endpoint_v;
    track_endpoint_v.clear();

    // Declare a vector that contains the 'BoundaryFlashIndex' information for the flash that the endpoints correspond to.
    std::vector< std::vector< BoundaryFlashIndex > > track_endpoint_flash_v;
    track_endpoint_flash_v.clear();

    // Declare a vector for the indices of the boundary of the detector that this flash corresponds to.
    std::vector< std::vector< BoundaryEnd_t >  > track_endpoint_boundary_type_idx_v;
    track_endpoint_boundary_type_idx_v.clear();

    // Declare a set of the 'BoundaryFlashIndex' information of the flashes that have already been well-matched to a track.
    std::vector< BoundaryFlashIndex > already_matched_flash_v;
    already_matched_flash_v.clear();

    // Declare a vector for the well-matched tracks 'well_matched_tracks_idx_v'.
    std::vector < int > well_matched_tracks_idx_v;
    well_matched_tracks_idx_v.clear();

    const int nendpts = (int)spacepts.size();
    used_endpoints_indices.resize( nendpts, 0 );

    // pair up containers
    int ntotsearched = 0;
    int npossible = 0;

    // compress image for astar. store in member class. we'll use these again in runAStarTagger
    m_img_compressed_v.clear();
    m_badch_compressed_v.clear();
    int downsampling_factor = m_config.downsampling_factor;
    for (size_t p=0; p<img_v.size(); p++) {
      larcv::Image2D img_compressed( img_v[p] );
      larcv::Image2D badch_compressed( badchimg_v[p] );
      img_compressed.compress( img_v[p].meta().rows()/downsampling_factor,
                                img_v[p].meta().cols()/downsampling_factor,
                                (larcv::Image2D::CompressionModes_t)m_config.compression_mode );
      badch_compressed.compress( img_v[p].meta().rows()/downsampling_factor,
                                img_v[p].meta().cols()/downsampling_factor,
                                (larcv::Image2D::CompressionModes_t)m_config.compression_mode );
      m_img_compressed_v.emplace_back( std::move(img_compressed) );
      m_badch_compressed_v.emplace_back( std::move(badch_compressed) );
    }

    // tagged image
    if ( tagged_v.size()==0 )  {
      for (size_t p=0; p<img_v.size(); p++) {
        larcv::Image2D tagged( img_v[p].meta() );
        tagged.paint(0);
        tagged_v.emplace_back( std::move(tagged) );
      }
    }


    // poor-man's profiling
    const clock_t begin_time = clock();

    // Set a boolean, 'single_pass', to 'true'.
    bool single_pass = true;

    // Set a boolean, 'anode_and_cathode' only, to 'true'.
    bool anode_and_cathode_only = true;

    std::vector<int> tracks_per_pass;
    for (int ipass=0; ipass<m_config.num_passes; ipass++) {
      
      // Clear 'well_matched_tracks_idx_v' at the start of the pass.
      well_matched_tracks_idx_v.clear();
      
      const ThruMuTrackerConfig::ThruMuPassConfig& passcfg = m_config.pass_configs.at(ipass);
      runPass( ipass, passcfg, spacepts, img_v, badchimg_v, tagged_v, used_endpoints_indices, trackclusters,
	       track_endpoint_flash_v, track_endpoint_boundary_type_idx_v, track_endpoint_v );

      int tracks_in_pass = trackclusters.size();
      for ( int i = int( tracks_per_pass.size() - 1 ); i > -1; --i ) {
        tracks_in_pass -= tracks_per_pass.at( i );
      }
      tracks_per_pass.push_back( tracks_in_pass );

      std::cout << "The number of tracks reconstructed in Pass #" << ipass << " before flashmatching = " << tracks_in_pass << "." << std::endl;
      std::cout << "The size of 'well_matched_flash_v' before flashmatching in Pass #" << ipass << " = " << already_matched_flash_v.size() << "." << std::endl;

      // Depending if this is greater than the first pass, then set 'anode_and_cathode_only' to 'false'.
      if ( ipass > 0 ) 
	anode_and_cathode_only = false;

      // Use the flag from the config file to determine if you want to flashmatch the tracks that survive this pass of the thrumu tracker.
      if (m_config.thrumu_flashmatch == true ) {

        // Call a function that uses all of the flashmatching infrastructure developed in 'GeneralFlashMatchAlgo'.
	flashMatchTracks( flash_match_config, spacepts, opflashsets, trackclusters, impossible_match_endpoint_v, 
			  already_matched_flash_v, well_matched_tracks_idx_v, tracks_in_pass, track_endpoint_flash_v, 
			  track_endpoint_boundary_type_idx_v, track_endpoint_v, anode_and_cathode_only, single_pass );	
	sortOutBadTracksAndEndpoints( trackclusters, track_endpoint_v, track_endpoint_flash_v, track_endpoint_boundary_type_idx_v, well_matched_tracks_idx_v, tracks_per_pass, tracks_per_pass.at( tracks_per_pass.size() - 1 ), single_pass ); 

	//std::cout<< "The number of tracks reconstructed in Pass #" << ipass << " after flashmatching = " << tracks_per_pass.at( tracks_per_pass.size() - 1 ) << "." << std::endl;

      }
     
      std::cout << "The number of flashes considered well-matched after pass #" << ipass << " = " << already_matched_flash_v.size() << "." << std::endl;
    }

    // Reset 'single_pass' to 'false', 'anode_and_cathode' only to 'false',  and reuse the 'flashMatchTracks' function.
    single_pass            = false;
    anode_and_cathode_only = false;

    // Set a dummy variable for 'num_of_tracks_added_in_pass', since that number will not be used.
    int num_of_tracks_added_in_pass = 0;

    if (m_config.thrumu_flashmatch == true ) {
      flashMatchTracks( flash_match_config, spacepts, opflashsets, trackclusters, impossible_match_endpoint_v,
			already_matched_flash_v, well_matched_tracks_idx_v, num_of_tracks_added_in_pass, track_endpoint_flash_v, 
			track_endpoint_boundary_type_idx_v, track_endpoint_v, anode_and_cathode_only, single_pass );
      sortOutBadTracksAndEndpoints( trackclusters, track_endpoint_v, track_endpoint_flash_v, track_endpoint_boundary_type_idx_v, well_matched_tracks_idx_v, tracks_per_pass, trackclusters.size(), single_pass ); 
    }

    std::cout<< "The number of flashes considered well-matched at the end of the function = " << already_matched_flash_v.size() << "." << std::endl;

    // Declare an iterator for the number of reconstructed tracks at the end.
    int num_reco_tracks_at_end = 0;

    // tag the pixels
    for (int itrack=0; itrack<(int)trackclusters.size(); itrack++ ) {
      // At this point, in addition to tagging the pixels, I can remove those tracks from consideration that have an empty path.  I no longer rely on the symmetry between the number of tracks in 'trackclusters' and the number of tracks reconstructed in the pass.
      if (trackclusters.at( itrack ).path3d.size() == 0 ) {
	trackclusters.erase( trackclusters.begin() + itrack );
	continue;
      }
      
      ++num_reco_tracks_at_end;
      
      BMTrackCluster3D& track3d = trackclusters[itrack];
      track3d.markImageWithTrack( img_v, badchimg_v, m_config.pixel_threshold, m_config.tag_neighborhood, tagged_v, 0.3, itrack+1 );
    }
    
    std::cout << "Number of tracks reconstructed at the end = " << num_reco_tracks_at_end << "." << std::endl;

    //if ( m_config.verbosity>0 ) {
    //  for (int ipass=0; ipass<m_config.num_passes; ipass++) {
    //    std::cout << "Number of tracks found in ThruMu pass #" << ipass << ": " << tracks_per_pass.at(ipass) << std::endl;
    //  }
    //}

  }

  // Update this function to pass the index of the flash that is matched to the endpoints of the track along with the indices themselves.
  // New arguments: 
  //   track_endpoint_flash_idx_v: This vector of ints corresponds to the index of the flash that is matched to the track endpoints at this point in the vector.
  //   track_endpoint_boundary_type_idx_v: This vector of ints contains the producer of the flash that corresponds to the endpoints of the track in the 'trackclusters' vector.
  void ThruMuTracker::runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
			       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
			       std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters,
			       std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v,
			       std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v, std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v) {

    // Make a pass to try and connect some end points.
    // A pass consists of running the linear3d tagger and the astar3d tagger.
    // Cuts on the quality of each determine which, if any, of the taggers are run and produce the track.
    // input:
    //   passid: id for the current pass
    //   passcfg: configuration for the pass
    //   spacepts: boundary end points
    //   img_v: original images
    //   tagged_v: tagged images
    // input/output:
    //   used_endpoints_indices: marks which endpoints have been used
    //   trackclusters: holds thrumu tracks made

    if ( m_config.verbosity>0 ) {
      std::cout << "Run Pass #" << passid << ": "
		            << " radialfilter=" << passcfg.run_radial_filter
		            << " linear=" << passcfg.run_linear_tagger
		            << " astar=" << passcfg.run_astar_tagger
		            << std::endl;
    }

    const int nendpts = (int)spacepts.size();
    std::vector<larlitecv::BMTrackCluster3D> pass_track_candidates;
    std::vector< std::vector<int> > pass_end_indices;

    // Declare the vector for the matched track indices up here.
    std::vector< BoundaryFlashIndex > single_track_endpoint_flash_v;
    single_track_endpoint_flash_v.clear();

    std::vector< larlitecv::BoundaryEnd_t > single_track_boundary_type_idx_v;
    single_track_boundary_type_idx_v.clear();

    std::vector< BoundarySpacePoint > single_track_endpoint_v;
    single_track_endpoint_v.clear();

    RadialEndpointFilter radialfilter;

    int num_tracked = 0;

    for (int i=0; i<nendpts; i++) {
      if ( used_endpoints_indices.at(i)==1 ) continue;
      const BoundarySpacePoint& pts_a = *(spacepts[i]);

      bool use_a = true;
      if ( passcfg.run_radial_filter ) {
        int num_segs_a = 0;
        bool within_line_a =  radialfilter.isWithinStraightSegment( pts_a.pos(), img_v, badchimg_v, passcfg.radial_cfg, num_segs_a );
        bool pass_num_segs =  passcfg.radial_cfg.min_segments<=num_segs_a && passcfg.radial_cfg.max_segments>=num_segs_a;
        if ( m_config.verbosity> 1 ) {
          std::cout << "Endpoint #" << i << ": "
                                    << " " << pts_a.printImageCoords( img_v.front().meta() ) << " "
                                    << " within_line=" << within_line_a << " num_segs=" << num_segs_a << " pass_num_segs=" << pass_num_segs
                                    << " run=" << (!within_line_a && pass_num_segs)
                                    << std::endl;
        }

        if ( !within_line_a &&  pass_num_segs )
          use_a = true;
        else
          use_a = false;
      }

      if ( !use_a ) {
        if ( m_config.verbosity>1 ) {
          std::cout << "[ Pass " << passid << " (" << i << ", X). Skip #" << i << " ]" << std::endl;
        }
        continue;
      }

      for (int j=i+1; j<nendpts; j++) {
        if ( used_endpoints_indices.at(j)==1) continue;
        const BoundarySpacePoint& pts_b = *(spacepts[j]);

        if ( pts_a.type()==pts_b.type() ) {
          if ( m_config.verbosity>1 )
            std::cout << "[ Pass " << passid << ":  endpoints (" << i << "," << j << ") skip same type ]" << std::endl;
          continue; // don't connect same type
        }


        float a2b = 0.;
        for (int v=0; v<3; v++) {
          a2b += (pts_a.pos()[v]-pts_b.pos()[v])*(pts_a.pos()[v]-pts_b.pos()[v]);
        }
        a2b = sqrt(a2b);

        if ( a2b<passcfg.min_point_separation ) {
          if ( m_config.verbosity>1 )
            std::cout << "[ Pass " << passid << ":  endpoints (" << i << "," << j << ") below min separation. ]" << std::endl;
          continue;
        }

        bool use_b = true;
        if ( passcfg.run_radial_filter ) {
          int num_segs_b = 0;
          bool within_line_b = radialfilter.isWithinStraightSegment( pts_b.pos(), img_v, badchimg_v, passcfg.radial_cfg, num_segs_b );
          bool pass_num_segs = passcfg.radial_cfg.min_segments<=num_segs_b && passcfg.radial_cfg.max_segments>=num_segs_b;
          if ( !within_line_b && passcfg.radial_cfg.min_segments<=num_segs_b && passcfg.radial_cfg.max_segments>=num_segs_b )
            use_b = true;
          else
            use_b = false;
        }
        if (!use_b) {
          if ( m_config.verbosity>1 ) {
            std::cout << "[ Pass " << passid << " (" << i << "," << j << "). Skip #" << j << " ]" << std::endl;
          }
          continue;
        }

        if ( m_config.verbosity>1 ) {
          std::cout << "[ Pass " << passid << ": path-finding for endpoints (" << i << "," << j << ") "
                    << "of type (" << pts_a.type() << ") -> (" << pts_b.type() << ") ]" << std::endl;
          std::cout << "  start: (" << pts_a.pos()[0] << "," << pts_a.pos()[1] << "," << pts_a.pos()[2] << ") "
                    << "  end: (" << pts_b.pos()[0] << "," << pts_b.pos()[1] << "," << pts_b.pos()[2] << ") " << std::endl;
          std::cout << "  start: " << pts_a.printImageCoords( img_v.front().meta() ) << " end: " << pts_b.printImageCoords( img_v.front().meta() ) << std::endl;
        }

        // empty candidate
        BMTrackCluster3D track3d;
        LinearTaggerInfo linear_result;
        AStarTaggerInfo astar_result;
        FoxTrotExtenderInfo extender_result;
        bool tracked = false;

        // first run the linear charge tagger, if so chosen
        if ( passcfg.run_linear_tagger ) {
          if ( m_config.verbosity>1 ) {
            std::cout << "  Running linear3dchargetagger." << std::endl;
          }
          BMTrackCluster3D linear_track = runLinearChargeTagger( passcfg, pts_a, pts_b, img_v, badchimg_v, linear_result );
          if ( m_config.verbosity>1 ) {
            std::cout << "  good fraction: " << linear_result.goodfrac << std::endl;
            std::cout << "  majority planes w/ charge fraction: " << linear_result.majfrac << std::endl;
	    if ( linear_track.path3d.size()>0 ) {
	      const std::vector<double>& retstart = linear_track.path3d.front();
	      const std::vector<double>& retend = linear_track.path3d.back();
	      std::cout << "  " << retstart.size() << " " << retend.size() << std::endl;
	      std::cout << "  returned start=(" << retstart[0] << "," << retstart[1] << "," << retstart[2] << ")" << std::endl;
	      std::cout << "  returned end=(" << retend[0] << "," << retend[1] << "," << retend[2] << ")" << std::endl;
	    }
          }
          if ( linear_result.isgood ) {
            if ( m_config.verbosity>1 ) std::cout << "  Linear Result is good. length=" << linear_track.path3d.size() << std::endl;
            std::swap(track3d,linear_track);
          }
          else {
            if ( m_config.verbosity>1 ) std::cout << "  Result is bad." << std::endl;
          }
          tracked = true;
        }//end of if run_linear

        // next run astar tagger
        if ( passcfg.run_astar_tagger
             && (linear_result.goodfrac>passcfg.astar3d_min_goodfrac || linear_result.majfrac>passcfg.astar3d_min_majfrac ) ) {
          if ( m_config.verbosity>1 )
            std::cout << "  Running astar tagger." << std::endl;
          BMTrackCluster3D astar_track = runAStarTagger( passcfg, pts_a, pts_b, img_v, badchimg_v, astar_result );
          tracked = true;
          if ( astar_result.isgood ) {
            if ( m_config.verbosity>1 )  {
              std::vector<double>& astar_end = astar_track.path3d.back();
              std::cout << "  AStar end point: (" << astar_end[0] << "," << astar_end[1] << "," << astar_end[2] << ")" << std::endl;
              std::cout << "  AStar Result is good." << std::endl;
              if ( astar_result.goal_reached )
                std::cout << "  AStar reached goal." << std::endl;
              else
                std::cout << "  AStar did not reach goal." << std::endl;
            }
            std::swap(track3d,astar_track);
          }
          else {
            if ( m_config.verbosity>1 ) std::cout << "  AStar Result is bad." << std::endl;
          }
        }

        if ( tracked )
          num_tracked++;

        // run bezier fitter (doesn't exist yet. write this?)

        // could add criteria to filter here
        bool accept_track = false;
        if ( (passcfg.run_linear_tagger && !passcfg.run_astar_tagger) && linear_result.isgood ) {
          // if we only run the linear tagger, then we save
          accept_track = true;
        }
        else if ( (passcfg.run_linear_tagger && passcfg.run_astar_tagger) && astar_result.isgood && astar_result.goal_reached ) {
          // when we run both, we use the linear tagger as a prefilter for astar, which can be expensive
          accept_track = true;
        }

        if ( accept_track ) {
          if ( m_config.verbosity>1 ) {
            std::cout << "  track found. "
                      << " length: " << track3d.path3d.size()
                      << " empty=" << track3d.isempty()
                      << " linear-good=" << linear_result.isgood
                      << " astar-good=" << astar_result.isgood
                      << " astar-reached=" << astar_result.goal_reached
                      << std::endl;
          }
          // extend the good track
          runFoxTrotExtender( passcfg, track3d.path3d, img_v, badchimg_v, tagged_v, extender_result);
        }

        if ( accept_track && !track3d.isempty() ) {
          std::vector<int> indices(2);
          indices[0] = i;
          indices[1] = j;

	  // Declare a vector for the flash index and producer index at the ith and jth points in the 'flash_idx_v' and 'boundary_type_idx_v' vectors.
	  // 'i' necessarily must be less than j because 'j' begins to iterate at 'i+1'.
	  single_track_endpoint_flash_v.resize(2);
	  single_track_endpoint_flash_v[0]          = spacepts.at( i )->getFlashIndex();
	  single_track_endpoint_flash_v[1]          = spacepts.at( j )->getFlashIndex();
	  track_endpoint_flash_v.emplace_back( std::move(single_track_endpoint_flash_v) );

	  single_track_boundary_type_idx_v.resize(2);
	  single_track_boundary_type_idx_v[0] = spacepts.at( i )->type();
	  single_track_boundary_type_idx_v[1] = spacepts.at( j )->type();
	  track_endpoint_boundary_type_idx_v.emplace_back( std::move(single_track_boundary_type_idx_v) );

	  // Save the information for the boundary points in this vector.
	  single_track_endpoint_v.resize(2);
	  single_track_endpoint_v[0] = *spacepts.at( i );
	  single_track_endpoint_v[1] = *spacepts.at( j );
	  track_endpoint_v.emplace_back( std::move(single_track_endpoint_v) );
	  
          if ( m_config.verbosity>1 ) {
            std::cout << "#### Storing track. size=" << track3d.path3d.size() << ". indices (" << indices[0] << "," << indices[1] << ") ####" << std::endl;
          }
          pass_end_indices.emplace_back( std::move(indices) );
          pass_track_candidates.emplace_back( std::move(track3d) );
	  
        }

      }// second loop over end points
    }// first loop over end points

    if ( m_config.verbosity>0 ) {
      std::cout << "Pass #" << passid << ": "
                << " number of pairs tracked: " << num_tracked
                << " number of tracks passed: " << pass_track_candidates.size()
                << " number of indice pairs: " << pass_end_indices.size()
                << std::endl;
    }


    // post-process tracks (nothing yet)
    for ( auto &track : pass_track_candidates ) {
      trackclusters.emplace_back( std::move(track) );
    }
    for ( auto&indices : pass_end_indices ) {
      used_endpoints_indices.at(indices[0]) = 1;
      used_endpoints_indices.at(indices[1]) = 1;
      
    }

    return;
  }


  larlitecv::BMTrackCluster3D ThruMuTracker::runLinearChargeTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
								    const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
								    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
								    ThruMuTracker::LinearTaggerInfo& result_info ) {

    // linear 3D track
    Linear3DChargeTagger linetrackalgo( pass_cfg.linear3d_cfg ); // finds charge along a line in 3D space

    // get pixel start and end points
    std::vector<int> a_cols(img_v.size(),0);
    std::vector<int> b_cols(img_v.size(),0);
    for (int p=0; p<(int)img_v.size(); p++) {
      a_cols[p] = pts_a.at(p).col;
      b_cols[p] = pts_b.at(p).col;
    }

    PointInfoList straight_track = linetrackalgo.findpath( img_v, badchimg_v, pts_a.at(0).row, pts_b.at(0).row, a_cols, b_cols );
    result_info.numpts   = straight_track.size();
    result_info.goodfrac = straight_track.fractionGood();
    result_info.majfrac  = straight_track.fractionHasChargeOnMajorityOfPlanes();
    if ( result_info.numpts   > pass_cfg.linear3d_min_tracksize
          && (result_info.goodfrac > pass_cfg.linear3d_min_goodfraction
	        || result_info.majfrac  > pass_cfg.linear3d_min_majoritychargefraction) ) {
      result_info.isgood = true;
    }
    else {
      result_info.isgood = false;
    }

    if ( !result_info.isgood ) {
      larlitecv::BMTrackCluster3D emptytrack;
      return emptytrack;
    }

    std::vector< std::vector<double> > path3d;
    for ( auto const& ptinfo : straight_track ) {
      std::vector<double> xyz(3);

      // Declare a float value for the same vector to use that.
      std::vector<float> xyz_float(3,0.0);
      
      for (int i=0; i<3; i++) {
      	xyz[i] = ptinfo.xyz[i];
        xyz_float[i] = ptinfo.xyz[i];
      }

      // Ensure that the 3 coordinates give a value inside the image before appending this value to the track.
      // Convert the 3D position to a 2D pixel image.
      std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( xyz_float, img_v.front().meta(), img_v.size() );

      // Declare an instance of type 'PushBoundarySpacePoint'.
      PushBoundarySpacePoint point_push;

      bool pixel_in_image = point_push.isPixelWithinImage(img_v, imgcoords);

      if (pixel_in_image)
	path3d.emplace_back( std::move(xyz) );

    }

    larlitecv::BMTrackCluster3D track3d( pts_a, pts_b, path3d );
    return track3d;
  }

  larlitecv::BMTrackCluster3D ThruMuTracker::runAStarTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                  const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                  ThruMuTracker::AStarTaggerInfo& result_info ) {

    // collect meta/translate start/goal tick/wires to row/col in compressed image
    std::vector< const larcv::ImageMeta* > meta_compressed_v;
    std::vector<int> start_cols(m_img_compressed_v.size(),0);
    std::vector<int> start_rows(m_img_compressed_v.size(),0);
    std::vector<int> goal_cols(m_img_compressed_v.size(),0);
    std::vector<int> goal_rows(m_img_compressed_v.size(),0);
    std::vector< const larcv::Pixel2D* > start_pix;
    std::vector< const larcv::Pixel2D* > goal_pix;


    for (size_t p=0; p<m_img_compressed_v.size(); p++) {

      // get start/end point informatio in compressed image
      const larcv::ImageMeta* ptr_meta = &(m_img_compressed_v[p].meta()); // compressed image meta
      const larcv::ImageMeta& meta = img_v[p].meta(); // original image meta

      const larlitecv::BoundaryEndPt& start_endpt = pts_a[p];
      start_rows[p] =  ptr_meta->row( meta.pos_y( start_endpt.row ) );
      start_cols[p] =  ptr_meta->col( meta.pos_x( start_endpt.col ) );
      start_pix.push_back( new larcv::Pixel2D( start_endpt.getcol(), start_endpt.getrow() ) );

      const larlitecv::BoundaryEndPt& goal_endpt  = pts_b[p];
      goal_rows[p]  =  ptr_meta->row( meta.pos_y( goal_endpt.row ) );
      goal_cols[p]  =  ptr_meta->col( meta.pos_x( goal_endpt.col ) );
      goal_pix.push_back( new larcv::Pixel2D( goal_endpt.getcol(), goal_endpt.getrow() ) );

      meta_compressed_v.push_back( ptr_meta );
    }

    larcv::AStar3DAlgo algo( pass_cfg.astar3d_cfg );
    std::vector<larcv::AStar3DNode> path;
    result_info.goal_reached = false;
    result_info.isgood = true;
    int goalhit = 0;
    try {
      path = algo.findpath( m_img_compressed_v, m_badch_compressed_v, m_badch_compressed_v, // tagged_compressed_v
			    start_rows.front(), goal_rows.front(), start_cols, goal_cols, goalhit );

      if ( goalhit==1 ) {
        result_info.goal_reached = true;
      }
    }
    catch (const std::exception& e) {
      std::cout << "*** [ Exception running astar3dalgo::findpath: " << e.what() << "] ***" << std::endl;
      result_info.isgood = false;
      result_info.goal_reached = false;
      result_info.nbad_nodes = -1;
      result_info.total_nodes = -1;
      std::vector< std::vector<double> > empty_path3d;
      BMTrackCluster3D track3d( pts_a, pts_b, empty_path3d );
      for (int p=0; p<3; p++ ) {
        delete start_pix.at(p);
        delete goal_pix.at(p);
      }
      return track3d; // return empty track
    }

    result_info.nbad_nodes = 0;
    result_info.total_nodes = 0;
    for ( auto& node : path ) {
      if ( node.badchnode )
        result_info.nbad_nodes+=1.0;
      result_info.total_nodes+=1.0;
    }

    // if majority are bad ch nodes, reject this track
    result_info.frac_bad = float(result_info.nbad_nodes)/float(result_info.total_nodes);
    if ( result_info.frac_bad>0.5 || result_info.total_nodes<=3)
      result_info.isgood = false;

    BMTrackCluster3D track3d = AStarNodes2BMTrackCluster3D( path, img_v, pts_a, pts_b, 0.3 );

    for (int p=0; p<3; p++ ) {
      delete start_pix.at(p);
      delete goal_pix.at(p);
    }

    return track3d;
  }

  void ThruMuTracker::runFoxTrotExtender( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg, std::vector<std::vector<double> >& track,
                  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
                  ThruMuTracker::FoxTrotExtenderInfo& result_info ) {
    // try to extend the track.
    // make sure we don't go backwards

    if ( !pass_cfg.run_foxtrot_extender || track.size()<2 ) {
      result_info.isgood = false;
      return;
    }

    ThruMuFoxExtender extender_algo( pass_cfg.foxextend_cfg );

    // need a forward and backward extension
    result_info.isgood  = extender_algo.extendTrack( track, img_v, badch_v, tagged_v );

    return;
  }

  // Define a function that will use 'GeneralFlashMatchAlgo' to flashmatch a track with a reconstructed flash of light in the event.
  // Input parameters:
  // flash_match_larlite_pset: This is the configuration PSet for the larlite flash matching infrastrcture that is used.
  // flash_match_config: This is the configuration file for the 'GeneralFlashMatchAlgo' function that will perform the functionality needed at this point in the program.
  // img_v: The input wireplane images for the events that we are considering.
  // tagged_v: The images of the tracks that have already been tagged by the 'ThruMu' tracker.
  // spacepts: All possible boundary points that can be used to make tracks in ThruMu.
  // opflash_v: The vector of flashes from the event, both of type 'simpleFlashBeam' and 'simpleFlashCosmic'.
  // trackclusters: The 'BMTrackCluster3D' objects that have already been generated by ThruMu.
  // impossible_match_endpoint_idx_v: A vector of integers for the flashes that have found to not be able to generate a track.
  // already_matched_flash_idx_v: This a vector of the flashes that are already determined (in a previous pass) to have a good
  //                              match to another track and have therefore been removed from consideration for other tracks.
  // num_of_tracks_added_in_pass: The number of tracks added in the pass for which we are flashmatching information.
  // track_endpoint_flash_v: This vector contains information for which flash (if any) determined the endpoint for a track.  This is now a vector of vectors.
  // track_endpoint_boundary_type_idx_v: This vector contains information for the flash producer used to generate the flash that produced the endpoint.  This is now a vector of vectors. 
  // track_endpoint_v: The vector of track endpoints from the entire event.
  // opflashsets: This is a vector of the pointers to the flashes in the event.
  // anode_and_cathode_only: This only compares the endpoints for anode-piercing/cathode-piercing tracks, for which one flash is matched to the track.
  // single_pass: This says if we are flashmatching the entire event or if we are just considering one pass of the tracker ('true' if we are considering one pass and 'false' if we are considering the entire event.
  void ThruMuTracker::flashMatchTracks( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< const BoundarySpacePoint* >& spacepts,
					const std::vector< larlite::event_opflash* >& opflash_v, std::vector< BMTrackCluster3D >& trackclusters,
					std::vector< std::vector< BoundarySpacePoint > >& impossible_match_endpoint_v, std::vector< BoundaryFlashIndex >& already_matched_flash_v,
					std::vector< int >& well_matched_tracks_idx_v, const int& num_of_tracks_added_in_pass,
					std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v, std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v,
					std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v, bool anode_and_cathode_only, bool single_pass ) {

    //std::cout << "At the start of the 'flashMatchTracks' function." << std::endl;

    // Declare an object of type 'GeneralFlashMatchAlgo' using the configuration object from 'GeneralFlashMatchAlgoConfig'.
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Declare a flash manager for use in this function.
    // Just a test at first to see if these will work......
    ::flashana::FlashMatchManager       _mgr;
    // Configure this with the flashmatch manager available.
    _mgr.Configure( flash_match_config.m_flashmatch_config );
    _mgr.Reset();

    // Declare three vectors: one for the BMTrackCluster3D objects, one for the flash indices of the track endpoints (two per track),
    //  and one for the flash producer of the track endpoints (two per track).
    //  Clear each vector to ensure that there are no issues with space in the vector.
    std::vector< larlitecv::BMTrackCluster3D > trackclusters_from_pass;
    trackclusters_from_pass.clear();
    std::vector< std::vector< BoundarySpacePoint > > track_endpoint_v_from_pass;
    track_endpoint_v_from_pass.clear();
    std::vector< std::vector< BoundaryFlashIndex > > track_endpoint_flash_v_from_pass;
    track_endpoint_flash_v_from_pass.clear();
    std::vector< std::vector< BoundaryEnd_t > > track_endpoint_boundary_type_idx_v_from_pass;
    track_endpoint_boundary_type_idx_v_from_pass.clear();

    // Take a given track added in the last pass and compare it to all of the flashes in the event, finding one flash closer to it than any of the others.
    for ( int track_last_pass_iter = int( trackclusters.size() - num_of_tracks_added_in_pass); track_last_pass_iter < int( trackclusters.size() ); ++track_last_pass_iter ) {

      trackclusters_from_pass.push_back( trackclusters.at( track_last_pass_iter ) );
      track_endpoint_v_from_pass.push_back( track_endpoint_v.at( track_last_pass_iter ) );
      track_endpoint_flash_v_from_pass.push_back( track_endpoint_flash_v.at( track_last_pass_iter ) );                                                                                    
      track_endpoint_boundary_type_idx_v_from_pass.push_back( track_endpoint_boundary_type_idx_v.at( track_last_pass_iter ) );                                                                 
    }

    // Print out the length of 'well_matched_tracks_idx_v'.
    //std::cout << "The length of 'trackclusters_from_pass' = " << trackclusters_from_pass.size() << "." << std::endl;
    //std::cout << "The length of 'well_matched_track_idx_v' at the start of the pass = " << well_matched_tracks_idx_v.size() << "." << std::endl;

    // Turn this into a vector of larlite tracks using the the 'GeneralFlashMatchAlgo' functionality.
    std::vector < larlite::track > larlite_track_vector = flash_match_obj.generate_tracks_between_passes( trackclusters_from_pass );

    // Generate a single opflash vector using the functionality in GeneralFlashMatchAlgo.
    std::vector< BoundaryFlashIndex > boundary_flash_index_vector = flash_match_obj.generate_boundaryflashindex_vector_for_event( opflash_v );

    // Print out the number of flashes in this event.
    std::cout << "The number of opflash objects recorded in this event = " << boundary_flash_index_vector.size() << "." << std::endl;

    int num_beam_window_flashes = 0;

    // Loop through the vector and count how many objects have 'ivec' equal to 0.
    for ( size_t i = 0; i < boundary_flash_index_vector.size(); ++i ) {

      if ( boundary_flash_index_vector.at( i ).ivec == 0 ) {
	++num_beam_window_flashes;
      }

    }

    std::cout << "The number of flashes reconstructed in the beam window = " << num_beam_window_flashes << "." << std::endl;

    // Declare an 'if....else' loop for if we are considering a single pass or if we are considering the entire event.
    if ( single_pass ) {

      // Resize 'well_matched_tracks_idx_v' to the size of 'trackclusters_from_pass', which is the same size as 'larlite_track_vector'.
      // We will fill the entries with -1, meaning that no decision has (yet) been rendered on these tracks.
      well_matched_tracks_idx_v.resize( trackclusters_from_pass.size(), -1);

      // Check the boolean 'anode_and_cathode_only'.  If 'true', then we will only consider the anode-piercing/cathode-piercing tracks, which are determined by a single flash.
      if ( anode_and_cathode_only ) {
	flashMatchAC( flash_match_config, larlite_track_vector, trackclusters, track_endpoint_v, track_endpoint_flash_v, track_endpoint_boundary_type_idx_v, boundary_flash_index_vector, 
		      track_endpoint_flash_v_from_pass, track_endpoint_boundary_type_idx_v_from_pass, impossible_match_endpoint_v, track_endpoint_v_from_pass, well_matched_tracks_idx_v, 
		      already_matched_flash_v );
      }
      else {
	
	bool entire_event = false;
	flashMatchAC( flash_match_config, larlite_track_vector, trackclusters, track_endpoint_v, track_endpoint_flash_v, track_endpoint_boundary_type_idx_v, boundary_flash_index_vector, 
		      track_endpoint_flash_v_from_pass, track_endpoint_boundary_type_idx_v_from_pass, impossible_match_endpoint_v, track_endpoint_v_from_pass, well_matched_tracks_idx_v, 
		      already_matched_flash_v );
	
	// Use 'track_endpoint_v_from_pass' here instead of 'track_endpoint_v'.  They will only be used to 
	flashMatchYZFaceTracks( flash_match_config, trackclusters, larlite_track_vector, boundary_flash_index_vector, impossible_match_endpoint_v,
			      track_endpoint_v_from_pass, well_matched_tracks_idx_v, already_matched_flash_v, single_pass );


      }
    
    }

    else {
      
      //std::cout << "Entering the loop where we flash-match the YZ tracks from the full event." << std::endl;
      //std::cout << "The size of 'trackclusters' in this function, the 'flashMatchTracks' function,  = " << trackclusters.size() << "." << std::endl;

      // Resize 'well_matched_tracks_idx_v' to be the size of 'trackclusters'.  'larlite_track_vector' will be reset within the 'flashMatchYZFaceTracks' function.
      well_matched_tracks_idx_v.resize( trackclusters.size(), -1);

      // 'single_pass' is set to 'false' when we use the 'flashMatchYZFaceTracks' here.
      flashMatchYZFaceTracks( flash_match_config, trackclusters, larlite_track_vector, boundary_flash_index_vector, impossible_match_endpoint_v, 
			      track_endpoint_v, well_matched_tracks_idx_v, already_matched_flash_v, single_pass );

    }


    // Put the code to remove the impossible endpoints here.
    // Repeat the previous loop to look for any tracks that may have been determined by a set of impossible endpoints.                                                                                 
    // This is only for tracks that passed the selection cuts.                                                                                                                                         
    for ( size_t track_endpt_iter = 0; track_endpt_iter < track_endpoint_v_from_pass.size(); ++track_endpt_iter ) {

      // Declare a boolean for if the track has impossible endpoints or not.                                                                                                                            
      bool determined_by_set_of_impossible_match_endpoints = false;

      // Declare values for the starting and ending endpoints of the track.
      std::vector< float > first_track_endpt_pos(3,0.0);  
      std::vector< float > second_track_endpt_pos(3,0.0); 
      
      // Set these variables.
      track_endpoint_v_from_pass.at( track_endpt_iter )[0].pos( first_track_endpt_pos );
      track_endpoint_v_from_pass.at( track_endpt_iter )[1].pos( second_track_endpt_pos );

      // Start an inner loop over the sets of impossible match endpoints.                                                                                                                              
      for ( size_t impossible_match_iter = 0; impossible_match_iter < impossible_match_endpoint_v.size(); ++impossible_match_iter ) {

	// Declare objects for the first and second impossible match endpoints.
	std::vector< float > first_impossible_endpt_pos(3,0.0);  
	std::vector< float > second_impossible_endpt_pos(3,0.0); 

	impossible_match_endpoint_v.at( impossible_match_iter )[0].pos( first_impossible_endpt_pos );
	impossible_match_endpoint_v.at( impossible_match_iter )[1].pos( second_impossible_endpt_pos );

	// Loop through the endpoints in 'impossible_match_endpoint_v' to see if any of them match the information in 'impossible_match.  Use the endpoint information here.                           
	// Same point in index.               
	if ( fabs( first_track_endpt_pos[0] - first_impossible_endpt_pos[0] ) < 0.00001 && fabs( first_track_endpt_pos[1] - first_impossible_endpt_pos[1] ) < 0.00001 && fabs( first_track_endpt_pos[2] - first_impossible_endpt_pos[2] ) < 0.00001 && fabs( second_track_endpt_pos[0] - second_impossible_endpt_pos[0] ) < 0.00001 && fabs( second_track_endpt_pos[1] - second_impossible_endpt_pos[1] ) < 0.00001 && fabs( second_track_endpt_pos[2] - second_impossible_endpt_pos[2] ) < 0.00001 ) {
	  determined_by_set_of_impossible_match_endpoints = true;
	  break;
	}

	// Opposite point in index.                                                                                                                                                                   
	if ( fabs( second_track_endpt_pos[0] - first_impossible_endpt_pos[0] ) < 0.00001 && fabs( second_track_endpt_pos[1] - first_impossible_endpt_pos[1] ) < 0.00001 && fabs( second_track_endpt_pos[2] - first_impossible_endpt_pos[2] ) < 0.00001 && fabs( first_track_endpt_pos[0] - second_impossible_endpt_pos[0] ) < 0.00001 && fabs( first_track_endpt_pos[1] - second_impossible_endpt_pos[1] ) < 0.00001 && fabs( first_track_endpt_pos[2] - second_impossible_endpt_pos[2] ) < 0.00001 ){
	  determined_by_set_of_impossible_match_endpoints = true;
	  break;
	}

	// Set the value of 'well_match_tracks_idx_v' at 'track_endpt_iter' equal to 0 if the conditional is true.                                                                                     
	// This is contingent on 'track_endpoint_v_from_pass' having the same length as 'track_endpt_iter', which it has been shown to have.                                                            
	if ( determined_by_set_of_impossible_match_endpoints == true ) {
	  well_matched_tracks_idx_v.at( track_endpt_iter ) = 0;
	}

      }

    }

    return;

  }

  // Declare a function for matching the anode-piercing/cathode-piercing tracks to their corresponding flash.
  // Inputs:
  // 'flash_match_config': This object configures the 'GeneralFlashMatchAlgo' object which is used to perform the flash-matching functionality in this class.
  // 'larlite_track_vector' - This is the vector of larlite tracks that corresponds to the entire vector of reconstructed tracks from the event.
  // trackclusters: The 'BMTrackCluster3D' objects that have already been generated by ThruMu.
  // track_endpoint_v: The vector of track endpoints from the entire event. 
  // track_endpoint_flash_v: This vector contains information for which flash (if any) determined the endpoint for a track.  This is now a vector of vectors.                                            
  // track_endpoint_boundary_type_idx_v: This vector contains information for the flash producer used to generate the flash that produced the endpoint.  This is now a vector of vectors.                
  // 'boundary_flash_index_vector' - This is the vector of the flash information in the event stored in the BoundaryFlashIndex class.
  // 'track_endpoint_flash_v_from_pass': This is a vector of the flashes that correspond to the endpoints that make up the tracks in this pass of the tagger.
  // 'track_endpoint_boundary_type_idx_v': This is a vector of the indices of the flash producer that determined the endpoint of the track.
  //                                       This is created with the same scheme that the 'single_opflash_vector' was created.
  //                                       This will be used to place a cut on either anode-piercing tracks or cathode-piercing tracks.
  // 'impossible_match_endpoints': This is a set of the endpoints that cannot form a valid track based on the information shown here.
  //                               This will be added to based on the information in this function.
  // 'already_matched_flash_v': This is a set of 'BoundaryFlashIndex' objects that correspond to the flashes that have already been well-matched to an anode-piercing/cathode-piercing track.
  void ThruMuTracker::flashMatchAC( GeneralFlashMatchAlgoConfig& flash_match_config, std::vector< larlite::track >& larlite_track_vector, std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
		     std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v, std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v,
		     std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v,  const std::vector< BoundaryFlashIndex >& boundary_flash_index_vector,
		     const std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v_from_pass,
		     const std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v_from_pass,
		     std::vector< std::vector< BoundarySpacePoint > >& impossible_match_endpoint_v, const std::vector< std::vector< BoundarySpacePoint > > & track_endpoint_v_from_pass,
		     std::vector < int >& well_matched_tracks_idx_v, std::vector< BoundaryFlashIndex >&  already_matched_flash_v ) 
  {

    //std::cout << "At the start of the 'flashMatchAC' function." << std::endl;
    
    // Create a new object of class 'GeneralFlashMatchAlgo'.
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Include the size of the vector of the track endpoints and the flash information for each of those track endpoints.
    //std::cout << "Size of the vector of track endpoints = " << track_endpoint_v_from_pass.size() << "." << std::endl;
    //std::cout << "Size of the vector of flash info corresponding to the track endpoints = " << track_endpoint_flash_v_from_pass.size() << "." << std::endl;

    //std::cout << "Size of the 'larlite_track_vector' in the 'flashMatchAC' function = " << larlite_track_vector.size() << "." << std::endl;
 
    // Loop through the qclusters added in the pass.
    //   If both of the entries in 'track_endpoint_flash_idx_v_from_pass' are set to -1,
    //   then that means that the track was not made from an anode-piercing or cathode-piercing endpoint.
    // The tracks were already checked for impossible match endpoints above in 'flashMatchTracks'.
    for ( size_t qcluster_iter = 0; qcluster_iter < larlite_track_vector.size(); ++qcluster_iter ) {  // Keep the name of the iterator as 'qcluster_iter' even though we are looking at larlite tracks.

      // Place a print statement telling where you are within this loop.
      //std::cout << "At the #" << qcluster_iter << " point iterating through 'larlite_track_vector'." << std::endl;

      //std::cout << "First Point Type: " << track_endpoint_v_from_pass.at( qcluster_iter)[0].type() << "." << std::endl;
      //std::cout << "Second Point Type: " << track_endpoint_v_from_pass.at( qcluster_iter)[1].type() << "." << std::endl;


      // If the track passes that selection, then it is anode-piercing or cathode-piercing.
      //   Declare an object for the opflash that determined the track endpoints based on which index in
      //   'track_endpoint_flash_idx_v' is greater than 0.
      // I could have also gotten the track producer information directly from the list of opflashes........an improvement later.
      BoundaryFlashIndex boundary_flash_idx_for_endpoint;

      bool               loop_entered = false;

      // The first track endpoint is A/C and corresponds to a flash.
      if ( track_endpoint_flash_v_from_pass.at( qcluster_iter )[0].ivec > -0.001 && track_endpoint_v_from_pass.at( qcluster_iter)[0].type() != larlitecv::kImageEnd ) {
	//std::cout << "Saving the 'BoundaryFlashIndex' object from the first endpoint." << std::endl;
	boundary_flash_idx_for_endpoint = track_endpoint_flash_v_from_pass.at( qcluster_iter )[0];

	//std::cout << "Type of First Boundary Point = " << track_endpoint_v_from_pass.at( qcluster_iter)[0].type() << "." << std::endl;

	loop_entered = true;
      }

      // The second track endpoint is A/C and corresponds to a flash.
      if ( track_endpoint_flash_v_from_pass.at( qcluster_iter )[1].ivec > -0.001 && track_endpoint_v_from_pass.at( qcluster_iter)[1].type() != larlitecv::kImageEnd ) {
	//std::cout << "Saving the 'BoundaryFlashIndex' object from the second endpoint." << std::endl;
	boundary_flash_idx_for_endpoint = track_endpoint_flash_v_from_pass.at( qcluster_iter )[1];

	//std::cout << "Type of Second Boundary Point = " << track_endpoint_v_from_pass.at( qcluster_iter)[1].type() << "." << std::endl;

	loop_entered = true;
      } 

      // Continue if 'loop_entered' is equal to false.
      if ( loop_entered == false ) {
	//std::cout << "Neither of the points are A/C piercing so I am not rendering a decision and continuing!" << std::endl;
	well_matched_tracks_idx_v.at( qcluster_iter ) = -1;
	continue;
      }
      
      // Check to see if that flash that defines these endpoints have already been well-matched to a track.
      bool well_matched_flash = false;

      //std::cout << "Size of 'already_matched_flash_v' = " << already_matched_flash_v.size() << "." << std::endl;

      // Loop through all of 'BoundaryFlashIndex' entries in 'already_matched_flash_v', comparing the 'ivec' and 'idx' values of each of the entries to those for 'boundary_flash_idx_for_endpoint'.
      for ( size_t already_matched_flash_iter = 0; already_matched_flash_iter < already_matched_flash_v.size(); ++already_matched_flash_iter ) {

	if ( boundary_flash_idx_for_endpoint.idx == already_matched_flash_v.at( already_matched_flash_iter ).idx && boundary_flash_idx_for_endpoint.ivec == already_matched_flash_v.at(already_matched_flash_iter ).ivec ) {
	  well_matched_flash = true;
	  break;
	}

      }

      // If this track is matched to a well-matched flash, then set its entry in 'well_matched_tracks_idx_v' equal to 0 and continue.
      if ( well_matched_flash == true ) {

	//std::cout << "Well-matched flash equals true so I am performing more functions on this track!" << std::endl;

	bool matched_to_other_AC_track = compareToOtherACTrackMatchedToThisFlash( flash_match_config, larlite_track_vector.at( qcluster_iter ), boundary_flash_idx_for_endpoint, trackclusters, track_endpoint_v, track_endpoint_flash_v, track_endpoint_boundary_type_idx_v );
	
	// The case in which the other track that the flash is matched to is not A/C piercing - too difficult to compare them in this prototype of the algorithm - or in which the other track has a better match than this track does.
	if ( matched_to_other_AC_track == false )  {
	  well_matched_tracks_idx_v.at( qcluster_iter ) = 0;
	  continue;
	}
      }

      // Move this up just to see that it works.
      //std::cout << "About to save the flash producer of the endpoint." << std::endl;
      int opflash_producer_idx           = boundary_flash_idx_for_endpoint.ivec;
      //std::cout << "Done saving the flash producer of the endpoint." << std::endl;

      //if ( boundary_flash_idx_for_endpoint.popflash == NULL ) {
      //	std::cout << "This flash pointer is NULL." << std::endl;
      //}

      // Print the information for the indices to the opflash pointer.
      //std::cout << "boundary_flash_idx_for_endpoint.ivec = " << boundary_flash_idx_for_endpoint.ivec << "." << std::endl;
      //std::cout << "boundary_flash_idx_for_endpoint.idx = " << boundary_flash_idx_for_endpoint.idx << "." << std::endl;

      // Try to make an opflash pointer and then obtain the flash object from that pointer.
      // First just try a 'const' modifier on the 'opflash_object' at this point in the function.
      const larlite::opflash* opflash_pointer = boundary_flash_idx_for_endpoint.popflash; // See if this works, and then take the 'event_opflash' info out of the arguments for the function.

      const larlite::opflash  opflash_object  = *opflash_pointer;

      // Declare a new larlite track object that is t0-tagged with the flash that determined its A/C endpoints.
      larlite::track t0_tagged_track;
      t0_tagged_track.clear_data();

      // Loop through the entries of the existing larlite track object and fill it with the track points, but t0-tagged.
      // This will allow the extensions to take into account the fact that the track should be t0-tagged.
      for ( size_t point_iter = 0; point_iter < larlite_track_vector.at( qcluster_iter ).NumberTrajectoryPoints(); ++point_iter ) {

	// Declare a TVector object with the trajectory points (t0-tagged) of this entry in 'larlite_track_vector' at this point in its trajectory.
	TVector3 larlite_trk_pt( larlite_track_vector.at( qcluster_iter ).LocationAtPoint( point_iter ).X() - opflash_object.Time()*0.1114, larlite_track_vector.at( qcluster_iter ).LocationAtPoint( point_iter ).Y(), larlite_track_vector.at( qcluster_iter ).LocationAtPoint( point_iter ).Z() );

	t0_tagged_track.add_vertex( larlite_trk_pt );
	t0_tagged_track.add_direction( larlite_track_vector.at( qcluster_iter ).DirectionAtPoint( point_iter ) );

      }

      //std::cout << "The number of points on the larlite t0-tagged track = " << t0_tagged_track.NumberTrajectoryPoints() << "." << std::endl;
      
      // Declare a qcluster object that will be extended while taking the t0-tagged x-coordinate into account.
      flashana::QCluster_t t0_tagged_qcluster;

      // Set the 'time' attribute of the qcluster equal to the time of this opflash.
      t0_tagged_qcluster.time = opflash_object.Time();
      
      // Use the 'FlashMatchInterface' functionality to expand this qcluster near the boundary.
      flash_match_obj.ExpandQClusterNearBoundaryFromLarliteTrack( t0_tagged_qcluster, t0_tagged_track, 10000.0, 10.0 );

      // Print out the number points on the qcluster.
      //std::cout << "The number of points on the qcluster = " << t0_tagged_qcluster.size() << "." << std::endl;

      // Determine the chi2 match between the qcluster and the opflash object.

      // Declare variables for the total PEs of the data flash and the total PEs of the hypothesis flash.
      float totpe_data = 0.0;
      float totpe_hypo = 0.0;

      float chi2 = flash_match_obj.generate_chi2_in_track_flash_comparison( t0_tagged_qcluster, opflash_object, totpe_data, totpe_hypo, opflash_producer_idx );

      std::cout << "The chi2 for the match of track #" << qcluster_iter << " among the tracks added in this pass to the qcluster vector = " << chi2 << "." << std::endl;

      // Put a chi2 cut on the match between the flash and the qcluster.
      if ( chi2 < flash_match_config.chi2_anode_cathode_cut ) { // to be defined somewhere... 

	//std::cout << "Track #" << qcluster_iter << " is A/C piercing and has a chi2 value below our cut and will be kept." << std::endl;
	//std::cout << "y-coordinate of the first track point = " << t0_tagged_track.LocationAtPoint( 0 )[1] << "." << std::endl;

	// Print out the cut that these tracks have to satisfy.
	//std::cout << "The chi2 cut for anode/cathode piercing tracks = " << flash_match_config.chi2_anode_cathode_cut << "." << std::endl;

	std::cout << "The chi2 value of the track passing the cuts in the 'flashMatchACTracks' function = " << chi2 << "." << std::endl;
	std::cout << "This track is matched to flash with ivec = " << boundary_flash_idx_for_endpoint.ivec << " and idx = " << boundary_flash_idx_for_endpoint.idx << "." << std::endl;
	std::cout << "The length of 'already_matched_flash_v' before this flash has been compared = " << already_matched_flash_v.size() << "." << std::endl;

	// If the flash has not yet been compared, then you should add it to 'already_matched_flash_v'.
	if ( well_matched_flash == false ) {
	  // Add to the end of this vector the flash object corresponding to the qcluster, 'boundary_flash_idx_for_endpoint', using the 'insert' functionality.
	  already_matched_flash_v.push_back( boundary_flash_idx_for_endpoint );
	}
	
	// Change the index of 'well_matched_tracks_idx_v' at this index to equal 1, because the track is well-reconstructed according to the chi2 cut.
	well_matched_tracks_idx_v.at( qcluster_iter ) = 1 ;

      }

      else { // If the track's chi2 was greater than the anode/cathode piercing chi2 cut value, then put the track's endpoints in the 'impossible_match_endpoint_idx_v' list.

	//std::cout << "Track #" << qcluster_iter << " is A/C piercing and has a chi2 iter above our cut and will be removed." << std::endl;
	//std::cout << "y-coordinate of the first track point = " << t0_tagged_track.LocationAtPoint( 0 )[1] << "." << std::endl;
	
	impossible_match_endpoint_v.push_back( track_endpoint_v_from_pass.at( qcluster_iter ) );

	// Change the index of 'well_matched_tracks_idx_v' at this index to equal 0, because the track is not well-reconstructed according to the chi2 cut.
	well_matched_tracks_idx_v.at( qcluster_iter ) = 0;

      }
	

    }
	
    return;

  }

  // Declare a function that will loop through the non-anode-piercing/cathode-piercing tracks to match them to the correct flash.
  // Input: 'flash_match_config' - This is the 'GeneralFlashMatchAlgoConfig' object that will be used to initialize the 'GeneralFlashMatchAlgo' object.
  //        'trackclusters'      - This is the total vector of remaining track clusters if we would like to look at the tracks in the entire event instead of just those 
  //        'larlite_track_vector' - This is the vector of larlite tracks that corresponds to the entire vector of reconstructed tracks from the event.
  //        'boundary_flash_index_vector' - This is the vector of the flash information in the event stored in the 'struct' BoundaryFlashIndex.
  //        'impossible_match_endpoint_idx_v' - This is a vector of vectors that contains the endpoints that are an impossible match with one another.
  //        'track_endpoint_v_pass_or_full'   - This is a vector of endpoints, either corresponding to the endpoints of all the tracks or only those reconstructed in the last pass.
  //        'well_matched_tracks_idx' - This vector contains the information from the tracks in the event that are well-matched already to a single flash.
  //        'already_matched_flash_idx_v' - This vector contains the indices of the flashes that are already well-matched to a track in the event, meaning that they can be passed over.
  //        'single_pass' - This boolean asks if we are looking at all of the tracks in the event or just at ones from a single pass (if this metric evaluates to 'true').
  void ThruMuTracker::flashMatchYZFaceTracks( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
					      std::vector < larlite::track >& larlite_track_vector, const std::vector< BoundaryFlashIndex >& boundary_flash_index_vector,
					      std::vector< std::vector< BoundarySpacePoint > >& impossible_match_endpoint_v,
					      const std::vector< std::vector< BoundarySpacePoint > > & track_endpoint_v_pass_or_full,
					      std::vector< int >& well_matched_tracks_idx_v, std::vector< BoundaryFlashIndex >& already_matched_flash_v, bool single_pass ) {

    // Declare a 'GeneralFlashMatchAlgo' object with the 'flash_match_config' object passed to 'flashMatchAllTracks'.
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Use the flag 'entire_event' to determine what the scope of flashes that we are looking at is.  This will allow for the definition of a vector of larlite tracks that we are interested in looking at within this event.
    // If '!entire_event', then we will only look at those qclusters in the input 'qcluster_vector', which only contains qclusters from this pass.
    // If 'entire_event', then we will look at those qclusters from the entire event, finding some flash that they are compatible with in the event.
    std::vector< larlite::track > larlite_tracks_being_checked;
    larlite_tracks_being_checked.clear();

    // Declare vectors of the best chi2 for each flash and the corresponding qcluster index for the best qcluster that it is matched to.
    std::vector< float > best_chi2_v;
    best_chi2_v.clear();
    std::vector< int > corresponding_qcluster_idx_v;
    corresponding_qcluster_idx_v.clear();

    if ( !single_pass ) {

      //std::cout << "At the start of the 'flashMatchYZFaceTracks' function while looping over the entire event." << std::endl;
      //std::cout << "Number of tracks that we are considering in total = " << trackclusters.size() << "." << std::endl;
      
      // Convert 'trackclusters' (the entire vector) to a vector of larlite tracks.
      larlite_tracks_being_checked = flash_match_obj.generate_tracks_between_passes( trackclusters );

      //std::cout << "Output of the 'generate_tracks_between_passes' function = " << larlite_tracks_being_checked.size() << "." << std::endl;

      // Spend more time on this later......working right now on getting the program to work for flashmatching just in two passes.                                                                        
      well_matched_tracks_idx_v.clear();
      well_matched_tracks_idx_v.resize( larlite_tracks_being_checked.size(), -1);

      int AC_tracks_already_passed_cuts = 0;

      // Loop through 'track_endpoint_v_pass_or_full' and fill the 'well_matched_track_idx_v' with a '1' if either of the tracks endpoints are A/C piercing.  That means that they've already been checked and they correspond to a well-matched track.
      for ( size_t track_endpoint_iter = 0; track_endpoint_iter < track_endpoint_v_pass_or_full.size(); ++track_endpoint_iter ) {

	for ( size_t pt_iter = 0; pt_iter < 2; ++pt_iter ) {

	  if ( track_endpoint_v_pass_or_full.at( track_endpoint_iter )[ pt_iter ].type() == larlitecv::kAnode || track_endpoint_v_pass_or_full.at( track_endpoint_iter )[ pt_iter ].type() == larlitecv::kCathode ) {
	    // These tracks will not be ranked because they are already set to '1', meaning that they are well-matched tracks.
	    well_matched_tracks_idx_v.at( track_endpoint_iter ) = 1;
	    ++AC_tracks_already_passed_cuts;
	    
	  }

	}

      }

      //std::cout << "The number of 'AC_tracks_already_passed_the_cuts' = " << AC_tracks_already_passed_cuts << "." << std::endl;

    }

    else {

	larlite_tracks_being_checked = larlite_track_vector;
    }

    //std::cout << "The size of 'larlite_tracks_already_checked' after the 'else' loop = " << larlite_tracks_being_checked.size() << "." << std::endl;

    // Declare an empty vector of vector of 'ThruMuTracker::FlashMatchCandidate' candidates.
    std::vector< std::vector< ThruMuTracker::FlashMatchCandidate > > all_opflash_candidate_list;
    all_opflash_candidate_list.clear();

    // Declare a vector for the flashes that are not well-matched already, for which a vector of tracks has to be compared.
    std::vector< BoundaryFlashIndex > opflash_with_track_lists_idx_v;
    opflash_with_track_lists_idx_v.clear();

    // Print statement to make sure that this is working the way that I think it is.
    //std::cout << "Starting to loop over the flashes to see which tracks should be matched to which flashes in general." << std::endl;

    // Loop through the eligible flashes, using the 'struct' declared above to create a 'list' ranking the tracks that are matched to each flash in an event.
    // The list of flashes is now contained in 'boundary_flash_index_vector'.
    for ( size_t opflash_i = 0; opflash_i < boundary_flash_index_vector.size(); ++opflash_i ) {
      
      // Check to see if that flash that defines these endpoints have already been well-matched to a track.                                                                                            
      bool well_matched_flash = false;

      // Loop through all of 'BoundaryFlashIndex' entries in 'already_matched_flash_v', comparing the 'ivec' and 'idx' values of each of the entries to those for 'boundary_flash_idx_for_endpoint'.       
      for ( size_t already_matched_flash_iter = 0; already_matched_flash_iter < already_matched_flash_v.size(); ++already_matched_flash_iter ) {

        if ( boundary_flash_index_vector.at( opflash_i).idx == already_matched_flash_v.at( already_matched_flash_iter ).idx && boundary_flash_index_vector.at( opflash_i).ivec == already_matched_flash_v.at(already_matched_flash_iter ).ivec ) {
          well_matched_flash = true;
          break;
        }

      }

      // If this flash is already well-matched, then continue on in the loop.
      if ( well_matched_flash == true ) {
	continue;
      }

      std::vector< ThruMuTracker::FlashMatchCandidate > opflash_track_match_list;
      opflash_track_match_list.clear();

      std::vector< ThruMuTracker::FlashMatchCandidate > ordered_opflash_track_match_list;
      ordered_opflash_track_match_list.clear();

      // Declare a variable for the opflash that you are comparing to the tracks in the event within this function.
      const larlite::opflash* opflash_pointer = boundary_flash_index_vector.at( opflash_i ).popflash;

      const larlite::opflash  opflash_object  = *opflash_pointer;

      // Print statement for the operation of the 'rankTrackFlashMatches' function.
      //std::cout << "Ranking the track flash matches for a flash." << std::endl;
      //std::cout << "The size of 'larlite_tracks_being_checked' sent into the 'rankTrackFlashMatches' function = " << larlite_tracks_being_checked.size() << "." << std::endl;

      // In this function I will do the naive thing and just find the flash that has the best match with the track, regardless of duplicate matches between the flashes.
      rankTrackFlashMatches( flash_match_config, larlite_tracks_being_checked, well_matched_tracks_idx_v,
			     opflash_object, boundary_flash_index_vector.at( opflash_i ).ivec, opflash_track_match_list, ordered_opflash_track_match_list );

      //std::cout << "The number of tracks that are considered good matches for this track = " << ordered_opflash_track_match_list.size() << "." << std::endl;

      // Append the result to 'all_opflash_candidate_list'.
      all_opflash_candidate_list.emplace_back( std::move(ordered_opflash_track_match_list) );
      // Push back the index of the track corresponding to this flash
      opflash_with_track_lists_idx_v.push_back( boundary_flash_index_vector.at( opflash_i ) );
      
    }

    //std::cout << "Finding a unique flash/track match for each track and flash." << std::endl;
    
    // Feed the 'opflash_track_match_list' to the 'findUniqueTrackFlashMatch' function to find a unique match between each track and flash.
    findUniqueTrackFlashMatch( all_opflash_candidate_list, best_chi2_v, corresponding_qcluster_idx_v );

    // Print statement.
    //std::cout << "The length of 'corresponding_qcluster_idx_v' that I'm getting from the 'findUniqueTrackFlashMatch' function = " << corresponding_qcluster_idx_v.size() << "." << std::endl;
    
    // With this information, you can update the information in 'well_matched_tracks_idx_v'.
    // Loop through the 'corresponding_qcluster_idx_v' list to reset whichever value is currently there
    //  (It will be either -1 if this is in a pass where the A/C tracks were set first, or it will be -1000 if we are looping over an entire pass to look at the tracks).
    for ( size_t qcluster_idx_v_iter = 0; qcluster_idx_v_iter < corresponding_qcluster_idx_v.size(); ++qcluster_idx_v_iter ) {

      // Continue if the chi2 value is negative. That means that the flash did not have any candidates that it could match to.
      if ( best_chi2_v.at( qcluster_idx_v_iter ) < 0.0 ) 
	continue;

      if ( best_chi2_v.at( qcluster_idx_v_iter ) < flash_match_config.chi2_yz_flash_cut ) {
	//std::cout << "qcluster_idx_v_iter of qcluster passing cut = " << qcluster_idx_v_iter << "." << std::endl;
	std::cout << "The chi2 value of the track passing the cuts in the 'flashMatchYZfunction' = " << best_chi2_v.at( qcluster_idx_v_iter ) << "." << std::endl;
	std::cout << "The ivec = " << opflash_with_track_lists_idx_v.at( qcluster_idx_v_iter ).ivec << " and idx = " << opflash_with_track_lists_idx_v.at( qcluster_idx_v_iter ).idx << "." << std::endl;
	well_matched_tracks_idx_v.at( corresponding_qcluster_idx_v.at( qcluster_idx_v_iter ) ) = 1; // corresponding to a good track.
	already_matched_flash_v.push_back( opflash_with_track_lists_idx_v.at( qcluster_idx_v_iter ) ); // append the index of the flash used to determine this track.
      }
      else {
	// Print out the qcluster iter.
	// Print out the flash cut that this qcluster had to satisfy.
	//std::cout << "The chi2 for tracks piercing a yz face = " << flash_match_config.chi2_yz_flash_cut << "." << std::endl;
	//std::cout << "qcluster_idx_v_iter of qcluster failing cut = " << qcluster_idx_v_iter << "." << std::endl;
	well_matched_tracks_idx_v.at( corresponding_qcluster_idx_v.at( qcluster_idx_v_iter ) ) = 0; // corresponding to a bad track.
	impossible_match_endpoint_v.push_back( track_endpoint_v_pass_or_full.at( corresponding_qcluster_idx_v.at( qcluster_idx_v_iter ) ) );
      } 
      
    }
    
    return;
    
  }

  // Declare a function that will create a vector of the struct 'ThruMuTracker::TrackMatchCandidate' for each flash.
  // This will be organized in order of the decreasing chi2 match.
  // Note that not this vector will not include tracks already well-matched to an anode/cathode flash and tracks that are deemed unmatchable to the flash based on the larlite machinery.
  // Inputs: 'flash_match_config' - The configuration object used to perform the flash matching and to initialize the manager for the larlite machinery.
  //         'larlite_tracks_being_checked' - The qclusters, either in the pass or in the entire event, that are being ranked for each of flashes.
  //         'well_matched_tracks_idx_v' - This vector contains information for if a track, in the pass or in the entire event, has already been well-matched to a flash.  If so, we skip it.
  //         'opflash' - This is the opflash product for which we are making the list.
  //         'opflash_producer' - This is the producer of the opflash, meaning either 'simpleFlashCosmic' or 'simpleFlashBeam', for which the flash was collected.
  //         'opflash_track_match_list' - This is the vector of the objects of struct 'FlashMatchCandidate' that is being constructed from the flash in question.
  void ThruMuTracker::rankTrackFlashMatches( GeneralFlashMatchAlgoConfig& flash_match_config, std::vector< larlite::track >& larlite_tracks_being_checked, const std::vector< int >& well_matched_tracks_idx_v, larlite::opflash opflash, int opflash_producer, std::vector< ThruMuTracker::FlashMatchCandidate >& opflash_track_match_list, std::vector< ThruMuTracker::FlashMatchCandidate >& ordered_opflash_track_match_list ) {

    // Declare an object of type 'GeneralFlashMatchAlgo' using the configuration object from 'GeneralFlashMatchAlgoConfig'.                                                                           
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Declare a flash manager for use in this function.                                                                                                                                              
    // Just a test at first to see if these will work......                                                                                                                                           
    ::flashana::FlashMatchManager       _mgr;
    // Configure this with the flashmatch manager available.                                                                                                                                               
    _mgr.Configure( flash_match_config.m_flashmatch_config );
    _mgr.Reset();

    // Convert the opflash into an unfitted data flash.
    flashana::Flash_t unfitted_data_flash = flash_match_obj.MakeDataFlash( opflash );

    // Clear the 'opflash_track_match_list' variable.
    opflash_track_match_list.clear();

    // Print out the size of 'larlite_tracks_being_checked.
    //std::cout << "The size of 'larlite_tracks_being_checked' = " << larlite_tracks_being_checked.size() << "." << std::endl;

    // Loop through the 'qclusters_being_checked' function to generate the 'opflash_track_match_list' that needs to be generated.
    // We will not consider: 1. tracks that already have been well-matched to A/C piercing flashes.
    //                       2. tracks that are impossible matches for this flash.
    for ( size_t qcluster_iter = 0; qcluster_iter < larlite_tracks_being_checked.size(); ++qcluster_iter ) {

      // Print out value of 'well_matched_tracks_idx_v' at this iterator.
      //std::cout << "The value of 'well_matched_track_idx_v' at this iterator = " << well_matched_tracks_idx_v.at( qcluster_iter ) << "." << std::endl;
      
      // Perform the two checks listed above first.
      if ( well_matched_tracks_idx_v.at( qcluster_iter ) == 1 ) {
	//std::cout << "This qcluster is already well-matched to a track.  Continuing!" << std::endl;
	continue;
      }

      // Declare a new larlite track object that is t0-tagged with the flash that determined its A/C endpoints.                                                                                        
      larlite::track t0_tagged_track;
      t0_tagged_track.clear_data();

      // Loop through the entries of the existing larlite track object and fill it with the track points, but t0-tagged.                                                                               
      // This will allow the extensions to take into account the fact that the track should be t0-tagged.                                                                                             
      for ( size_t point_iter = 0; point_iter < larlite_tracks_being_checked.at( qcluster_iter ).NumberTrajectoryPoints(); ++point_iter ) {

        // Declare a TVector object with the trajectory points (t0-tagged) of this entry in 'larlite_track_vector' at this point in its trajectory.                                                         
        TVector3 larlite_trk_pt( larlite_tracks_being_checked.at( qcluster_iter ).LocationAtPoint( point_iter ).X() - opflash.Time()*0.1114, larlite_tracks_being_checked.at( qcluster_iter ).LocationAtPoint( point_iter ).Y(), larlite_tracks_being_checked.at( qcluster_iter ).LocationAtPoint( point_iter ).Z() );

        t0_tagged_track.add_vertex( larlite_trk_pt );
	t0_tagged_track.add_direction( larlite_tracks_being_checked.at( qcluster_iter ).DirectionAtPoint( point_iter ) );

      }
      
      // Print out the number of points on the larlite track.
      //std::cout << "The number of points on the larlite t0-tagged track = " << t0_tagged_track.NumberTrajectoryPoints() << "." << std::endl;

      // Declare a qcluster object that will be extended while taking the t0-tagged x-coordinate into account.                                                                                        
      flashana::QCluster_t t0_tagged_qcluster;
      t0_tagged_qcluster.clear();

      // Use the 'FlashMatchInterface' functionality to expand this qcluster near the boundary.                                                                                                        
      flash_match_obj.ExpandQClusterNearBoundaryFromLarliteTrack( t0_tagged_qcluster, t0_tagged_track, 10000.0, 10.0 );

      // Set the 'time' attribute of the qcluster equal to the time of this opflash.                                                                                                                 
      t0_tagged_qcluster.time = opflash.Time();

      // Print out the number points on the qcluster.
      //std::cout << "The number of points on the qcluster = " << t0_tagged_qcluster.size() << "." << std::endl;

      // Declare a qcluster object that will not be extended when it is compared to the flash to see if they are an impossible match.
      flashana::QCluster_t t0_tagged_qcluster_not_extended;
      t0_tagged_qcluster_not_extended.clear();

      t0_tagged_qcluster_not_extended.time = opflash.Time();

      // Use the 'FlashMatchInterface' functionality to set this qcluster without extending it so that you can use the 'MatchCompatible' function.
      flash_match_obj.ExpandQClusterStartingWithLarliteTrack( t0_tagged_qcluster_not_extended, t0_tagged_track, 0.0, false, false );

      if ( ((flashana::TimeCompatMatch*)(_mgr.GetAlgo(flashana::kMatchProhibit)))->MatchCompatible(t0_tagged_qcluster_not_extended, unfitted_data_flash ) == false ) {
	//std::cout << "This flash and this qcluster are incompatible for the qcluster with index #" << qcluster_iter << " ." << std::endl;
	continue;
      }

      // Declare an object of type 'ThruMuTracker::FlashMatchCandidate' for the flash.
      ThruMuTracker::FlashMatchCandidate candidate;

      // Declare variables for the total PEs of the data flash and the total PEs of the hypothesis flash.
      float totpe_data = 0.0;
      float totpe_hypo = 0.0;

      // Set the chi2 and the index of the qcluster at this particular spot.
      candidate.idx     = qcluster_iter;
      candidate.chi2    = flash_match_obj.generate_chi2_in_track_flash_comparison( t0_tagged_qcluster, opflash, totpe_data, totpe_hypo, opflash_producer );
      
      std::cout << "chi2 = " << candidate.chi2 << "." << std::endl;
     
      // Append this list onto the 'opflash_track_match_list'.
      opflash_track_match_list.push_back( candidate );

    }

    //std::cout << "The size of this particular 'opflash_track_match_list' = " << opflash_track_match_list.size() << "." << std::endl;

    // Use the 'orderInAscendingChi2Order' function to put this list so that its chi2 matches are in descending order.
    // The indices will come along for the ride.
    orderInAscendingChi2Order( opflash_track_match_list, ordered_opflash_track_match_list );

  }

  // Declare a function that will order the chi2 matches of the tracks in descending order, keeping the index of the qcluster matched with the track with the correct chi2.
  // Input: 'opflash_track_match_list' - This is a vector of the 'FlashMatchCandidate' objects that contains both the chi2 and the index of the qclusters that are being matched.
  //        'ordered_opflash_track_match_list' - This is a bector of the 'FlashMatchCandidate' objects that contains the chi2 and the index of the qclusters that are being matched in the correct order.
  // The general idea behind this function is to reorganize all of the members into increasing chi2 order, while being careful to keep the same qcluster index with the flash that it started with.
  void ThruMuTracker::orderInAscendingChi2Order( const std::vector< ThruMuTracker::FlashMatchCandidate >& opflash_track_match_list, std::vector< ThruMuTracker::FlashMatchCandidate >& ordered_opflash_track_match_list ) {

    // This function will take the following form: 
    // 1. Make a vector of all chi2 values in the event.
    // 2. Place them in descending order.
    // 3. Match them with their correct index and place them in 'ordered_opflash_track_match_list'.

    // Declare a variable for the length of the 'opflash_track_match_list'.
    const int array_size = opflash_track_match_list.size();
    
    // Create 'ordered_chi', a vector that contains all of the chi2 values for the elements of 'opflash_track_match_list'.
    float ordered_chi2[ array_size ];  // Change this back to 'ordered_chi2' when you are done with this test.

    // Loop through 'opflash_track_match_list' and fill it with the 'chi2' components of each of the 'FlashMatchCandidate' objects in 'opflash_track_match_list'.
    for ( size_t match_idx = 0; match_idx < array_size; ++match_idx ) {

      float chi2 = opflash_track_match_list.at( match_idx).chi2;

      ordered_chi2[ match_idx ] = chi2;

    }

    // Print out the values in 'ordered_chi2' (change back to 'unordered_chi2' once this is complete).
    //for ( size_t ordered_chi2_iter = 0; ordered_chi2_iter < array_size; ++ordered_chi2_iter ) {

    //  std::cout << "Value of ordered_chi2 at Iter #" << ordered_chi2_iter << " before it is sorted = " << ordered_chi2[ ordered_chi2_iter ] << "." << std::endl; // Print out the value of 'ordered_chi2' here.
      
    //}

    // Use the 'sort()' function of C++ vectors to sort the elements of 'ordered_chi2' from least to greatest (ascending order because the chi2 least in number is the one with the best match to a flash.)
    std::sort( ordered_chi2 , ordered_chi2 + array_size );

    // Print out the values in 'ordered_chi2' after it is sorted.
    //for ( size_t iter = 0; iter < array_size; ++iter ) {

    //  std::cout << "Value of ordered_chi2 at Iter #" << iter << " after it is sorted = " << ordered_chi2[ iter ] << "." << std::endl;

    //}
    
    // Declare an array 

    // Clear 'ordered_opflash_track_match_list' with the intention of filling it with the same information as in 'opflash_track_match_list', but with the chi2 matches ordered from least to greatest.
    ordered_opflash_track_match_list.clear();

    // Loop through 'ordered_chi2' to match the ordered chi2 with its index in 'opflash_track_match_list', filling the information in 'ordered_opflash_track_match_list' as we go.
    // We will have to compare the two to a deal of precision (~10e-7) because some of the chi2 matches will be very close to one another.
    for ( size_t chi2_iter = 0; chi2_iter < array_size; ++chi2_iter ) {

      float ordered_chi2_val = ordered_chi2[ chi2_iter ];

      // Loop through 'opflash_track_match_list' in order to compare the chi2 at each index with the one at this index in 'ordered_chi2'.
      for ( size_t original_order_idx = 0; original_order_idx < opflash_track_match_list.size(); ++original_order_idx ) {

	// Compare the chi2 of the qcluster entry at this index and the value of the ordered chi2 contained above.
	if ( fabs( ordered_chi2_val - opflash_track_match_list.at( original_order_idx).chi2 ) < 0.0000001 ) {

	  //std::cout << "(After a match has been found between the 'ordered_chi2_value' and the 'opflash_track_match_list'." << std::endl;
	  //std::cout << "Ordered chi2 value that is currently being looked at = " << ordered_chi2_val << "." << std::endl;
	  //std::cout << "Opflash_track_match_list chi2 value that is currently being looked at = " << opflash_track_match_list.at( original_order_idx).chi2 << "." << std::endl;
	  //std::cout << "QCluster index of the qcluster with which the flash has this chi2 value = " << opflash_track_match_list.at( original_order_idx ).idx << "." << std::endl;
	  //std::cout << "\n" << std::endl;

	  // Append this value onto 'ordered_opflash_track_match_list' along with its index in a 'FlashMatchCandidate'.
	  ThruMuTracker::FlashMatchCandidate ordered_candidate;
	  ordered_candidate.idx  = opflash_track_match_list.at( original_order_idx).idx;
	  ordered_candidate.chi2 = opflash_track_match_list.at( original_order_idx).chi2;

	  // Make sure that this 'ordered_candidate' has not already been matched once.
	  // The order does not matter if the flash has the same chi2 with two different qclusters. The idea is just to keep two entries in 'ordered_opflash_track_match_list' from having the 
	  // same index value.
	  bool already_matched = false;

	  // Loop through 'ordered_opflash_track_match_list' and make sure an element with the same qcluster idx has not already been added.
	  for ( size_t iter = 0; iter < ordered_opflash_track_match_list.size(); ++iter ) {

	    if ( ordered_candidate.idx == ordered_opflash_track_match_list.at( iter ).idx ) {
	      already_matched = true;
	      break;
	    }
	  }

	  // Continue on with the loop over 'opflash_track_match_list' if 'already_matched' is true.
	  if ( already_matched == true )
	    continue;

	  // Print out the index of the 'ordered_candidate' that you are adding to 'ordered_opflash_track_match_list'.
	  //std::cout << "Index of the 'ordered_candidate' being added to 'ordered_opflash_track_match_list' = " << ordered_candidate.idx << "." << std::endl;

	  // Place 'ordered_candidate' onto the 'ordered_opflash_track_match_list' vector.
	  ordered_opflash_track_match_list.emplace_back( std::move( ordered_candidate ) );

	  // Break the loop at this point, because this chi2 value has been matched.
	  // break; // Take out this break statement for now.
	}

      }

    }

    // Return.  'ordered_opflash_track_match_list' is now filled with pairs of chi2/idx in increasing order.
    return;
      
  }
      
  // Declare a function that will loop through the entries in 'all_opflash_candidate_list' and return the best chi2 for each one and the corresponding qcluster index.
  // Inputs: 'all_opflash_candidate_list': This is a vector of 'FlashMatchCandidate' vectors, with each organized in order of increasing chi2 value.
  //         'best_chi2_v': This returns the best chi2 for each track.  The default value is -1 if the flash has no match, meaning that there is no track that does not have a better chi2 match with another track.  Impossible matches have already been removed.
  //         'corresponding_qcluster_idx_v': This gives the index of the corresponding qcluster.
  void ThruMuTracker::findUniqueTrackFlashMatch( const std::vector< std::vector< ThruMuTracker::FlashMatchCandidate > >& all_opflash_candidate_list, std::vector< float >& best_chi2_v, std::vector< int >& corresponding_qcluster_idx_v ) {

    // Claear the two input vectors for good practice.
    best_chi2_v.clear();
    corresponding_qcluster_idx_v.clear();

    // Declare a vector of type 'ThruMuTracker::UniqueTrackFlashMatch' to hold the matches.
    std::vector< ThruMuTracker::UniqueTrackFlashMatch > unique_track_flash_match_v;
    unique_track_flash_match_v.clear();

    // Make a vector of objects of type 'UniqueTrackFlashMatch'.
    // Initialize this vector with the first element of each vector.
    for ( size_t opflash_candidate_idx = 0; opflash_candidate_idx < all_opflash_candidate_list.size(); ++opflash_candidate_idx ) {

      ThruMuTracker::UniqueTrackFlashMatch first_element_of_list;

      // Continue if the entry of 'all_opflash_candidate_list' does not have any entries.
      if ( all_opflash_candidate_list.at( opflash_candidate_idx).size() == 0 ) { // This should point to the specific set of flashes and its candidates.

	// Fill 'unique_track_flash_match_v' with a dummy object that will maintain the correspondence between the elements in 'unique_track_flash_match_v' and this vector, 'all_opflash_candidate_list'.
	first_element_of_list.chi2              = -1;
	first_element_of_list.idx               = -1;
	first_element_of_list.spot_in_ranking   = -1; // You can set 'first_element_of_list.spot_in_ranking' just equal to the list here.

	// Append this element onto 'unique_track_flash_match_v'.
	unique_track_flash_match_v.emplace_back( std::move( first_element_of_list ) );

	continue;

      }

      // Set it equal to the first element of each vector within 'all_opflash_candidate_list'.
      first_element_of_list.chi2            = all_opflash_candidate_list.at( opflash_candidate_idx).at( 0 ).chi2;
      first_element_of_list.idx             = all_opflash_candidate_list.at( opflash_candidate_idx).at( 0 ).idx;
      first_element_of_list.spot_in_ranking = 0; // This has the best chi2 of the flash matched to any track, so it has the first index.

      // Append 'first_element_of_list' onto 'unique_track_flash_match_v'.
      unique_track_flash_match_v.emplace_back( std::move( first_element_of_list ) );

    }

    // Print the total length of 'all_opflash_candidate_list'.
    //std::cout << "The total length of 'all_opflash_candidate_list' = " << all_opflash_candidate_list.size() << "." << std::endl;

    // Play the 'juggling' game where we find a unique track match for each flash.
    while ( isSameQClusterMatchedToDifferentFlashes( unique_track_flash_match_v ) ) {

      // Print out the qcluster index that each track corresponds to.
      //for ( size_t unique_track_flash_match_iter = 0; unique_track_flash_match_iter < unique_track_flash_match_v.size(); ++unique_track_flash_match_iter ) {

      //std::cout << "Index of QCluster matched to Flash #" << unique_track_flash_match_iter << " = " << unique_track_flash_match_v.at( unique_track_flash_match_iter ).idx << "." << std::endl;
	
      //}

      // Place a print statement here to see if the algorithm is stuck within this loop.
      //std::cout << "Trying to match one qcluster to one flash each." << std::endl;

      // Find the first pair of flashes matched to the same qcluster, displace the flash with the greater (less optimal) chi2 value when matched to 
      int                                    unique_track_flash_match_idx; // index of the flash within 'unique_track_flash_match'.  This index is the same as the one within 'all_opflash_candidate_list'.
      ThruMuTracker::UniqueTrackFlashMatch   duplicate_flash_match_with_worse_chi2;

      // Find these values.
      findQClusterFlashMatchWithWorseChi2( unique_track_flash_match_v, unique_track_flash_match_idx, duplicate_flash_match_with_worse_chi2 );

      // Replace the entry in 'unique_track_flash_match' representing this flash with the next element in its vector in 'all_opflash_candidate_list' ONLY if it is not 
      // the last index in the corresponding vector in 'all_opflash_candidate_list'.
      if ( duplicate_flash_match_with_worse_chi2.spot_in_ranking < int( all_opflash_candidate_list.at( unique_track_flash_match_idx ).size() - 1 ) && duplicate_flash_match_with_worse_chi2.spot_in_ranking > -0.001) { // Make sure we do not look at values that are below our default value.

	// Print out the new flash information for this match.
	//std::cout << "The index in 'unique_track_flash_match_v' and 'all_opflash_candidate_list' corresponding to the flash that got bumped = " << unique_track_flash_match_idx << "." << std::endl;
	//std::cout << "The chi2 for the qcluster in the event with the best match to this flash = " << unique_track_flash_match_v.at( unique_track_flash_match_idx ).chi2 << "." << std::endl;
	//std::cout << "The qcluster index that corresponds to the best match for this flash = " << unique_track_flash_match_v.at( unique_track_flash_match_idx ).idx << "." << std::endl;
	//std::cout << "The spot in the ranking that this qcluster corresponds to for this flash = " << unique_track_flash_match_v.at( unique_track_flash_match_idx ).spot_in_ranking << "." << std::endl;

	// Change the entry in 'unique_track_flash_match' corresponding to this flash to the next element in its rank.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).chi2 = all_opflash_candidate_list.at( unique_track_flash_match_idx).at( duplicate_flash_match_with_worse_chi2.spot_in_ranking + 1 ).chi2; // The next chi2.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).idx  = all_opflash_candidate_list.at( unique_track_flash_match_idx).at( duplicate_flash_match_with_worse_chi2.spot_in_ranking + 1).idx; // The next index.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).spot_in_ranking = duplicate_flash_match_with_worse_chi2.spot_in_ranking + 1;

      }

      // Otherwise, there is no qcluster left for the flash to be matched to.  There could have been matches removed by the 'MatchCompatible' function, as well as every other qcluster having a better match with another flash.
      else {

	//std::cout << "Entering the 'else' loop within the 'findUniqueTrackMatch' function." << std::endl;

	// Set the 'unique_track_flash_match' entry equal to the default values, and set the 'spot_in_ranking' variable equal to the size of the all_opflash_candidate_list.at( unique_track_flash_match_idx) 'TrackFlashMatch' vector.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).chi2            = -1.0;
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).idx             = -1;
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).spot_in_ranking = -1;

	// Print out the index that this corresponds to and the information that has been reset.
	//std::cout << "Index in the 'else' loop with information that is being reset: " << unique_track_flash_match_idx << "." << std::endl;
	//std::cout << "chi2 of the 'unique_track_flash_match_v' in the vector = " << unique_track_flash_match_v.at( unique_track_flash_match_idx ).chi2 << "." << std::endl;
	//std::cout << "idx of the 'unique_track_flash_match_v' in the vector = " << unique_track_flash_match_v.at( unique_track_flash_match_idx).idx << "." << std::endl;
	//std::cout << "spot_in_ranking of the 'unique_track_flash_match_v' in the vector = " << unique_track_flash_match_v.at( unique_track_flash_match_idx ).spot_in_ranking << "." << std::endl;
	

      }

    }

    // After each of the flashes are matched to a unique qcluster or set to a default value (the 'else' loop above sets the default values), output the best chi2 value and the corresponding qcluster index value from the original list of tracks being looped over.
    for ( size_t unique_track_flash_match_iter = 0; unique_track_flash_match_iter < unique_track_flash_match_v.size(); ++unique_track_flash_match_iter ) {

      best_chi2_v.push_back( unique_track_flash_match_v.at( unique_track_flash_match_iter ).chi2 );
      corresponding_qcluster_idx_v.push_back( unique_track_flash_match_v.at( unique_track_flash_match_iter ).idx ); // This corresponds to the index of the qcluster in the original list.

    }

    return;

  }
      
  // Declare a function that will find out if two flashes in an object of type 'std::vector< ThruMuTracker::UniqueTrackFlashMatch >' are matched to qclusters with the same idx.
  // Inputs : 'unique_track_flash_match_v' - contains 'UniqueTrackFlashMatch' objects, which contain the index of the matched qcluster, the chi2 of the matched qcluster, and the spot in the flash's ranking where this qcluster is contained.
  bool ThruMuTracker::isSameQClusterMatchedToDifferentFlashes( const std::vector< ThruMuTracker::UniqueTrackFlashMatch >& unique_track_flash_match_v ) {

    // Loop through the 'unique_track_flash_match_v' in two separate loops to find if any two sets of qcluster indices correspond to one another.
    for ( size_t outer_idx = 0; outer_idx < unique_track_flash_match_v.size(); ++outer_idx ) {

      // Declare a variable for the outer index of the best qcluster matched to this flash.
      int outer_qcluster_idx = unique_track_flash_match_v.at( outer_idx).idx;

      for ( size_t inner_idx = outer_idx + 1; inner_idx < unique_track_flash_match_v.size(); ++inner_idx ) {
	
	int inner_qcluster_idx = unique_track_flash_match_v.at( inner_idx).idx;

	// Return 'false' if 'outer_qcluster_idx' equals 'inner_qcluster_idx'.
	if ( outer_qcluster_idx == inner_qcluster_idx && outer_qcluster_idx > -0.001 )  {// Eliminate the possibility that both flashes do not have a match. 
	  // Print the two indices here (they are probably both 0).
	  //std::cout << "outer_qcluster_idx = " << outer_qcluster_idx << "at flash index #" << outer_idx << "." << std::endl;
	  //std::cout << "inner_qcluster_idx = " << inner_qcluster_idx << "at flash index #" << inner_idx << "." << std::endl;
	  return true;
	}

      }

    }

    // If the qcluster indices for each of the 'ThruMuTracker::UniqueTrackFlashMatch' objects are distinct, then return 'false'.
    return false;

  }

  // Declare a function that will find one pair of flashes matched to the same qcluster, compare their chi2 values, and return the information for the worse flash/qcluster match so that its information can be reassigned.
  // This function will keep being called until none of the flashes are matched to the same qcluster.
  // Inputs: unique_track_flash_match_v            - This is the vector of 'ThruMuTracker::UniqueTrackFlashMatch' objects that contains the information to be studied.
  //         unique_track_flash_match_idx          - This is the index of the 'ThruMuTracker::UniqueTrackFlashMatch' information of the flash with the worse chi2 match to the same qcluster.
  //         duplicate_flash_match_with_worse_chi2 - This is the 'ThruMuTracker::UniqueTrackFlashMatch' object corresponding to the flash with the worse chi2 match.
  void ThruMuTracker::findQClusterFlashMatchWithWorseChi2( const std::vector< ThruMuTracker::UniqueTrackFlashMatch >& unique_track_flash_match_v, int& unique_track_flash_match_idx,
					    ThruMuTracker::UniqueTrackFlashMatch& duplicate_flash_match_with_worse_chi2 ) {

    // Loop through the 'unique_track_flash_match_v' vector to find the first two entries that are matched to the same qcluster, and return the information for the one with the worse (greater) chi2 value.
    for ( size_t outer_idx = 0; outer_idx < unique_track_flash_match_v.size(); ++outer_idx ) {

      int   outer_qcluster_idx = unique_track_flash_match_v.at( outer_idx ).idx;
      float outer_chi2         = unique_track_flash_match_v.at( outer_idx ).chi2;

      for (size_t inner_idx = outer_idx+1; inner_idx < unique_track_flash_match_v.size(); ++inner_idx ) {

	int   inner_qcluster_idx = unique_track_flash_match_v.at( inner_idx ).idx;
	float inner_chi2         = unique_track_flash_match_v.at( inner_idx ).chi2;

	// In this loop, make sure that the chi2 of on
	if ( outer_qcluster_idx == inner_qcluster_idx && outer_qcluster_idx > -0.001 ) { // Add the second part so there are not two flashes that are matched to the same qcluster.

	  //std::cout << "In the 'findQClusterFlashMatchWithWorseChi2', the two flash indices that are matched to one another are = " << inner_idx << " and = " << outer_idx << "." << std::endl;
	  
	  // Set 'unique_track_flash_match_idx' and 'duplicate_flash_match_with_worse_chi2' based on which track has a greater chi2 value, and then return.
	  if ( outer_chi2 < inner_chi2 ) {

	    unique_track_flash_match_idx           = inner_idx; // Make sure to set this to the LOOP INDEX, not to the index of the qcluster.
	    duplicate_flash_match_with_worse_chi2  = unique_track_flash_match_v.at(inner_idx );

	  }

	  // The outer chi2 is the worse chi2 value.
	  else {

	    unique_track_flash_match_idx           = outer_idx;
	    duplicate_flash_match_with_worse_chi2  = unique_track_flash_match_v.at( outer_idx ); // Make sure to set this to the LOOP INDEX, not to the index of the qcluster.
	  }

	  // This point necessarily must be reached, because if two of the qcluster indices were not indentical then the function would not have been called.
	  return;

	}

      }

    }

  }

  // Declare a function that will check to see which of the two A/C tracks matched to the same track makes a better match.
  // Input: 'larlite_track': This is the larlite track currently under investigation that has not yet been matched to the track.
  //        'flash_info': This is the flash info for the flash that was under threshold for both tracks.
  //        'trackclusters': This is the vector of 'BMTrackCluster3D' objects that will hold the information about the last track matched to this flash.
  //        'track_endpoint_v': This is the vector of track endpoints that will reveal which was the last track matched to this flash.
  //        'track_endpoint_flash_v': This vector contains information for which flash (if any) determined the endpoint for a track.  This is now a vector of vectors.                                   
  //        'track_endpoint_boundary_type_idx_v': This vector contains information for the flash producer used to generate the flash that produced the endpoint.  This is now a vector of vectors.  
  //        Note: I will fill 'trackclusters', 'track_endpoint_v', 'track_endpoint_flash_v', and 'track_boundary_type_idx_v' with empty objects if the new flash has a better match than the old flash.
  bool ThruMuTracker::compareToOtherACTrackMatchedToThisFlash( GeneralFlashMatchAlgoConfig& flash_match_config, const larlite::track& larlite_track, const larlitecv::BoundaryFlashIndex& flash_info, std::vector< larlitecv::BMTrackCluster3D >& trackclusters, std::vector< std::vector< larlitecv::BoundarySpacePoint > >& track_endpoint_v, std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v, std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v ) {

    //std::cout << "At the top of the 'compareToOtherACTrackMatchedToThisFlash'." << std::endl;

    // Create a new object of class 'GeneralFlashMatchAlgo'.                                                                                                                                          
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Get the 'ivec' and the 'idx' of the flash matched to both of these flashes.
    int flash_ivec = flash_info.ivec;
    int flash_idx  = flash_info.idx;

    // Declare a variable for the number of times that this flash is matched to endpoints in the 'track_endpoint_v' vector.  If it is less than 2, than the previous track that this flash was matched to was a YZ flash (non-A/C-piercing), so it is too difficult to compare the two of them.  If there are two or more, than we can compare the two tracks head-to-head and keep the one with the better match to the flash.
    int num_of_tracks_matched_to_this_flash = 0;
    std::vector< int > idx_of_tracks_matched_to_this_flash;
    idx_of_tracks_matched_to_this_flash.clear();

    for ( size_t i = 0; i < track_endpoint_v.size(); ++i ) {

      // Continue if the size of one of the elements is 0.  That means that this function was already used once to eliminate those endpoints from the the 'track_endpoint_v' vector.
      if ( track_endpoint_v.at( i ).size() == 0 ) 
	continue;

      // Declare variables for the track endpoint information.
      larlitecv::BoundaryFlashIndex first_endpoint_info  = track_endpoint_v.at( i )[0].getFlashIndex();
      larlitecv::BoundaryFlashIndex second_endpoint_info = track_endpoint_v.at( i )[1].getFlashIndex();

      // Compare the two endpoints (separately) to the flash information contained in the information above.
      if ( ( flash_ivec == first_endpoint_info.ivec && flash_idx == first_endpoint_info.idx ) || ( flash_ivec == second_endpoint_info.ivec && flash_idx == second_endpoint_info.idx ) ) {

	++num_of_tracks_matched_to_this_flash;
	idx_of_tracks_matched_to_this_flash.push_back( i );

      }

    }

    // Return false if 'num_of_tracks_matched_to_this_flash' is less than 2.
    if ( num_of_tracks_matched_to_this_flash < 2 ) {
      std::cout << "This track was previously matched to a non-A/C-piercing flash; it is too difficult to recover that track and make the comparison." << std::endl;
      std::cout << "'num_of_tracks_matched_to_this_flash' = " << num_of_tracks_matched_to_this_flash << ", continuing!!" << std::endl;
      return false;
    }

    // Otherwise, compare the element in 'trackclusters' to the input larlite track to this function to see which one has the better match to the flash.
    else {

      // Print the length of 'num_of_tracks_matched_to_this_flash' (should be 2, otherwise something is wrong).
      std::cout << "'num_of_tracks_matched_to_this_flash' = " << num_of_tracks_matched_to_this_flash << "." << std::endl;

      // Generate an opflash object out of the input 'BoundaryFlashIndex' object.
      const larlite::opflash* opflash_pointer = flash_info.popflash;
      const larlite::opflash  opflash         = *opflash_pointer;
      
      // Generate qclusters out of the input 'larlite_track' (the second track matched to this flash) and the track corresponding to the first index in 'trackclusters' (which was previously matched to this flash).

      // Declare 'larlite_track_input_to_this_function_t0_tagged', and fill it with the t0-tagged trajectory points of the points on the track input to this function.
      larlite::track larlite_track_input_to_this_function_t0_tagged;
      larlite_track_input_to_this_function_t0_tagged.clear_data();
      
      for ( size_t traj_iter = 0; traj_iter < larlite_track.NumberTrajectoryPoints(); ++traj_iter ) {

	TVector3 traj_point( larlite_track.LocationAtPoint( traj_iter )[0] - opflash.Time()*0.1114, larlite_track.LocationAtPoint( traj_iter )[1], larlite_track.LocationAtPoint( traj_iter )[2] );

	larlite_track_input_to_this_function_t0_tagged.add_vertex( traj_point );
	larlite_track_input_to_this_function_t0_tagged.add_direction( larlite_track.DirectionAtPoint( traj_iter ) );  // Add the direction as well, which stays the same with a translation in any Cartesian direction.

      }

      // Declare a qcluster object that will be extended while taking the t0-tagged x-coordinate into account.                                                                                        
      flashana::QCluster_t t0_tagged_qcluster_input_to_this_function;
      t0_tagged_qcluster_input_to_this_function.clear();

      // Use the 'FlashMatchInterface' functionality to expand this qcluster near the boundary.                                                                                                       
      flash_match_obj.ExpandQClusterNearBoundaryFromLarliteTrack(t0_tagged_qcluster_input_to_this_function, larlite_track_input_to_this_function_t0_tagged, 10000.0, 10.0 );

      // Set the 'time' attribute of the qcluster equal to the time of this opflash.                                                                                                                       
      t0_tagged_qcluster_input_to_this_function.time = opflash.Time();

      // Declare 'first_track_matched_to_this_flash_t0_tagged', and fill it with the t0-tagged trajectory points of the first track matched to this flash.
      larlite::track first_track_matched_to_this_flash_t0_tagged;
      first_track_matched_to_this_flash_t0_tagged.clear_data();

      // Declare an object for the track in the trackclusters vector that corresponds to the first track matched to this flash.
      larlite::track first_track_matched_to_this_flash = trackclusters.at( idx_of_tracks_matched_to_this_flash.at( 0 ) ).makeTrack();

      // Loop through 'first_track_matched_to_this_flash' to fill 'first_track_matched_to_this_flash_t0_tagged' with the t0-tagged trajectory points of this flash.
      for ( size_t traj_iter = 0; traj_iter < first_track_matched_to_this_flash.NumberTrajectoryPoints(); ++traj_iter ) {

	TVector3 traj_point( first_track_matched_to_this_flash.LocationAtPoint( traj_iter )[0] - opflash.Time()*0.1114, first_track_matched_to_this_flash.LocationAtPoint(traj_iter )[1], first_track_matched_to_this_flash.LocationAtPoint(traj_iter )[2] );

	first_track_matched_to_this_flash_t0_tagged.add_vertex( traj_point );
	first_track_matched_to_this_flash_t0_tagged.add_direction( first_track_matched_to_this_flash.DirectionAtPoint( traj_iter ) ); // Add the direction as well, which stays the same with a translation in any Cartesian direction.

      }

      // Declare a qcluster object that will be extended while taking the t0-tagged x-coordinate into account.
      flashana::QCluster_t qcluster_first_track_matched_to_this_flash_t0_tagged;
      qcluster_first_track_matched_to_this_flash_t0_tagged.clear();

      flash_match_obj.ExpandQClusterNearBoundaryFromLarliteTrack( qcluster_first_track_matched_to_this_flash_t0_tagged, first_track_matched_to_this_flash_t0_tagged, 10000.0, 10.0 );

      // Set the 'time' attribute of the qcluster equal to the time of this opflash.
      qcluster_first_track_matched_to_this_flash_t0_tagged.time = opflash.Time();

      // Now, compare the chi2 matches for both of these qclusters to the flash.

      // Declare variables for the total PEs given by each of these flashes.
      float data_PE_input_track                  = 0.0;
      float hypo_PE_input_track                  = 0.0;

      float data_PE_first_track_matched_to_flash = 0.0;
      float hypo_PE_first_track_matched_to_flash = 0.0;
      
      float chi2_of_track_input_to_this_function = flash_match_obj.generate_chi2_in_track_flash_comparison( t0_tagged_qcluster_input_to_this_function, opflash, data_PE_input_track, hypo_PE_input_track, flash_ivec );
      float chi2_of_first_track_matched_to_flash = flash_match_obj.generate_chi2_in_track_flash_comparison( qcluster_first_track_matched_to_this_flash_t0_tagged, opflash, data_PE_first_track_matched_to_flash, hypo_PE_first_track_matched_to_flash, flash_ivec );

      // Print these values out.
      std::cout << "The chi2 of the track input to the 'compareToOtherACTrackMatchedToThisFlash' = " << chi2_of_track_input_to_this_function << "." << std::endl;
      std::cout << "The chi2 of the track previously matched to this flash in the 'compareToOtherACTrackMatchedToThisFlash' = " << chi2_of_first_track_matched_to_flash << "." << std::endl;

      // Compare the two chi2 values.  If the current track has a lower chi2, then empty out element of 'trackclusters' at the index located at 'idx_of_tracks_matched_to_this_flash.at( 0 )'.  This allows for two A/C tracks matched to the same flash to be compared between passes in the same way.  The element corresponding to the track previously matched to this flash ('trackclusters'), the element corresponding to the endpoints of the track previously matched to this flash ('track_endpoint_v'), the element corresponding to the flash information for this track ('track_endpoint_flash_v'), and the element corresponding to the boundary endpoint information for this track ('track_endpoint_boundary_type_idx_v') will all be set to empty values to keep the symmetry of those four lists with respect to each pass.
      if ( chi2_of_track_input_to_this_function < chi2_of_first_track_matched_to_flash ) {

	std::cout << "The new input track has a better chi2 match.  Removing the information in the four vectors corresponding to the previous track." << std::endl;

	// Go through each of those lists and declare an object of the list, empty it, and replace the element in the list at index 'idx_of_tracks_matched_to_this_flash.at( 0 )' to this empty element.
	
	// 'BMTrackCluster3D' vector.
	larlitecv::BMTrackCluster3D empty_trackcluster;
	empty_trackcluster.path3d.clear();

	trackclusters.at( idx_of_tracks_matched_to_this_flash.at( 0 ) ) = empty_trackcluster;

	// 'Clear' the element at this entry of the track endpoint vector, 'track_endpoint_v'.
	track_endpoint_v.at( idx_of_tracks_matched_to_this_flash.at( 0 ) ).clear();

	// 'Clear' the element at this entry of the flash information corresponding to the track endpoints, 'track_endpoint_flash_v'.
	track_endpoint_flash_v.at( idx_of_tracks_matched_to_this_flash.at( 0 ) ).clear();

	// 'Clear' the element at this entry of the flash information corresponding to the types of the boundary endpoints, 'track_endpoint_boundary_type_idx_v'. 
	track_endpoint_boundary_type_idx_v.at( idx_of_tracks_matched_to_this_flash.at( 0 ) ).clear();

	// Return 'true', because the new flash is better than the old one.
	return true;

      }

      else {

	std::cout << "The old track had the better chi2 match.  Set 'well_matched_tracks_v' at the new tracks' index to 0." << std::endl;

	// Return 'false', because the originally matched track is better matched to the flash than this one.
	return false;

      }

    }

  }	
      
       

  // Declare a function that will get rid of bad tracks found by ThruMu and tagged in each of the events.
  // Input: trackclusters: the list of track clusters from the event to be parsed through and limited.
  //        well_matched_tracks_idx: This is a list of the tracks that are well-matched to a flash and therefore should be marked in an image.
  //        tracks_under_consideration: This determines how many of the tracks from the back of 'trackclusters' you want to be considered for removal from the vector based on a chi2 match.
  //        single_pass: This boolean determines whether you want to remove bad tracks from a single pass or from the entire vector.  If 'true', then we will remove from the last pass of the tracks ( from 'tracks_per_pass') the number of tracks that are poorly reconstructed. 

  void ThruMuTracker::sortOutBadTracksAndEndpoints( std::vector < larlitecv::BMTrackCluster3D >& trackclusters, std::vector< std::vector< larlitecv::BoundarySpacePoint > >& track_endpoint_v, std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v, std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v, const std::vector< int >& well_matched_tracks_idx_v, std::vector< int >& tracks_per_pass, int tracks_under_consideration, bool single_pass ) {

    // Print out the lengths of 'trackclusters' and 'track_endpoint_v' to make sure that they are the same length.
    //std::cout << "The length of 'trackclusters' at the start of the 'sortOutBadTracksAndEndpoints' function = " << trackclusters.size() << "." << std::endl;
    //std::cout << "The length of 'track_endpoint_v' at the start of the 'sortOutBadTracksAndEndpoints' function = " << track_endpoint_v.size() << "." << std::endl;
    //std::cout << "The length of 'track_endpoint_flash_v' at the start of 'sortOutBadTracksAndEndpoints' function = " << track_endpoint_flash_v.size() << "." << std::endl;
    //std::cout << "The length of 'track_endpoint_boundary_type_idx_v' at the start of the 'sortOutBadTracksAndEndpoints' function = " << track_endpoint_boundary_type_idx_v.size() << "." << std::endl;

    int num_of_tracks_kept           = 0;
    int num_of_tracks_kept_but_empty = 0;

    // Start a loop through 'trackclusters', looping from 'trackclusters.size() - 1' through 'trackclusters.size() - 1 - ( well_matched_tracks_idx.size() - 1)'.
    for ( int  well_matched_tracks_index = 0; well_matched_tracks_index < tracks_under_consideration; ++well_matched_tracks_index ) {

      // Loop through the tracks to see if any of these tracks have a '0' in their index for the 'well_matched_tracks' vector.  That means they should be removed from the 'trackclusters'
      // vector.
      // tracks_per_pass.at( tracks_per_pass.size() - 1 ) == well_matched_tracks_idx.size().
      // This means that the track was already found to be a poorly-matched track.  If it was a good track, then it was labelled '1'.  If the test was inconclusive, then it was labelled with '-1'.
	if ( well_matched_tracks_idx_v.at( well_matched_tracks_index ) == 0 ) {

	  // Remove that item from 'trackclusters'.
	  // Convert the 'int' that iterates in this loop to the correct type.

	  // Make the element of 'trackclusters' at this index into a track and compare its information to that in the function above.
	  larlite::track track = trackclusters.at( trackclusters.size() - tracks_under_consideration + well_matched_tracks_index ).makeTrack();

	  // Print out the information of the track being removed.
	  // std::cout << "Removing the track with y-coordinate = " << track.LocationAtPoint( 0 )[1] << " on the first coordinate of its trajectory at index #" << well_matched_tracks_index << "." << std::endl;

	  // Take these tracks off the END of 'trackclusters', NOT from the beginning; some of these tracks have already been included in the previous loop.
	  trackclusters.erase( trackclusters.end() - tracks_under_consideration + well_matched_tracks_index ); // This corresponds to the number element that 'well_matched_tracks_index' is located at.
	  track_endpoint_v.erase( track_endpoint_v.end() - tracks_under_consideration + well_matched_tracks_index ); // This should be fine if the lengths of 'trackclusters' and 'track_endpoint_v' remain consistent.
	  track_endpoint_flash_v.erase( track_endpoint_flash_v.end() - tracks_under_consideration + well_matched_tracks_index );
	  track_endpoint_boundary_type_idx_v.erase( track_endpoint_boundary_type_idx_v.end() - tracks_under_consideration + well_matched_tracks_index );

	  if ( single_pass ) {
	    // Subtract '1' from the number of tracks reconstructed per event.
	    tracks_per_pass.at( tracks_per_pass.size() - 1 ) -= 1;
	  }

	}

	else {

	  ++num_of_tracks_kept;
	  
	  if (trackclusters.at( trackclusters.size() - tracks_under_consideration + well_matched_tracks_index ).path3d.size() == 0 ) {
	    ++num_of_tracks_kept_but_empty;
	  }

	}
	  
    
    }

    // Print out the variables for the number of tracks kept and the number of tracks kept that are empty.
    std::cout << "The number of tracks passing the flash-matching cuts = " << num_of_tracks_kept << "." << std::endl;
    std::cout << "The number of tracks passing the flash-matching cuts that are empty = " << num_of_tracks_kept_but_empty << "." << std::endl;
    
    return;

  }

}





	
	
    
