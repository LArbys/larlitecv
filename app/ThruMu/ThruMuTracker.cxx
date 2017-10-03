#include "ThruMuTracker.h"

#include "RadialEndpointFilter.h"
#include "AStar3DAlgo.h"
#include "AStarNodes2BMTrackCluster3D.h"
#include "ThruMuFoxExtender.h"
#include "PushBoundarySpacePoint.h"

// larcv 
#include "UBWireTool/UBWireTool.h"


namespace larlitecv {

  ThruMuTracker::ThruMuTracker( const ThruMuTrackerConfig& config )
    : m_config(config)
  {}

  // Add additional arguments: 'flash_producer_idx_v' - a map that  will contain '0' if the flash was produced using 'simpleFlashBeam' and '1' if the flash was produced using 'simpleFlashCosmic'.
  //                           'flash_match_config'         - an object that sets the 'GeneralFlashMatchAlgo' object within the 'flashMatchTracks' function.
  void ThruMuTracker::makeTrackClusters3D( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector<larcv::Image2D>& img_v,  
					   const std::vector<larcv::Image2D>& badchimg_v, const std::vector< const BoundarySpacePoint* >& spacepts,
					   std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
					   std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices, const std::vector< larlite::event_opflash* >& opflashsets, std::vector< int >                                           & flash_idx_v, std::vector< int >& boundary_type_idx_v ) {
    // This method takes in the list of boundaryspacepoints and pairs them up in order to try to find through-going muons
    // input:
    //   img_v: vector images, one for each plane. assumed to be in U,V,Y order
    //   spacepts: all possible boundary crossing points
    // output:
    //   trackclusters: thrumu tracks. collection of space points and pixel2dclusters
    //   tagged_v: pixels tagged as thrumu.  the value indicates the pass when the muon was tagged.
    //   used_endpoints_indices: indicates which endpoints in spacepts was used. 0=not used. 1=used.

    // Declare a vector, 'impossible_match_endpoints', meant to save the indices of endpoints that together constituted an impossible match.  This removes from consideration endpoints that were found to fail the flashmatching stage of reconstruction.                                                                                                                                                
    std::vector< std::vector< int > > impossible_match_endpoint_idx_v;
    impossible_match_endpoint_idx_v.clear();

    // Declare a new vector, 'track_endpoint_indices', based off the indices in 'spacepts', for the endpoints that are used to make tracks.
    std::vector < std::vector < int > > track_endpoint_indices;
    track_endpoint_indices.clear();

    // Declare a vector for the indices of track endpoints that are oriented according to which flash that they correspond to.
    std::vector< std::vector< int > > track_endpoint_flash_idx_v;
    track_endpoint_flash_idx_v.clear();

    // Declare a vector for the indices of the flash producer that generated the flashes: either 'simpleFlashBeam' or 'simpleFlashCosmic'.
    std::vector< std::vector< int >  > track_endpoint_boundary_type_idx_v;
    track_endpoint_boundary_type_idx_v.clear();

    std::vector < int > already_matched_flash_idx_v;
    already_matched_flash_idx_v.clear();

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

    std::vector<int> tracks_per_pass;
    for (int ipass=0; ipass<m_config.num_passes; ipass++) {
      
      // Clear 'well_matched_tracks_idx_v' at the start of the pass.
      well_matched_tracks_idx_v.clear();
      
      const ThruMuTrackerConfig::ThruMuPassConfig& passcfg = m_config.pass_configs.at(ipass);
      runPass( ipass, passcfg, spacepts, img_v, badchimg_v, tagged_v, used_endpoints_indices, trackclusters, flash_idx_v, boundary_type_idx_v, track_endpoint_flash_idx_v, track_endpoint_boundary_type_idx_v, track_endpoint_indices );

      // Print out the information for the 'track_endpoint_flash_idx_v' here in this function.
      std::cout << "Number of filtered endpoints in 'spacepts' = " << spacepts.size() << "." << std::endl;
      std::cout << "Size of 'flash_idx_v' in 'makeTrackClusters3D' = " << flash_idx_v.size() << "." << std::endl;
      std::cout << "Size of 'boundary_type_idx_v' in 'makeTrackClusters3D' = " << boundary_type_idx_v.size() << "." << std::endl;

      bool anode_and_cathode_only = true;
      
      int tracks_in_pass = trackclusters.size();
      for ( int i = int( tracks_per_pass.size() - 1 ); i > -1; --i ) {
        tracks_in_pass -= tracks_per_pass.at( i );
      }
      tracks_per_pass.push_back( tracks_in_pass );

      // Depending if this is greater than the first pass, then set 'anode_and_cathode_only' to 'false'.
      if ( ipass > 0 ) 
	anode_and_cathode_only = false;

      // Use the flag from the config file to determine if you want to flashmatch the tracks that survive this pass of the thrumu tracker.
      if (m_config.thrumu_flashmatch == true ) {

	bool single_pass = true;

        // Call a function that uses all of the flashmatching infrastructure developed in 'GeneralFlashMatchAlgo'.
	flashMatchTracks( flash_match_config, img_v, tagged_v, spacepts, opflashsets, trackclusters, impossible_match_endpoint_idx_v, already_matched_flash_idx_v, well_matched_tracks_idx_v, tracks_in_pass, track_endpoint_flash_idx_v, track_endpoint_boundary_type_idx_v, track_endpoint_indices, anode_and_cathode_only );	
	sortOutBadTracks( trackclusters, well_matched_tracks_idx_v, tracks_per_pass, tracks_per_pass.at( tracks_per_pass.size() - 1 ), single_pass ); 

      }
      
    }

    // tag the pixels
    for (int itrack=0; itrack<(int)trackclusters.size(); itrack++ ) {
      BMTrackCluster3D& track3d = trackclusters[itrack];
      track3d.markImageWithTrack( img_v, badchimg_v, m_config.pixel_threshold, m_config.tag_neighborhood, tagged_v, 0.3, itrack+1 );
    }

    if ( m_config.verbosity>0 ) {
      for (int ipass=0; ipass<m_config.num_passes; ipass++) {
        std::cout << "Number of tracks found in ThruMu pass #" << ipass << ": " << tracks_per_pass.at(ipass) << std::endl;
      }
    }

  }

  // Update this function to pass the index of the flash that is matched to the endpoints of the track along with the indices themselves.
  // New arguments: flash_idx_v: This vector of ints contains the indices of the flash (organized in absolute order using the flash producer) corresponding to the endpoint that is being looped over.
  //                boundary_type_idx_v: This vector of ints contains the producer of the flash that corresponds to the endpoints at this point in the vector.
  //                track_endpoint_flash_idx_v: This vector of ints corresponds to the index of the flash that is matched to the track endpoints at this point in the vector.
  //                track_endpoint_boundary_type_idx_v: This vector of ints contains the producer of the flash that corresponds to the endpoints of the track in the 'trackclusters' vector.
  void ThruMuTracker::runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
			       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
			       std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters, const std::vector< int >& flash_idx_v,
			       const std::vector< int >& boundary_type_idx_v, std::vector< std::vector< int > >& track_endpoint_flash_idx_v, 
			       std::vector< std::vector< int > >& track_endpoint_boundary_type_idx_v, std::vector< std::vector< int > >& track_endpoint_indices ) {

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

    // Declare analogous vectors for the indices of the track as they correspond to the flash.
    std::vector< std::vector<int> > endpoint_flash_idx_v;
    std::vector< std::vector<int> > endpoint_boundary_type_idx_v;

    // Declare the vector for the matched track indices up here.
    std::vector< int > track_flash_idx_v;
    std::vector< int > track_boundary_type_idx_v;

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
	  track_flash_idx_v.resize(2);
	  track_flash_idx_v[0]          = flash_idx_v.at( i );
	  track_flash_idx_v[1]          = flash_idx_v.at( j );

	  track_boundary_type_idx_v.resize(2);
	  track_boundary_type_idx_v[0] = boundary_type_idx_v.at( i );
	  track_boundary_type_idx_v[1] = boundary_type_idx_v.at( j );

	  // Print out the flash information for the two endpoints.                                                                                                                                   
	  std::cout << "boundary_type_idx_v.at( i ) = " << boundary_type_idx_v.at( i ) << "." << std::endl;
	  std::cout << "boundary_type_idx_v.at( j ) = " << boundary_type_idx_v.at( j ) << "." << std::endl;
	  
          if ( m_config.verbosity>1 ) {
            std::cout << "#### Storing track. size=" << track3d.path3d.size() << ". indices (" << indices[0] << "," << indices[1] << ") ####" << std::endl;
          }
          pass_end_indices.emplace_back( std::move(indices) );
          pass_track_candidates.emplace_back( std::move(track3d) );
	  
	  endpoint_flash_idx_v.emplace_back( std::move(track_flash_idx_v) );
	  endpoint_boundary_type_idx_v.emplace_back( std::move(track_boundary_type_idx_v) );
        }

      }// second loop over end points
    }// first loop over end points

    // Print out the length of the two sets of endpoint producers.
    std::cout << "Pairs of flash endpoint indices = " << endpoint_flash_idx_v.size() << "." << std::endl;
    std::cout << "Pairs of flash endpoint producer indices = " << endpoint_boundary_type_idx_v.size() << "." << std::endl;

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
      
      // Push back 'indices' into 'track_endpoint_indices'.
      track_endpoint_indices.push_back( indices );
      
    }

    // Use the same logic to place the indices of the flash corresponding to the endpoints and the flash producer corresponding to the flash that determined the endpoints.
    // Flash index values.
    // Use a loop so this will not segfault.
    for ( auto& track_flash_idx_v: endpoint_flash_idx_v ) {
      std::cout << "Appending track endpoints to the 'track_endpoint_flash_idx_v' vector." << std::endl;
      track_endpoint_flash_idx_v.emplace_back( std::move( track_flash_idx_v ) );
      // Print out the length of this vector.
      std::cout << "track_endpoint_flash_idx_v vector size = " << track_endpoint_flash_idx_v.size() << "." << std::endl;
    }
    // Flash producer index values.
    for( auto& track_boundary_type_idx_v: endpoint_boundary_type_idx_v) {
      std::cout << "Appending track endpoints to the 'track_endpoint_boundary_type_idx_v' vector." << std::endl;
      track_endpoint_boundary_type_idx_v.emplace_back( std::move( track_boundary_type_idx_v ) );
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

    larlitecv::AStar3DAlgo algo( pass_cfg.astar3d_cfg );
    std::vector<larlitecv::AStar3DNode> path;
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
  // already_matched_flash_idx_v: This a vector of the flashes that are already determined (in a previous pass) to have a good match to another track and have therefore been removed from consideration for other tracks.
  // num_of_tracks_added_in_pass: The number of tracks added in the pass for which we are flashmatching information.
  // track_endpoint_flash_idx_v: This vector contains information for which flash (if any) determined the endpoint for a track.  This is now a vector of vectors.
  // track_endpoint_boundary_type_idx_v: This vector contains information for the flash producer used to generate the flash that produced the endpoint.  This is now a vector of vectors. 
  // anode_and_cathode_only: This only compares the endpoints for anode-piercing/cathode-piercing tracks, for which one flash is matched to the track.
  void ThruMuTracker::flashMatchTracks( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& tagged_v, const std::vector< const BoundarySpacePoint* >& spacepts, const std::vector< larlite::event_opflash* >& opflash_v, std::vector< BMTrackCluster3D >& trackclusters, std::vector< std::vector< int >  >& impossible_match_endpoint_idx_v, std::vector< int >& already_matched_flash_idx_v, std::vector< int >& well_matched_tracks_idx_v, const int& num_of_tracks_added_in_pass, std::vector< std::vector< int > >& track_endpoint_flash_idx_v, std::vector< std::vector< int > >& track_endpoint_boundary_type_idx_v, std::vector< std::vector< int > >& track_endpoint_indices, bool anode_and_cathode_only ) {

    // Declare an object of type 'GeneralFlashMatchAlgo' using the configuration object from 'GeneralFlashMatchAlgoConfig'.
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Declare a flash manager for use in this function.
    // Just a test at first to see if these will work......
    ::flashana::FlashMatchManager       _mgr;
    // Configure this with the flashmatch manager available.
    _mgr.Configure( flash_match_config.m_flashmatch_config );
    _mgr.Reset();

    // Declare three vectors: one for the BMTrackCluster3D objects, one for the flash indices of the track endpoints (two per track), and one for the flash producer of the track endpoints (two per track). Clear each vector to ensure that there are no issues with space in the vector.
    std::vector< larlitecv::BMTrackCluster3D > trackclusters_from_pass;
    trackclusters_from_pass.clear();
    std::vector< std::vector< int > > track_endpoint_indices_from_pass;
    track_endpoint_indices_from_pass.clear();
    std::vector< std::vector< int > > track_endpoint_flash_idx_v_from_pass;
    track_endpoint_flash_idx_v_from_pass.clear();
    std::vector< std::vector< int > > track_endpoint_boundary_type_idx_v_from_pass;
    track_endpoint_boundary_type_idx_v_from_pass.clear();

    // Take a given track added in the last pass and compare it to all of the flashes in the event, finding one flash closer to it than any of the others.
    for ( int track_last_pass_iter = int( trackclusters.size() - num_of_tracks_added_in_pass); track_last_pass_iter < int( trackclusters.size() ); ++track_last_pass_iter ) {

      trackclusters_from_pass.push_back( trackclusters.at( track_last_pass_iter ) );
      track_endpoint_indices_from_pass.push_back( track_endpoint_indices.at( track_last_pass_iter ) );
      
      // The trajectory points corresponding to this track are located at indices '2i' and '2i + 1' within the 'track_endpoint_flash_idx_v' and 'track_endpoint_boundary_type_idx_v'.               
      // Do not invert the order of the two points as they are within the original lists over the flash indices.                                                                                      
      track_endpoint_flash_idx_v_from_pass.push_back( track_endpoint_flash_idx_v.at( track_last_pass_iter ) );                                                                                    
      track_endpoint_boundary_type_idx_v_from_pass.push_back( track_endpoint_boundary_type_idx_v.at( track_last_pass_iter ) );                                                                 
    }

    // Turn this into a vector of larlite tracks using the the 'GeneralFlashMatchAlgo' functionality.
    std::vector < larlite::track > larlite_track_vector = flash_match_obj.generate_tracks_between_passes( trackclusters_from_pass );

    // Turn this into a vector of qclusters and expand it.
    std::vector < flashana::QCluster_t > qcluster_vector;
    qcluster_vector.clear();

    // For now, assume that all flashes are reconstructed in the cosmic discriminator window.
    // Loop through the tracks.
    for ( size_t i = 0; i < larlite_track_vector.size(); ++i ) {
      
      // Make each track into an expanded qcluster.
      flashana::QCluster_t qcluster;
      flash_match_obj.ExpandQClusterStartingWithLarliteTrack( qcluster, larlite_track_vector.at( i ), 10000., true, true ); 

      qcluster_vector.push_back( qcluster );

    }

    // Generate a single opflash vector using the functionality in GeneralFlashMatchAlgo.
    std::vector< larlite::opflash > single_opflash_vector = flash_match_obj.generate_single_opflash_vector_for_event( opflash_v );

    // Generate a vector of the producer corresponding to these flashes as well.
    std::vector< int > single_opflash_producer_idx_v = flash_match_obj.generate_single_opflash_idx_vector_for_event( opflash_v );

    // Check the boolean 'anode_and_cathode_only'.  If 'true', then we will only consider the anode-piercing/cathode-piercing tracks, which are determined by a single flash.
    if ( anode_and_cathode_only ) {
      flashMatchAC( flash_match_config, qcluster_vector, single_opflash_vector, single_opflash_producer_idx_v, track_endpoint_flash_idx_v_from_pass, track_endpoint_boundary_type_idx_v_from_pass, impossible_match_endpoint_idx_v, track_endpoint_indices_from_pass, well_matched_tracks_idx_v, already_matched_flash_idx_v );
      return;
    }

    else {

      bool entire_event = false;

      flashMatchAC( flash_match_config, qcluster_vector, single_opflash_vector, single_opflash_producer_idx_v, track_endpoint_flash_idx_v_from_pass, track_endpoint_boundary_type_idx_v_from_pass, impossible_match_endpoint_idx_v, track_endpoint_indices_from_pass, well_matched_tracks_idx_v, already_matched_flash_idx_v );
      flashMatchYZFaceTracks( flash_match_config, trackclusters, qcluster_vector, single_opflash_vector, single_opflash_producer_idx_v, impossible_match_endpoint_idx_v, well_matched_tracks_idx_v,  already_matched_flash_idx_v, entire_event );
	return;

    }

  }

  // Declare a function for matching the anode-piercing/cathode-piercing tracks to their corresponding flash.
  // Inputs:
  // qcluster_vector: The vector of qclusters, fully extended outside the TPC and all.
  // single_opflash_vector: This is the vector of 'opflash' products, which are the flashes in the event organized first with those reconstructed using 'simpleFlashBeam' and second with those reconstructed using 'simpleFlashCosmic'.
  // 'track_endpoint_flash_idx_v: This is a vector of the indices of the flash that determined the endpoint of the track.  This is created with the same scheme that the 'single_opflash_vector' was created.
  // 'track_endpoint_boundary_type_idx_v': This is a vector of the indices of the flash producer that determined the endpoint of the track.  This is created with the same scheme that the 'single_opflash_vector' was created.
  // 'impossible_match_endpoints': This is a vector of the endpoints that cannot form a valid track based on the information shown here.
  // 'already_matched_flash_idx_v': This is a vector of the indices of the flashes that have already been well-matched to an anode-piercing/cathode-piercing track.
  void ThruMuTracker::flashMatchAC( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< flashana::QCluster_t >& qcluster_vector, 
				    const std::vector< larlite::opflash >& single_opflash_vector, const std::vector< int >& single_opflash_producer_idx_v, 
				    const std::vector< std::vector< int > >& track_endpoint_flash_idx_v_from_pass, const std::vector< std::vector< int > >& track_endpoint_boundary_type_idx_v_from_pass, 
				    std::vector< std::vector< int > >& impossible_match_endpoint_idx_v, std::vector< std::vector < int > > & track_endpoint_indices_from_pass, 
				    std::vector < int >& well_matched_tracks_idx_v, std::vector< int >& already_matched_flash_idx_v ) {

    // Create a new object of class 'GeneralFlashMatchAlgo'.
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );
 
    // Loop through the qclusters added in the pass.  If both of the entries in 'track_endpoint_flash_idx_v_from_pass' are set to -1, then that means that the track was not made from an anode-piercing or cathode-piercing endpoint.
    for ( size_t qcluster_iter = 0; qcluster_iter < qcluster_vector.size(); ++qcluster_iter ) { 

      // This means that the track was not determined by an anode-piercing/cathode-piercing flash and so we have to continue.
      if ( track_endpoint_flash_idx_v_from_pass.at( qcluster_iter )[0] < 0 && track_endpoint_flash_idx_v_from_pass.at( qcluster_iter )[1] < 0 ) {
	well_matched_tracks_idx_v.push_back( -1 );
	continue;
      }

      // If the track passes that selection, then it is anode-piercing or cathode-piercing.  Declare an object for the opflash that determined the track endpoints based on which index in
      // 'track_endpoint_flash_idx_v' is greater than 0.
      // I could have also gotten the track producer information directly from the list of opflashes........an improvement later.
      int endpoint_flash_idx;
      int endpoint_boundary_type_idx;

      if ( track_endpoint_flash_idx_v_from_pass.at( qcluster_iter )[0] > -0.001 ) {
	endpoint_flash_idx                   = track_endpoint_flash_idx_v_from_pass.at( qcluster_iter )[0];
	endpoint_boundary_type_idx           = track_endpoint_boundary_type_idx_v_from_pass.at( qcluster_iter )[0];
      }

      if ( track_endpoint_flash_idx_v_from_pass.at( qcluster_iter )[1] > -0.001 ) {
	endpoint_flash_idx                   = track_endpoint_flash_idx_v_from_pass.at( qcluster_iter )[1];
	endpoint_boundary_type_idx           = track_endpoint_boundary_type_idx_v_from_pass.at( qcluster_iter )[1];
      } 

      bool impossible_match_endpoints = false;

      // Make sure that the endpoints corresponding to this track do not belong to the endpoints already found to not constitute a real track.
      for ( size_t impossible_match = 0; impossible_match < impossible_match_endpoint_idx_v.size(); ++impossible_match ) {

	// Check either possibility of the endpoints being identical to those found impossible to connect tracks in this stage of the tagger.
	if ( ( track_endpoint_indices_from_pass.at( qcluster_iter )[0] == impossible_match_endpoint_idx_v.at( impossible_match)[0] && track_endpoint_indices_from_pass.at( qcluster_iter )[1] == impossible_match_endpoint_idx_v.at( impossible_match)[1] ) || ( track_endpoint_indices_from_pass.at( qcluster_iter )[0] == impossible_match_endpoint_idx_v.at( impossible_match)[1] && track_endpoint_indices_from_pass.at( qcluster_iter )[1] == impossible_match_endpoint_idx_v.at( impossible_match)[0] ) ) {
	  impossible_match_endpoints = true;
	}
      }

      // Append a '0' onto 'well_matched_tracks_idx_v' because these are a set of points already found to contain a bad track, so the track must be bad.
      if ( impossible_match_endpoints == true ) {
	well_matched_tracks_idx_v.push_back( 0 );
	continue;
      }

      bool determined_by_well_matched_flash = false;

      // If 'enpoint_flash_idx' is equal to any of the entries of 'already_matched_flash_idx', then continue.
      for ( size_t well_matched_flash_iter = 0; well_matched_flash_iter < already_matched_flash_idx_v.size(); ++well_matched_flash_iter ) {
	if ( already_matched_flash_idx_v.at( well_matched_flash_iter ) == endpoint_flash_idx ) {
	  determined_by_well_matched_flash = true;
	}
      }
      
      // This flash was already well-matched to another flash, so here we have to append a '0' onto the vector.
      if ( determined_by_well_matched_flash == true ) {
	well_matched_tracks_idx_v.push_back( 0 );
	continue;
      }

      // This will have to be taken out when we do performance studies.....
      // Only look at the cathode-piercing tracks right now.
      //if ( endpoint_boundary_type_idx == larlitecv::kCathode ) continue;

      // Set the opflash object.
      larlite::opflash opflash_object    = single_opflash_vector.at( endpoint_flash_idx );

      // Set the producer that made the opflash object.
      int opflash_producer_idx           = single_opflash_producer_idx_v.at( endpoint_flash_idx );

      // Create a copy of the qcluster that is t0-tagged with the endpoint of its track.
      flashana::QCluster_t qcluster = qcluster_vector.at( qcluster_iter );
      
      for ( size_t qpt = 0; qpt < qcluster.size(); ++qpt ) {
	qcluster.at( qpt ).x -= opflash_object.Time()*0.1114;
      }

      // Determine the chi2 match between the qcluster and the opflash object.
      // Now we are getting the opflash producer from the entire list.
      float chi2 = flash_match_obj.generate_chi2_in_track_flash_comparison( qcluster, opflash_object, opflash_producer_idx );

      std::cout << "The chi2 for the match of track #" << qcluster_iter << " among the tracks added in this pass to the qcluster vector = " << chi2 << "." << std::endl;

      // Put a chi2 cut on the match between the flash and the qcluster.
      if ( chi2 < flash_match_config.chi2_anode_cathode_cut ) { // to be defined somewhere... 
	
	// Add the flash to the list of 'matched' flashes 
	already_matched_flash_idx_v.push_back( endpoint_flash_idx );
	
	// Append the index of the qcluster to the list of well-matched tracks.
	well_matched_tracks_idx_v.push_back( 1 );

      }

      else { // If the track's chi2 was greater than the anode/cathode piercing chi2 cut value, then put the track's endpoints in the 'impossible_match_endpoint_idx_v' list.
	
	impossible_match_endpoint_idx_v.push_back( track_endpoint_indices_from_pass.at( qcluster_iter ) );

	well_matched_tracks_idx_v.push_back( 0 );

      }
	

    }
	
    return;

  }

  // Declare a function that will loop through the non-anode-piercing/cathode-piercing tracks to match them to the correct flash.
  // Input: 'flash_match_config' - This is the 'GeneralFlashMatchAlgoConfig' object that will be used to initialize the 'GeneralFlashMatchAlgo' object.
  //        'trackclusters'      - This is the total vector of remaining track clusters if we would like to look at the tracks in the entire event instead of just those 
  //        'qcluster_vector' - This is the vector of 'qclusters' from the pass that will be compared to the flashes to return the best match.
  //        'single_opflash_vector' - This is vector of opflashes in the event contained in a single list, those produced by 'simpleFlashBeam' and those produced by 'simpleFlashCosmic'.  
  //        'single_opflash_producer_idx_v' - This is a vector of the producer of the flashes - '0' for 'simpleFlashBeam' and '1' for 'simpleFlashCosmic'.
  //        'impossible_match_endpoint_idx_v' - This is a vector of vectors that contains the endpoints that are an impossible match with one another.
  //        'well_matched_tracks_idx' - This vector contains the information from the tracks in the event that are well-matched already to a single flash.
  //        'already_matched_flash_idx_v' - This vector contains the indices of the flashes that are already well-matched to a track in the event, meaning that they can be passed over.
  //        'entire_event' - This boolean indicator asks if we want to loop through all of the information contained in the event or in a pass after the 'ACflashMatch'.
  void ThruMuTracker::flashMatchYZFaceTracks( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< larlitecv::BMTrackCluster3D >& trackclusters, const std::vector < flashana::QCluster_t > qcluster_vector, const std::vector< larlite::opflash >& single_opflash_vector, const std::vector< int >& single_opflash_producer_idx_v, std::vector< std::vector< int > >& impossible_match_endpoint_idx_v, std::vector< int >& well_matched_tracks_idx_v, std::vector< int >& already_matched_flash_idx_v, bool entire_event ) {

    // Declare a 'GeneralFlashMatchAlgo' object with the 'flash_match_config' object passed to 'flashMatchAllTracks'.
    larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_config );

    // Use the flag 'entire_event' to determine what the scope of flashes that we are looking at is.  This will allow for the definition of a vector of qclusters, called 'qclusters_being_checked'.
    // If '!entire_event', then we will only look at those qclusters in the input 'qcluster_vector', which only contains qclusters from this pass.
    // If 'entire_event', then we will look at those qclusters from the entire event, finding some flash that they are compatible with in the event.
    std::vector< flashana::QCluster_t > qclusters_being_checked;
    qclusters_being_checked.clear();

    // Declare vectors of the best chi2 for each flash and the corresponding qcluster index for the best qcluster that it is matched to.
    std::vector< float > best_chi2_v;
    best_chi2_v.clear();
    std::vector< int > corresponding_qcluster_idx_v;
    corresponding_qcluster_idx_v.clear();

    if ( entire_event) {
      
      // Convert 'trackclusters' (the entire vector) to a vector of larlite tracks.
      std::vector< larlite::track > larlite_track_vector = flash_match_obj.generate_tracks_between_passes( trackclusters );

      // Loop through the 'larlite_track_vector', convert each track to an extended qcluster, and append to 'qclusters_being_checked'.
      for ( size_t trk = 0; trk < larlite_track_vector.size(); ++trk ) {

	flashana::QCluster_t qcluster;
	flash_match_obj.ExpandQClusterStartingWithLarliteTrack( qcluster, larlite_track_vector.at( trk ), 10000., true, true );
	qclusters_being_checked.emplace_back( std::move(qcluster) );

      }

    }

    else {

	qclusters_being_checked = qcluster_vector;
	well_matched_tracks_idx_v.clear();
    }

    // Declare an empty vector of vector of 'ThruMuTracker::FlashMatchCandidate' candidates.
    std::vector< std::vector< ThruMuTracker::FlashMatchCandidate > > all_opflash_candidate_list;
    all_opflash_candidate_list.clear();

    // Declare a vector for the flashes that are not well-matched already, for which a vector of tracks has to be compared.
    std::vector< int > opflash_with_track_lists_idx_v;
    opflash_with_track_lists_idx_v.clear();

    // Loop through the eligible flashes, using the 'struct' declared above to create a 'list' ranking the tracks that are matched to each flash in an event.
    for ( size_t opflash_i = 0; opflash_i < single_opflash_vector.size(); ++opflash_i ) {

      bool already_matched_idx = false;

      // Make sure that this flash does not have the same index as one of the indices in the 'already_matched_flash_idx_v' machine.
      for ( size_t already_matched_iter = 0; already_matched_iter < already_matched_flash_idx_v.size(); ++already_matched_iter ) {

	if ( int(opflash_i) == already_matched_flash_idx_v.at( already_matched_iter ) ) 
	  already_matched_idx = true;
      }

    // Continue if 'already_matched_idx' is true.
    if ( already_matched_idx == true ) continue;

    std::vector< ThruMuTracker::FlashMatchCandidate > opflash_track_match_list;
    opflash_track_match_list.clear();

    // In this function I will do the naive thing and just find the flash that has the best match with the track, regardless of duplicate matches between the flashes.
    rankTrackFlashMatches( flash_match_config, qclusters_being_checked, well_matched_tracks_idx_v, single_opflash_vector.at(opflash_i) , single_opflash_producer_idx_v.at( opflash_i ), opflash_track_match_list );

    // Append the result to 'all_opflash_candidate_list'.
    all_opflash_candidate_list.emplace_back( std::move(opflash_track_match_list) );
    // Push back the index of the track corresponding to this flash
    opflash_with_track_lists_idx_v.push_back( opflash_i );
 
    }
 
    // Feed the 'opflash_track_match_list' to the 'findUniqueTrackFlashMatch' function to find a unique match between each track and flash.
    findUniqueTrackFlashMatch( all_opflash_candidate_list, best_chi2_v, corresponding_qcluster_idx_v );

  // With this information, you can update the information in 'well_matched_tracks_idx_v'.
  // Loop through the 'corresponding_qcluster_idx_v' list to reset whichever value is currently there (It will be either -1 if this is in a pass where the A/C tracks were set first, or it will be -1000 if we are looping over an entire pass to look at the tracks).
  for ( size_t qcluster_idx_v_iter = 0; qcluster_idx_v_iter < corresponding_qcluster_idx_v.size(); ++qcluster_idx_v_iter ) {

    if ( best_chi2_v.at( qcluster_idx_v_iter ) < flash_match_config.chi2_yz_flash_cut ) {
      well_matched_tracks_idx_v.at( qcluster_idx_v_iter ) = 1; // corresponding to a good track.
      already_matched_flash_idx_v.push_back( opflash_with_track_lists_idx_v.at( qcluster_idx_v_iter ) ); // append the index of the flash used to determine this track.
    }
    else {
      well_matched_tracks_idx_v.at( qcluster_idx_v_iter ) = 0; // corresponding to a bad track.
    } // We won't worry about the endpoints for now, but maybe later that can be added in.
      
  }

  return;

  }

  // Declare a function that will create a vector of the struct 'ThruMuTracker::TrackMatchCandidate' for each flash.
  // This will be organized in order of the decreasing chi2 match.
  // Note that not this vector will not include tracks already well-matched to an anode/cathode flash and tracks that are deemed unmatchable to the flash based on the larlite machinery.
  // Inputs: 'flash_match_config' - The configuration object used to perform the flash matching and to initialize the manager for the larlite machinery.
  //         'qclusters_being_checked' - The qclusters, either in the pass or in the entire event, that are being ranked for each of flashes.
  //         'well_matched_tracks_idx_v' - This vector contains information for if a track, in the pass or in the entire event, has already been well-matched to a flash.  If so, we skip it.
  //         'opflash' - This is the opflash product for which we are making the list.
  //         'opflash_producer' - This is the producer of the opflash, meaning either 'simpleFlashCosmic' or 'simpleFlashBeam', for which the flash was collected.
  //         'opflash_track_match_list' - This is the vector of the objects of struct 'FlashMatchCandidate' that is being constructed from the flash in question.
  void ThruMuTracker::rankTrackFlashMatches( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< flashana::QCluster_t >& qclusters_being_checked, const std::vector< int >& well_matched_tracks_idx_v, larlite::opflash opflash, int opflash_producer, std::vector< ThruMuTracker::FlashMatchCandidate >& opflash_track_match_list ) {

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

    // Loop through the 'qclusters_being_checked' function to generate the 'opflash_track_match_list' that needs to be generated.
    // We will not consider: 1. tracks that already have been well-matched to A/C piercing flashes.
    //                       2. tracks that are impossible matches for this flash.
    for ( size_t qcluster_iter = 0; qcluster_iter < qclusters_being_checked.size(); ++qcluster_iter ) {

      // Perform the two checks listed above first.
      if ( well_matched_tracks_idx_v.at( qcluster_iter ) == 1 ) continue;

      flashana::QCluster_t qcluster_copy = qclusters_being_checked.at( qcluster_iter );

      // make a copy of the qcluster, and t0-tag that copy with information from the flash.
      for ( size_t qpt = 0; qpt < qcluster_copy.size(); ++qpt ) {
	qcluster_copy.at( qpt ).x -= opflash.Time()*0.1114;
      }


      if ( ((flashana::TimeCompatMatch*)(_mgr.GetAlgo(flashana::kMatchProhibit)))->MatchCompatible(qcluster_copy, unfitted_data_flash ) == false )
	continue;

      // Declare an object of type 'ThruMuTracker::FlashMatchCandidate' for the flash.
      ThruMuTracker::FlashMatchCandidate candidate;

      // Set the chi2 and the index of the qcluster at this particular spot.
      candidate.idx     = qcluster_iter;
      candidate.chi2    = flash_match_obj.generate_chi2_in_track_flash_comparison( qcluster_copy, opflash, opflash_producer );
     
      // Append this list onto the 'opflash_track_match_list'.
      opflash_track_match_list.push_back( candidate );

    }

    // Declare an object for the ordered chi2/flash match list.
    std::vector< ThruMuTracker::FlashMatchCandidate > ordered_opflash_track_match_list;
    ordered_opflash_track_match_list.clear();

    // Use the 'orderInDescendingChi2Order' function to put this list so that its chi2 matches are in descending order.
    // The indices will come along for the ride.
    orderInDescendingChi2Order( opflash_track_match_list, ordered_opflash_track_match_list );

  }

  // Declare a function that will order the chi2 matches of the tracks in descending order, keeping the index of the qcluster matched with the track with the correct chi2.
  // Input: 'opflash_track_match_list' - This is a vector of the 'FlashMatchCandidate' objects that contains both the chi2 and the index of the qclusters that are being matched.
  //        'ordered_opflash_track_match_list' - This is a bector of the 'FlashMatchCandidate' objects that contains the chi2 and the index of the qclusters that are being matched in the correct order.
  // The general idea behind this function is to reorganize all of the members into increasing chi2 order, while being careful to keep the same qcluster index with the flash that it started with.
  void ThruMuTracker::orderInDescendingChi2Order( const std::vector< ThruMuTracker::FlashMatchCandidate >& opflash_track_match_list, std::vector< ThruMuTracker::FlashMatchCandidate >& ordered_opflash_track_match_list ) {

    // This function will take the following form: 
    // 1. Make a vector of all chi2 values in the event.
    // 2. Place them in descending order.
    // 3. Match them with their correct index and place them in 'ordered_opflash_track_match_list'.
    
    // Create 'ordered_chi', a vector that contains all of the chi2 values for the elements of 'opflash_track_match_list'.
    std::vector< float > unordered_chi2;
    unordered_chi2.clear();

    // Loop through 'opflash_track_match_list' and fill it with the 'chi2' components of each of the 'FlashMatchCandidate' objects in 'opflash_track_match_list'.
    for ( size_t match_idx = 0; match_idx < opflash_track_match_list.size(); ++match_idx ) {

      float chi2 = opflash_track_match_list.at( match_idx).chi2;

      unordered_chi2.emplace_back( std::move( chi2 ) );

    }

    // Use the 'sort()' function of C++ vectors to sort the elements of 'ordered_chi2' from least to greatest (ascending order because the chi2 least in number is the one with the best match to a flash.)
    // Sort the 'ordered_chi2' function in descending order with this homemade function - C++ gave me too much trouble in setting it up.
    // Declare a vector for values already used.
    std::vector< float > values_already_used;
    values_already_used.clear();

    // Declare the vector of ordered values.
    std::vector< float > ordered_chi2;
    ordered_chi2.clear();

    // Declare a float to be used within the loop.
    float min_remaining;

    // Run the outer loop 'unordered_chi2.size()' number of times to correctly place the items in ascending order.
    for ( size_t outer = 0; outer < unordered_chi2.size(); ++outer ) {

      // Loop through each element of 'ordered_chi2' and place them in ascending order.
      for ( size_t i = 0; i < unordered_chi2.size(); ++i ) {

	// This is larger than any of the chi2 values.
	min_remaining = 1e20;

	bool already_ordered = false;

	// See if this value of 'unordered_chi2' has been set yet by looping over 'values_already_used'.
	for ( size_t j = 0; j < values_already_used.size(); ++j ) {

	  if ( fabs(unordered_chi2.at( i ) - values_already_used.at( j ) ) < 0.0000001 ) {
	    already_ordered = true;
	    break;
	  }

	}

	// Continue if the chi2 value has already been ordered.
	if ( already_ordered == true ) 
	  continue;

	// If it hasn't, compare it to 'min_remaining'.
	if ( unordered_chi2.at( i ) < min_remaining ) {

	  // Reset 'min_remaining'.
	  min_remaining = unordered_chi2.at( i );

	}

      }

      // After you find the value of 'min_remaining', append it to both 'ordered_chi2' and 'values_already_used'.
      ordered_chi2.push_back( min_remaining );
      values_already_used.push_back( min_remaining );

    }
      

    // Declare an array 

    // Clear 'ordered_opflash_track_match_list' with the intention of filling it with the same information as in 'opflash_track_match_list', but with the chi2 matches ordered from least to greatest.
    ordered_opflash_track_match_list.clear();

    // Loop through 'ordered_chi2' to match the ordered chi2 with its index in 'opflash_track_match_list', filling the information in 'ordered_opflash_track_match_liat' as we go.
    // We will have to compare the two to a deal of precision (~10e-7) because some of the chi2 matches will be very close to one another.
    for ( size_t chi2_iter = 0; chi2_iter < ordered_chi2.size(); ++chi2_iter ) {

      float ordered_chi2_val = ordered_chi2.at( chi2_iter );

      // Loop through 'opflash_track_match_list' in order to compare the chi2 at each index with the one at this index in 'ordered_chi2'.
      for ( size_t original_order_idx = 0; original_order_idx < opflash_track_match_list.size(); ++original_order_idx ) {

	// Compare the chi2 of the qcluster entry at this index and the value of the ordered chi2 contained above.
	if ( fabs( ordered_chi2_val - opflash_track_match_list.at( original_order_idx).chi2 ) < 0.0000001 ) {

	  // Append this value onto 'ordered_opflash_track_match_list' along with its index in a 'FlashMatchCandidate'.
	  ThruMuTracker::FlashMatchCandidate ordered_candidate;
	  ordered_candidate.idx  = opflash_track_match_list.at( original_order_idx).idx;
	  ordered_candidate.chi2 = opflash_track_match_list.at( original_order_idx).chi2;

	  // Place 'ordered_candidate' onto the 'ordered_opflash_track_match_list' vector.
	  ordered_opflash_track_match_list.emplace_back( std::move( ordered_candidate ) );

	  // Break the loop at this point, because this chi2 value has been matched.
	  break;
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

      // Initialize a 'UniqueTrackFlashMatch' object to be used for this object.
      ThruMuTracker::UniqueTrackFlashMatch first_element_of_list;

      // Set it equal to the first element of each vector within 'all_opflash_candidate_list'.
      first_element_of_list.chi2            = all_opflash_candidate_list.at( opflash_candidate_idx).at( 0 ).chi2;
      first_element_of_list.idx             = all_opflash_candidate_list.at( opflash_candidate_idx).at( 0 ).idx;
      first_element_of_list.spot_in_ranking = 0; // This has the best chi2 of the flash matched to any track, so it has the first index.

      // Append 'first_element_of_list' onto 'unique_track_flash_match_v'.
      unique_track_flash_match_v.emplace_back( std::move( first_element_of_list ) );

    }

    // Play the 'juggling' game where we find a unique track match for each flash.
    while ( isSameQClusterMatchedToDifferentFlashes( unique_track_flash_match_v ) ) {

      // Find the displaced flash, and continue onto the next spot in its ranked list for one of them.
      int                                unique_track_flash_match_idx; // index of the flash within 'unique_track_flash_match'.  This index is the same as the one within 'all_opflash_candidate_list'.
      ThruMuTracker::UniqueTrackFlashMatch   duplicate_flash_match_with_worse_chi2;

      // Find these values.
      findQClusterFlashMatchWithWorseChi2( unique_track_flash_match_v, unique_track_flash_match_idx, duplicate_flash_match_with_worse_chi2 );

      // Replace the entry in 'unique_track_flash_match' representing this flash with the next element in its vector in 'all_opflash_candidate_list' ONLY if it is not 
      // the last index in the corresponding vector in 'all_opflash_candidate_list'.
      if ( duplicate_flash_match_with_worse_chi2.spot_in_ranking < int(all_opflash_candidate_list.at( unique_track_flash_match_idx ).size() - 1) ) {

	// Change the entry in 'unique_track_flash_match' corresponding to this flash to the next element in its rank.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).chi2 = all_opflash_candidate_list.at( unique_track_flash_match_idx).at( duplicate_flash_match_with_worse_chi2.spot_in_ranking + 1 ).chi2; // The next chi2.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).idx  = all_opflash_candidate_list.at( unique_track_flash_match_idx).at( duplicate_flash_match_with_worse_chi2.spot_in_ranking + 1).idx; // The next index.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).spot_in_ranking = duplicate_flash_match_with_worse_chi2.spot_in_ranking + 1;

      }

      // Otherwise, there is no qcluster left for the flash to be matched to.
      else {

	// Set the 'unique_track_flash_match' entry equal to the default values, and set the 'spot_in_ranking' variable equal to the size of the all_opflash_candidate_list.at( unique_track_flash_match_idx) 'TrackFlashMatch' vector.
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).chi2            = -1.0;
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).idx             = -1.0;
	unique_track_flash_match_v.at( unique_track_flash_match_idx ).spot_in_ranking = all_opflash_candidate_list.at( unique_track_flash_match_idx).size(); // One greater than the maximum index.
      }

    }

    // After each of the flashes are matched to a unique qcluster or set to a default value (the 'else' loop above sets the default values), output the best chi2 value and the corresponding qcluster index value from the original list of tracks being looped over.
    for ( size_t unique_track_flash_match_iter = 0; unique_track_flash_match_iter < unique_track_flash_match_v.size(); ++unique_track_flash_match_iter ) {

      best_chi2_v.push_back( unique_track_flash_match_v.at( unique_track_flash_match_iter ).chi2 );
      corresponding_qcluster_idx_v.push_back( unique_track_flash_match_v.at( unique_track_flash_match_iter ).idx );

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
	if ( outer_qcluster_idx == inner_qcluster_idx ) 
	  return true;

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
  void ThruMuTracker::findQClusterFlashMatchWithWorseChi2( const std::vector< ThruMuTracker::UniqueTrackFlashMatch >& unique_track_flash_match_v, int& unique_track_flash_match_idx, ThruMuTracker::UniqueTrackFlashMatch& duplicate_flash_match_with_worse_chi2 ) {

    // Loop through the 'unique_track_flash_match_v' vector to find the first two entries that are matched to the same qcluster, and return the information for the one with the worse (greater) chi2 value.
    for ( size_t outer_idx = 0; outer_idx < unique_track_flash_match_v.size(); ++outer_idx ) {

      int   outer_qcluster_idx = unique_track_flash_match_v.at( outer_idx ).idx;
      float outer_chi2         = unique_track_flash_match_v.at( outer_idx ).chi2;

      for (size_t inner_idx = 0; inner_idx < unique_track_flash_match_v.size(); ++inner_idx ) {

	int   inner_qcluster_idx = unique_track_flash_match_v.at( inner_idx ).idx;
	float inner_chi2         = unique_track_flash_match_v.at( inner_idx ).chi2;

	if ( outer_qcluster_idx == inner_qcluster_idx ) {
	  
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
      


  // Declare a function that will get rid of bad tracks found by ThruMu and tagged in each of the events.
  // Input: trackclusters: the list of track clusters from the event to be parsed through and limited.
  //        well_matched_tracks_idx: This is a list of the tracks that are well-matched to a flash and therefore should be marked in an image.
  //        tracks_under_consideration: This determines how many of the tracks from the back of 'trackclusters' you want to be considered for removal from the vector based on a chi2 match.
  //        single_pass: This boolean determines whether you want to remove bad tracks from a single pass or from the entire vector.  If 'true', then we will remove from the last pass of the tracks ( from 'tracks_per_pass') the number of tracks that are poorly reconstructed. 

  void ThruMuTracker::sortOutBadTracks( std::vector < larlitecv::BMTrackCluster3D >& trackclusters, const std::vector< int >& well_matched_tracks_idx_v, std::vector< int >& tracks_per_pass, int tracks_under_consideration, bool single_pass ) {

    // Start a loop through 'trackclusters', looping from 'trackclusters.size() - 1' through 'trackclusters.size() - 1 - ( well_matched_tracks_idx.size() - 1)'.
    for ( size_t well_matched_tracks_index = 0; well_matched_tracks_index < tracks_under_consideration; ++well_matched_tracks_index ) {

      // Loop through the tracks to see if any of these tracks have a '0' in their index for the 'well_matched_tracks' vector.  That means they should be removed from the 'trackclusters'
      // vector.
      // tracks_per_pass.at( tracks_per_pass.size() - 1 ) == well_matched_tracks_idx.size().
      // This means that the track was already found to be a poorly-matched track.  If it was a good track, then it was labelled '1'.  If the test was inconclusive, then it was labelled with '-1'.
	if ( well_matched_tracks_idx_v.at( well_matched_tracks_index ) == 0 ) {
	  
	  // Remove that item from 'trackclusters'.
	  // Convert the 'int' that iterates in this loop to the correct type.
	  trackclusters.erase( trackclusters.begin() + (well_matched_tracks_index + 1) ); // This corresponds to the number element that 'well_matched_tracks_index' is located at.

	  if ( single_pass ) {
	    // Subtract '1' from the number of tracks reconstructed per event.
	    tracks_per_pass.at( tracks_per_pass.size() - 1 ) -= 1;
	  }

	}
    
      }
    
    return;

  }

}





	
	
    
