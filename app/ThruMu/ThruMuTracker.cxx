#include "ThruMuTracker.h"

namespace larlitecv {

  ThruMuTracker::ThruMuTracker( const ThruMuTrackerConfig& config )
    : m_config(config)
  {}

  void ThruMuTracker::makeTrackClusters3D( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
					  const std::vector< const BoundarySpacePoint* >& spacepts,
					  std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
					  std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices) {
    // This method takes in the list of boundaryspacepoints and pairs them up in order to try to find through-going muons
    // input:
    //   img_v: vector images, one for each plane. assumed to be in U,V,Y order
    //   spacepts: all possible boundary crossing points
    // output:
    //   trackclusters: thrumu tracks. collection of space points and pixel2dclusters
    //   tagged_v: pixels tagged as thrumu.  the value indicates the pass when the muon was tagged.
    //   used_endpoints_indices: indicates which endpoints in spacepts was used. 0=not used. 1=used.

    const int nendpts = (int)spacepts.size();
    used_endpoints_indices.resize( nendpts, 0 );

    // pair up containers
    int ntotsearched = 0;
    int npossible = 0;
    // poor-man's profiling
    const clock_t begin_time = clock();

    std::vector<int> tracks_per_pass;
    for (int ipass=0; ipass<m_config.num_passes; ipass++) {
      const ThruMuTrackerConfig::ThruMuPassConfig& passcfg = m_config.pass_configs.at(ipass);
      runPass( ipass, passcfg, spacepts, img_v, badchimg_v, tagged_v, used_endpoints_indices, trackclusters );
      int tracks_in_pass = trackclusters.size();
      if (ipass>0)
        tracks_in_pass -= tracks_per_pass.back();
      tracks_per_pass.push_back( tracks_in_pass );
    }

    if ( m_config.verbosity>0 ) {
      for (int ipass=0; ipass<m_config.num_passes; ipass++) {
        std::cout << "Number of tracks found in ThruMu pass #" << ipass << ": " << tracks_per_pass.at(ipass) << std::endl;
      }
    }

  }

  void ThruMuTracker::runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
			       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
			       std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters ) {

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


    std::cout << "Run Pass #" << passid << ": linear=" << passcfg.run_linear_tagger << " astar=" << passcfg.run_astar_tagger << std::endl;

    const int nendpts = (int)spacepts.size();
    std::vector<larlitecv::BMTrackCluster3D> pass_track_candidates;
    std::vector< std::vector<int> > pass_end_indices;

    for (int i=0; i<nendpts; i++) {
      if ( used_endpoints_indices.at(i) ) continue;
      const BoundarySpacePoint& pts_a = *(spacepts[i]);
      for (int j=i+1; j<nendpts; j++) {
        if ( used_endpoints_indices.at(j)) continue;
        const BoundarySpacePoint& pts_b = *(spacepts[j]);

        if ( pts_a.type()==pts_b.type() ) continue; // don't connect same type

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

        // first run the linear charge tagger, if so chosen
        if ( passcfg.run_linear_tagger ) {
          if ( m_config.verbosity>1 ) {
            std::cout << "  Running linear3dchargetagger." << std::endl;
          }
          BMTrackCluster3D linear_track = runLinearChargeTagger( passcfg, pts_a, pts_b, img_v, badchimg_v, linear_result );
          if ( m_config.verbosity>1 ) {
            std::cout << "  good fraction: " << linear_result.goodfrac << std::endl;
            std::cout << "  majority planes w/ charge fraction: " << linear_result.majfrac << std::endl;
          }
          if ( linear_result.isgood ) {
            if ( m_config.verbosity>1 ) std::cout << "  Result is good. length=" << linear_track.path3d.size() << std::endl;
            std::swap(track3d,linear_track);
          }
          else {
            if ( m_config.verbosity>1 ) std::cout << "  Result is bad." << std::endl;
          }
        }//end of if run_linear

        // next run astar tagger
        if ( passcfg.run_astar_tagger
             && (linear_result.goodfrac>passcfg.astar3d_min_goodfrac || linear_result.majfrac>passcfg.astar3d_min_majfrac ) ) {
          if ( m_config.verbosity>1 )
            std::cout << "  Running astar tagger." << std::endl;
          BMTrackCluster3D astar_track = runAStarTagger( passcfg, pts_a, pts_b, img_v, badchimg_v, astar_result );
          if ( astar_result.isgood ) {
            if ( m_config.verbosity>1 ) std::cout << "  Result is good." << std::endl;
            std::swap(track3d,astar_track);
          }
          else {
            if ( m_config.verbosity>1 ) std::cout << "  Result is bad." << std::endl;
          }
        }

        // run bezier fitter (doesn't exist yet. write this?)

        // could add criteria to filter here

        if ( (linear_result.isgood || astar_result.isgood) && m_config.verbosity>1 )
          std::cout << "  track found. length: " << track3d.path3d.size() << " empty=" << track3d.isempty() << std::endl;

        if ( !track3d.isempty() ) {
          std::vector<int> indices(2);
          indices[0] = i;
          indices[1] = j;
          pass_end_indices.emplace_back( std::move(indices) );
          pass_track_candidates.emplace_back( std::move(track3d) );
        }

      }// second loop over end points
    }// first loop over end points


    // post-process tracks (nothing yet)
    for ( auto &track : pass_track_candidates ) {
      trackclusters.emplace_back( std::move(track) );
    }
    for ( auto&indices : pass_end_indices ) {
      used_endpoints_indices.at(indices[0]) = 1;
      used_endpoints_indices.at(indices[1]) = 1;
    }
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
      && result_info.goodfrac > pass_cfg.linear3d_min_goodfraction
      && result_info.majfrac  > pass_cfg.linear3d_min_majoritychargefraction ) {
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
      for (int i=0; i<3; i++)
      	xyz[i] = ptinfo.xyz[i];
      path3d.emplace_back( std::move(xyz) );
    }

    larlitecv::BMTrackCluster3D track3d( pts_a, pts_b, path3d );
    return track3d;
  }

  larlitecv::BMTrackCluster3D ThruMuTracker::runAStarTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                  const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                  ThruMuTracker::AStarTaggerInfo& result_info ) {
    larlitecv::BMTrackCluster3D track3d;
    result_info.isgood = false;
    return track3d;
  }


  /*
    // post-processors. basically a filter for the tracks created.
    Linear3DPostProcessor linear_postprocess;

    // compress images for AStar3D
    std::vector< larcv::Image2D > img_compressed_v;
    std::vector< larcv::Image2D > badch_compressed_v;
    int downsampling_factor = 4;
    for (int p=0; p<3; p++) {
      larcv::Image2D img_compressed( img_v.at(p) );
      larcv::Image2D badch_compressed( badchimg_v.at(p) );
      img_compressed.compress( img_v.at(p).meta().rows()/downsampling_factor, img_v.at(p).meta().cols()/downsampling_factor );
      badch_compressed.compress( img_v.at(p).meta().rows()/downsampling_factor, img_v.at(p).meta().cols()/downsampling_factor );
      img_compressed_v.emplace_back( std::move(img_compressed) );
      badch_compressed_v.emplace_back( std::move(badch_compressed) );
    }

    // Do multiple passes in trying to make 3D tracks
    const int NumPasses = 2;
    // passes
    // (1) straight line fitter
    // (2) A* on remainder

    // make a vector of bools to tag end points that have already been used to find a track
    std::vector<bool> space_point_used( spacepts.size(), false );
    // we make blank images to tag pixels that have been assigned to a track
    if ( tagged_v.size()!=img_v.size() ) {
      tagged_v.clear();
      for (int p=0; p<(int)img_v.size(); p++) {
        larcv::Image2D tagged( img_v.at(p).meta() );
        tagged.paint(0.0);
        tagged_v.emplace_back( std::move(tagged) );
      }
    }

    // during pass 0 (line track), we store information about combinations of points
    // it will tell us if we should try astar during the next pass
    typedef std::pair<int,int> Combo_t;
    struct ComboInfo_t {
      float goodstart;
      float goodend;
      float fracGood;
      float fracMajCharge;
      int startext_majcharge;
      int endext_majcharge;
      int startext_ngood;
      int endext_ngood;
      bool track_made;
      ComboInfo_t() { goodstart = 0; goodend = 0; track_made = false; fracGood=0.; fracMajCharge=0.; };
    };
    typedef std::pair<Combo_t, ComboInfo_t> Pass0Pair;
    std::map< Combo_t, ComboInfo_t > pass0_combos;
    std::vector< Combo_t > pass0_endpts_connected;

    // do the track-finding passes
    for (int pass=0; pass<NumPasses; pass++) {

      // for each pass, we continue tag the images
      std::vector< larcv::Image2D > pass_tagged_v;
      for (size_t p=0; p<img_v.size(); p++ ) {
        larcv::Image2D tagged( img_v.at(p).meta() );
        tagged.paint(0.0);
        pass_tagged_v.emplace_back( std::move(tagged) );
      }
      // we compress tagged information from previous passes
      std::vector< larcv::Image2D > past_tagged_compressed_v;
      for (size_t p=0; p<tagged_v.size(); p++ ) {
        larcv::Image2D tagged_compressed( tagged_v.at(p) );
        tagged_compressed.compress( img_v.at(p).meta().rows()/downsampling_factor, img_v.at(p).meta().cols()/downsampling_factor );
        past_tagged_compressed_v.emplace_back( std::move(tagged_compressed) );
      }

      std::vector< BMTrackCluster3D > pass_tracks;

      for (int i=0; i<nendpts; i++) {
        //if ( space_point_used.at(i) ) continue; // a little too crude as control
        const BoundarySpacePoint& pts_a = *(spacepts[i]);
        for (int j=i+1; j<nendpts; j++) {
          //if ( space_point_used.at(j)) continue;
          const BoundarySpacePoint& pts_b = *(spacepts[j]);

          if ( pts_a.type()==pts_b.type() ) continue; // don't connect same type
          npossible++;

          if ( _config.verbosity>1 ) {
            std::cout << "[ Pass " << pass << ": path-finding for endpoints (" << i << "," << j << ") "
                      << "of type (" << pts_a.at(0).type << ") -> (" << pts_b.at(0).type << ") ]" << std::endl;

           for (int p=0; p<3; p++) {
              int row_a = pts_a.at(p).row;
              int col_a = pts_a.at(p).col;
              int col_b = pts_b.at(p).col;
              int row_b = pts_b.at(p).row;
              std::cout << "  plane=" << p << ": "
                        << " (w,t): (" << img_v.at(p).meta().pos_x( col_a ) << ", " << img_v.at(p).meta().pos_y( row_a ) << ") ->"
                        << " (" << img_v.at(p).meta().pos_x( col_b ) << "," << img_v.at(p).meta().pos_y( row_b ) << ")"
                        << std::endl;
            }
          }


          // don't try to connect points that are further apart in time than a full drift window
          bool within_drift = true;
          for (int p=0; p<3; p++) {
            int row_a = pts_a.at(p).row;
            int row_b = pts_b.at(p).row;
            if ( fabs( img_v.at(p).meta().pos_y( row_b )-img_v.at(p).meta().pos_y( row_a ) )>_config.ticks_per_full_drift ) { // ticks
              within_drift = false;
            }
          }
          if ( !within_drift ) {
            if ( _config.verbosity>1 )
              std::cout << "  time separation longer than drift window" << std::endl;
            continue;
          }

          // ====================================================================================
          // PASSES
          BMTrackCluster3D track3d;

          bool track_made = false;
          if ( pass==0 ) {
            // straight line fitter (Chris)
            std::vector<int> a_cols(img_v.size(),0);
            std::vector<int> b_cols(img_v.size(),0);
            for (int p=0; p<(int)img_v.size(); p++) {
              a_cols[p] = pts_a.at(p).col;
              b_cols[p] = pts_b.at(p).col;
            }

            PointInfoList straight_track = linetrackalgo.findpath( img_v, badchimg_v, pts_a.at(0).row, pts_b.at(0).row, a_cols, b_cols );

            // store info about the combination
            Combo_t thiscombo(i,j);
            ComboInfo_t badcombo;
            badcombo.goodstart = 0;
            badcombo.goodend   = 0;

            // we count the number of good points at start and end of track
            if ( straight_track.size()>0 ) {

              int nstart = straight_track.size()/4;
              int nend   = straight_track.size()/4;
              for (int ipt=0; ipt<nstart; ipt++) {
                if ( straight_track.at(ipt).planeswithcharge>=2 )
                  badcombo.goodstart+=1.0/float(nstart);
              }
             for (int ipt=(int)straight_track.size()-nend; ipt<(int)straight_track.size(); ipt++) {
                if ( straight_track.at(ipt).planeswithcharge>=2 )
                  badcombo.goodend += 1.0/float(nend);
              }
            }

            track_made = false;
            std::cout << "straight-track result: "
                      << "size=" << straight_track.size() << " "
                      << " fracGood=" << straight_track.fractionGood() << " "
                      << " fracAllCharge=" << straight_track.fractionHasChargeWith3Planes() << " "
                      << " fracMajCharge=" << straight_track.fractionHasChargeOnMajorityOfPlanes() << " "
                      << std::endl;
            badcombo.fracGood      = straight_track.fractionGood();
            badcombo.fracMajCharge = straight_track.fractionHasChargeOnMajorityOfPlanes();


            if ( (int)straight_track.size() > _config.linear3d_min_tracksize
                  && straight_track.fractionGood() > _config.linear3d_min_goodfraction
                  && straight_track.fractionHasChargeOnMajorityOfPlanes() > _config.linear3d_min_majoritychargefraction ) {
              std::cout << "  - straight-track accepted." << std::endl;
              // // cool, it is mostly a good track. let's check if end point
              // Endpoint checking didn't seem to contirbute much. We skip it for now.
              // PointInfoList start_ext;
              // PointInfoList end_ext;
              // linetrackalgo.getTrackExtension( straight_track, img_v, badchimg_v, 30.0, start_ext, end_ext );
              // //std::cout << "good track by linear fit. nstart=" << nstart << " nend=" << nend << std::endl;
              // badcombo.startext_majcharge = start_ext.num_pts_w_majcharge;
              // badcombo.endext_majcharge   = end_ext.num_pts_w_majcharge;
              // badcombo.startext_ngood     = start_ext.num_pts_good;
              // badcombo.endext_ngood       = end_ext.num_pts_good;

              // // use extensions to reject the track
              // if ( ( start_ext.num_pts_w_majcharge>=1000)
              //   || ( end_ext.num_pts_w_majcharge>=1000) ) {
              //   // this above cut has to be loose, because we often tag inside the track a bit
              //   // rejected
              //   track_made = false;
              // }
              // else {

              track_made = true;
              track3d = linetrackalgo.makeTrackCluster3D( img_v, badchimg_v, pts_a, pts_b, straight_track );
	      track3d.start_index = i;
	      track3d.end_index   = j;
              pass0_endpts_connected.push_back( thiscombo );
              //}

            }

            badcombo.track_made = track_made;
            pass0_combos.insert( Pass0Pair(thiscombo,badcombo) );

          }
          else if ( pass==1 ) {
            // A* star

            // because this can be an expensive algorithm we use test heuristic to see if we should run astar
            //bool shallwe = passTrackTest( pts_a, pts_b, img_v, badchimg_v );
            // for debugging specific tracks

            // We use information from the line tests to determine to run A*
            Combo_t thiscombo(i,j);

            // We already connect the two?
            auto it_combo = pass0_combos.find( thiscombo );
            if ( it_combo!=pass0_combos.end() ) {

              if ( (*it_combo).second.track_made )
              continue; // we have connected these points,  move on

              // end points much have way to travel to it
              if ( (*it_combo).second.goodend<0.10 || (*it_combo).second.goodstart<0.10 ) {
                if ( _config.verbosity>1 ) {
                  std::cout << "combo (" << i << "," << j << ") failed pass0 end heuristic "
                    << " fracstart=" << (*it_combo).second.goodstart << " "
                    << " fracend=" << (*it_combo).second.goodend
                    << std::endl;
                }
                continue;
              }

              // use pass0 calculation to determine if we run A*
              if ( (*it_combo).second.fracGood<_config.astar3d_min_goodfrac || (*it_combo).second.fracMajCharge<_config.astar3d_min_majfrac ) {
                if ( _config.verbosity>1)
                  std::cout << "failed pass0 heuristic [" << i << "," << j << "]: "
                            << " fracgood=" << (*it_combo).second.fracGood << " "
                            << " fracmajcharge=" << (*it_combo).second.fracMajCharge << " "
                            << std::endl;
                continue;
              }

            }
            else {
              // combo not found
              std::cout << "Pass0 combo not found." << std::endl;
              std::cin.get();
              continue;
            }

            if ( _config.verbosity>1)
                std::cout << "passed pass0 heuristic [" << i << "," << j << "]: "
                          << " fracgood=" << (*it_combo).second.fracGood << " "
                          << " fracmajcharge=" << (*it_combo).second.fracMajCharge << " "
                          << std::endl;

            // 3D A* path-finding
            bool goal_reached = false;
            track3d = runAstar3D( pts_a, pts_b, img_v, badchimg_v, tagged_v, img_compressed_v, badch_compressed_v, past_tagged_compressed_v, goal_reached );
	    track3d.start_index = i;
	    track3d.end_index   = j;
            std::vector< BMTrackCluster2D >& planetracks = track3d.plane_paths;

            if ( goal_reached ) {
              track_made = true;
            }
          } //end of pass 2 (a*)
          // ====================================================================================

          ntotsearched++;

          if ( track_made ) {
            // if we made a good track, we mark the end points as used. we also tag the path through the image
            track3d.markImageWithTrack( img_v, badchimg_v, _config.thresholds, _config.tag_neighborhood, pass_tagged_v );
            pass_tracks.emplace_back( std::move(track3d) );
          }
        }//loop over pt b
      }//loop over pt a

      // for combos where tracks were made, we remove them from the search if for some conditions
      for ( auto &combo : pass0_endpts_connected ) {

	// we tried to remove those tracks in the middle, but they were not easily discernable from those tagged some distance from the end
	// could explore a more conservative cut here in the future...

        // if we connected them. we no longer need the points
        space_point_used.at( combo.first )  = true;
        space_point_used.at( combo.second ) = true;
      }

      // merge tagged images from this pass to final tagged images
      for (size_t p=0; p<tagged_v.size(); p++ ) {
        larcv::Image2D& final_tagged = tagged_v.at(p);
        larcv::Image2D& pass_tagged  = pass_tagged_v.at(p);
        final_tagged += pass_tagged;
        final_tagged.binary_threshold( 1, 0, 255 );
      }

      // POST-PROCESS TRACKS
      if ( pass==0 ) {
	// linear post-processor
        std::vector< BMTrackCluster3D > filtered_tracks = linear_postprocess.process( pass_tracks );
        // fill the output tracks
        for ( auto &trk : filtered_tracks ) {
          trackclusters.emplace_back( std::move(trk) );
        }
      }
      else if ( pass==1 ) {
	// astar post-processor
        for ( auto &trk : pass_tracks )
          trackclusters.emplace_back( std::move(trk) );
      }

    }//end of loop over passes

    // boring book-keeping stuff ...
    for ( auto& track : trackclusters ) {
      used_endpoints_indices.at( track.start_index ) = 1;
      used_endpoints_indices.at( track.end_index ) = 1;
    }

    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "total paths searched: " << ntotsearched << " of " << npossible << " possible combinations. time=" << elapsed_secs << " secs" << std::endl;
    std::cout << "number of tracks created: " << trackclusters.size() << std::endl;

    return 0;
  }
  */


}
