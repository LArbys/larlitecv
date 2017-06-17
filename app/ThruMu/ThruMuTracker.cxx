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
      const ThruMuTrackerConfig::ThruMuPassConfig& passcfg = m_config.pass_configs.at(ipass);
      runPass( ipass, passcfg, spacepts, img_v, badchimg_v, tagged_v, used_endpoints_indices, trackclusters );
      int tracks_in_pass = trackclusters.size();
      if (ipass>0)
        tracks_in_pass -= tracks_per_pass.back();
      tracks_per_pass.push_back( tracks_in_pass );
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


}
