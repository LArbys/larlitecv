#include "BoundaryMuonTaggerAlgo.h"

// std
#include <vector>
#include <cmath>
#include <assert.h>
#include <ctime>

// larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "BasicTool/GeoAlgo/GeoAlgo.h"

// larcv
#include "UBWireTool/UBWireTool.h"

#include "Linear3DChargeTagger.h"
#include "Linear3DPostProcessor.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

namespace larlitecv {

  BoundaryMuonTaggerAlgo:: ~BoundaryMuonTaggerAlgo() {
    delete matchalgo_tight;
    delete matchalgo_loose;
  }

  void BoundaryMuonTaggerAlgo::run() {
    if (true)
      return;
    
    // this is still a work in progress. keeping the order of operations here as notes for now
    // eventually this will be the main function run by the user
    // searchforboundarypixels
    // clusterBoundaryPixels
  }

  int BoundaryMuonTaggerAlgo::searchforboundarypixels3D( const std::vector< larcv::Image2D >& imgs, 
                                                         const std::vector< larcv::Image2D >& badchs,
                                                         std::vector< BoundarySpacePoint >& end_points,
                                                         std::vector< larcv::Image2D >& matchedpixels,
                                                         std::vector< larcv::Image2D >& matchedspacepts ) {
    // this checks for pixels consistent with boundary crossings at the top/bottom/upstream/downstream portions of the TPC
    // returns an image that marks the location of such pixels
    // the vector returned has N*4 elements from N planes x 4 crossing types
    // the output images are in a different space:
    //  top/bottom: wire axis -> z-axis
    //  upstream/downstream: wire-axis -> y-axis
    // 
    // inputs
    // ------
    // imgs: ADC images
    // badchs: bad channel marked images
    //
    // outputs
    // -------
    // end_points: output spacepoints
    // matchedpixels:  marks ADC image where boundary tags found
    // matchedspacepts: marks in real space where boundary tags found

    TRandom3 rand( time(NULL) );
    //TRandom3 rand( 12345 );    

    std::cout << "BoundaryMuonTaggerAlgo::searchforboundarypixels3D" << std::endl;

    const int ncrossings = (int)BoundaryMuonTaggerAlgo::kNumCrossings; // 4, it's 4
    const int nplanes = (int)imgs.size();
    clock_t begin_time = clock();

    if ( !_config.checkOK() )  {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Not configured." << std::endl;
      return kErr_NotConfigured;
    }
    
    if ( nplanes<3 ) {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Expecting 3 planes currently." << std::endl;
      return kErr_BadInput;
    }

    // clear the output container
    matchedpixels.clear();
    matchedspacepts.clear();

    // get meta for input image
    const larcv::ImageMeta& meta = imgs.at(0).meta();

    // now we make dilated images
    std::vector<larcv::Image2D> binblur_v;
    int p=-1;
    for (auto const& img : imgs ) {
      p++;
      larcv::Image2D bbimg( img );
      bbimg.paint(0);
      for (int row=0; row<(int)img.meta().rows(); row++) {
	for (int col=0; col<(int)img.meta().cols(); col++) {
	  if ( img.pixel(row,col)>_config.thresholds[p] ) {
	    for (int dr=-_config.fKernelRadius; dr<=_config.fKernelRadius; dr++) {
	      int r = row+dr;
	      if ( r<0 || r>=(int)img.meta().rows() )
		continue;
	      for (int dc=-_config.fKernelRadius; dc<=_config.fKernelRadius; dc++) {
		int c = col+dc;
		if ( c<0 || c>=(int)img.meta().cols() )
		  continue;
		bbimg.set_pixel( r, c, 255 );
	      }
	    }
	  }
	}
      }//end of row loop
      binblur_v.emplace_back( std::move(bbimg) );
    }//end of img loop
    
    // create Image2D objects to store the result of boundary matching
    // we use meta to make the same sized image as the input
    if ( _config.save_endpt_images ) {
      for (int i=0; i<ncrossings; i++) { // top, bottom, upstream, downstream x 3 planes
        larcv::Image2D matchimage( meta );
        matchimage.paint(0.0);
        matchedspacepts.emplace_back( std::move(matchimage) ); // uses the same meta... [I should control this. this is fragile as is.]
        for (int p=0; p<nplanes; p++) {
          larcv::Image2D matchimage2( imgs.at(p).meta() );
          matchimage2.paint(0.0);
          matchedpixels.emplace_back( std::move(matchimage2) );
        }
      }
    }
    
    // storage for boundary combination
    clock_t begin_pixel_search = clock();
    std::cout << "Begin Boundary Pixel Search..." << std::endl;
    std::vector< dbscan::dbPoints > combo_points(ncrossings); // these points are in detector space
    std::vector< std::vector< std::vector<int> > > combo_cols(ncrossings); // [ncrossings][number of combos][column combination]
    CollectCandidateBoundaryPixels( binblur_v, badchs, combo_points, combo_cols, matchedspacepts );
    int total_combos = 0;
    for ( auto& combo_col : combo_cols ) {
      total_combos += combo_col.size();
    }
    float elapsed_hitsearch = float( clock()-begin_pixel_search )/CLOCKS_PER_SEC;
    std::cout << "... hit search time: " << elapsed_hitsearch << " secs" << std::endl;

    // cluster each boundary type pixel (cluster of points in detector space)
    clock_t begin_clustering = clock();
    std::cout << "Begin Clustering..." << std::endl;
    std::vector< BoundarySpacePoint > candidate_endpts;
    ClusterBoundaryHitsIntoEndpointCandidates( imgs, badchs, combo_points, combo_cols, candidate_endpts, matchedpixels );
    float elapsed_clustering = float( clock()-begin_clustering )/CLOCKS_PER_SEC;
    std::cout << "... clustering time: " << elapsed_clustering << " secs" << std::endl;


    for ( int endpt_idx=0; endpt_idx<(int)candidate_endpts.size(); endpt_idx++ ) {
      end_points.emplace_back( std::move( candidate_endpts.at(endpt_idx) ) );
    }

    // generate the meta data we will use to filter end points
    //std::vector< std::vector<ClusterExtrema_t> > candidate_metadata; 
    //int rmax_window = 10;
    //int rmin_window = 10;
    //int col_window  = 10;
    //GenerateEndPointMetaData( candidate_endpts, imgs, rmax_window, rmin_window, col_window, candidate_metadata );

    //std::vector<int> passes_check;
    //CheckClusterExtrema(  candidate_endpts, candidate_metadata, imgs, passes_check);

    //std::vector< int > cluster_passed;
    //SelectOnlyTrackEnds( candidate_endpts, imgs, 10, 10, 20, cluster_passed );

    //for ( int endpt_idx=0; endpt_idx<(int)candidate_endpts.size(); endpt_idx++ ) {
    //  if ( cluster_passed[endpt_idx]>=0 )
    //    end_points.emplace_back( std::move( candidate_endpts.at(endpt_idx) ) );
    //}
    
    
    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "boundary pixel search took " << elapsed_secs << " secs" << std::endl;
    std::cout << "  hit collecting: " << elapsed_hitsearch << " secs" << std::endl;
    std::cout << "  clustering time: " << elapsed_clustering << " secs" << std::endl;
    std::cout << "  total number of combos found: " << total_combos << std::endl;

    
    return kOK;
  }

  void BoundaryMuonTaggerAlgo::loadGeoInfo() {

    TFile fGeoFile( Form("%s/app/PMTWeights/dat/geoinfo.root",getenv("LARCV_BASEDIR")), "OPEN" );

    // Get the PMT Info
    fNPMTs = 32;
    TTree* fPMTTree  = (TTree*)fGeoFile.Get( "imagedivider/pmtInfo" );
    int femch;
    float pos[3];
    fPMTTree->SetBranchAddress( "femch", &femch );
    fPMTTree->SetBranchAddress( "pos", pos );
    for (int n=0; n<fNPMTs; n++) {
      fPMTTree->GetEntry(n);
      for (int i=0; i<3; i++) {
        pmtpos[femch][i] = pos[i];
      }
      //std::cout << "[POS " << femch << "] " << " (" << pmtpos[femch][0] << "," << pmtpos[femch][1] << "," << pmtpos[femch][2] << ")" << std::endl;
    }
    
    fGeoFile.Close();
    
    matchalgo_tight = new larlitecv::BoundaryMatchAlgo( larlitecv::BoundaryMatchArrays::kTight );
    matchalgo_loose = new larlitecv::BoundaryMatchAlgo( larlitecv::BoundaryMatchArrays::kLoose );
    
  }
  
//   int BoundaryMuonTaggerAlgo::makeTrackClusters3D( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
//                                                    const std::vector< const BoundarySpacePoint* >& spacepts,
//                                                    std::vector< larlitecv::BMTrackCluster3D >& trackclusters, 
//                                                    std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices, const std::vector <larlite::event_opflash>& opflash_v ) {
//     // This method takes in the list of boundaryspacepoints and pairs them up in order to try to find through-going muons
//     int nendpts = (int)spacepts.size();
//     used_endpoints_indices.resize( nendpts, 0 );

//     // pair up containers
//     int ntotsearched = 0;
//     int npossible = 0;
//     // poor-man's profiling
//     const clock_t begin_time = clock();
//     std::vector<int> modimg(nendpts,0);

//     // linear 3D track
//     Linear3DChargeTagger linetrackalgo( _config.linear3d_cfg ); // finds charge along a line in 3D space

//     // post-processors. basically a filter for the tracks created.
//     Linear3DPostProcessor linear_postprocess;

//     // compress images for AStar3D
//     std::vector< larcv::Image2D > img_compressed_v;
//     std::vector< larcv::Image2D > badch_compressed_v;
//     int downsampling_factor = 4;
//     for (int p=0; p<3; p++) {
//       larcv::Image2D img_compressed( img_v.at(p) );
//       larcv::Image2D badch_compressed( badchimg_v.at(p) );
//       img_compressed.compress( img_v.at(p).meta().rows()/downsampling_factor, img_v.at(p).meta().cols()/downsampling_factor );
//       badch_compressed.compress( img_v.at(p).meta().rows()/downsampling_factor, img_v.at(p).meta().cols()/downsampling_factor );
//       img_compressed_v.emplace_back( std::move(img_compressed) );
//       badch_compressed_v.emplace_back( std::move(badch_compressed) );
//     }

//     // Do multiple passes in trying to make 3D tracks
//     const int NumPasses = 2;
//     // passes
//     // (1) straight line fitter
//     // (2) A* on remainder

//     // make a vector of bools to tag end points that have already been used to find a track
//     std::vector<bool> space_point_used( spacepts.size(), false );
//     // we make blank images to tag pixels that have been assigned to a track
//     if ( tagged_v.size()!=img_v.size() ) {
//       tagged_v.clear();
//       for (int p=0; p<(int)img_v.size(); p++) {
//         larcv::Image2D tagged( img_v.at(p).meta() );
//         tagged.paint(0.0);
//         tagged_v.emplace_back( std::move(tagged) );
//       }
//     }

//     // during pass 0 (line track), we store information about combinations of points
//     // it will tell us if we should try astar during the next pass
//     typedef std::pair<int,int> Combo_t;
//     struct ComboInfo_t {
//       float goodstart;
//       float goodend;
//       float fracGood;
//       float fracMajCharge;
//       int startext_majcharge;
//       int endext_majcharge;
//       int startext_ngood;
//       int endext_ngood;
//       bool track_made;
//       ComboInfo_t() { goodstart = 0; goodend = 0; track_made = false; fracGood=0.; fracMajCharge=0.; };
//     };
//     typedef std::pair<Combo_t, ComboInfo_t> Pass0Pair;
//     std::map< Combo_t, ComboInfo_t > pass0_combos;
//     std::vector< Combo_t > pass0_endpts_connected;

//     // do the track-finding passes
//     for (int pass=0; pass<NumPasses; pass++) {

//       // for each pass, we continue tag the images
//       std::vector< larcv::Image2D > pass_tagged_v;
//       for (size_t p=0; p<img_v.size(); p++ ) {
//         larcv::Image2D tagged( img_v.at(p).meta() );
//         tagged.paint(0.0);
//         pass_tagged_v.emplace_back( std::move(tagged) );
//       }
//       // we compress tagged information from previous passes
//       std::vector< larcv::Image2D > past_tagged_compressed_v; 
//       for (size_t p=0; p<tagged_v.size(); p++ ) {
//         larcv::Image2D tagged_compressed( tagged_v.at(p) );
//         tagged_compressed.compress( img_v.at(p).meta().rows()/downsampling_factor, img_v.at(p).meta().cols()/downsampling_factor );
//         past_tagged_compressed_v.emplace_back( std::move(tagged_compressed) );
//       }

//       std::vector< BMTrackCluster3D > pass_tracks;      

//       for (int i=0; i<nendpts; i++) {
//         //if ( space_point_used.at(i) ) continue; // a little too crude as control
//         const BoundarySpacePoint& pts_a = *(spacepts[i]);        
//         for (int j=i+1; j<nendpts; j++) {
//           //if ( space_point_used.at(j)) continue;
//           const BoundarySpacePoint& pts_b = *(spacepts[j]);

//           if ( pts_a.type()==pts_b.type() ) continue; // don't connect same type
//           npossible++;
          
//           if ( _config.verbosity>1 ) {
//             std::cout << "[ Pass " << pass << ": path-finding for endpoints (" << i << "," << j << ") "
//                       << "of type (" << pts_a.at(0).type << ") -> (" << pts_b.at(0).type << ") ]" << std::endl;

//            for (int p=0; p<3; p++) {
//               int row_a = pts_a.at(p).row;
//               int col_a = pts_a.at(p).col;
//               int col_b = pts_b.at(p).col;
//               int row_b = pts_b.at(p).row;
//               std::cout << "  plane=" << p << ": "
//                         << " (w,t): (" << img_v.at(p).meta().pos_x( col_a ) << ", " << img_v.at(p).meta().pos_y( row_a ) << ") ->"
//                         << " (" << img_v.at(p).meta().pos_x( col_b ) << "," << img_v.at(p).meta().pos_y( row_b ) << ")"
//                         << std::endl;           
//             }
//           }
        

//           // don't try to connect points that are further apart in time than a full drift window
//           bool within_drift = true;
//           for (int p=0; p<3; p++) {
//             int row_a = pts_a.at(p).row;
//             int row_b = pts_b.at(p).row;
//             if ( fabs( img_v.at(p).meta().pos_y( row_b )-img_v.at(p).meta().pos_y( row_a ) )>_config.ticks_per_full_drift ) { // ticks
//               within_drift = false;
//             }
//           }
//           if ( !within_drift ) {
//             if ( _config.verbosity>1 )
//               std::cout << "  time separation longer than drift window" << std::endl;
//             continue;
//           }

//           // ====================================================================================
//           // PASSES
//           BMTrackCluster3D track3d;

//           bool track_made = false;
//           if ( pass==0 ) {
//             // straight line fitter (Chris)
//             std::vector<int> a_cols(img_v.size(),0);
//             std::vector<int> b_cols(img_v.size(),0);
//             for (int p=0; p<(int)img_v.size(); p++) {
//               a_cols[p] = pts_a.at(p).col;
//               b_cols[p] = pts_b.at(p).col;
//             }

//             PointInfoList straight_track = linetrackalgo.findpath( img_v, badchimg_v, pts_a.at(0).row, pts_b.at(0).row, a_cols, b_cols );

//             // store info about the combination
//             Combo_t thiscombo(i,j);
//             ComboInfo_t badcombo;
//             badcombo.goodstart = 0;
//             badcombo.goodend   = 0;

//             // we count the number of good points at start and end of track
//             if ( straight_track.size()>0 ) {

//               int nstart = straight_track.size()/4;
//               int nend   = straight_track.size()/4;
//               for (int ipt=0; ipt<nstart; ipt++) {
//                 if ( straight_track.at(ipt).planeswithcharge>=2 )
//                   badcombo.goodstart+=1.0/float(nstart);
//               }
//              for (int ipt=(int)straight_track.size()-nend; ipt<(int)straight_track.size(); ipt++) {
//                 if ( straight_track.at(ipt).planeswithcharge>=2 )
//                   badcombo.goodend += 1.0/float(nend);
//               }
//             }

//             track_made = false;
//             std::cout << "straight-track result: "
//                       << "size=" << straight_track.size() << " "
//                       << " fracGood=" << straight_track.fractionGood() << " "
//                       << " fracAllCharge=" << straight_track.fractionHasChargeWith3Planes() << " "
//                       << " fracMajCharge=" << straight_track.fractionHasChargeOnMajorityOfPlanes() << " "
//                       << std::endl;
//             badcombo.fracGood      = straight_track.fractionGood();
//             badcombo.fracMajCharge = straight_track.fractionHasChargeOnMajorityOfPlanes();


//             if ( (int)straight_track.size() > _config.linear3d_min_tracksize
//                   && straight_track.fractionGood() > _config.linear3d_min_goodfraction
//                   && straight_track.fractionHasChargeOnMajorityOfPlanes() > _config.linear3d_min_majoritychargefraction ) {
//               std::cout << "  - straight-track accepted." << std::endl;
//               // // cool, it is mostly a good track. let's check if end point
//               // Endpoint checking didn't seem to contirbute much. We skip it for now.
//               // PointInfoList start_ext;
//               // PointInfoList end_ext;
//               // linetrackalgo.getTrackExtension( straight_track, img_v, badchimg_v, 30.0, start_ext, end_ext );
//               // //std::cout << "good track by linear fit. nstart=" << nstart << " nend=" << nend << std::endl;
//               // badcombo.startext_majcharge = start_ext.num_pts_w_majcharge;
//               // badcombo.endext_majcharge   = end_ext.num_pts_w_majcharge;
//               // badcombo.startext_ngood     = start_ext.num_pts_good;
//               // badcombo.endext_ngood       = end_ext.num_pts_good;

//               // // use extensions to reject the track
//               // if ( ( start_ext.num_pts_w_majcharge>=1000)
//               //   || ( end_ext.num_pts_w_majcharge>=1000) ) {
//               //   // this above cut has to be loose, because we often tag inside the track a bit
//               //   // rejected
//               //   track_made = false;
//               // }
//               // else {

//               track_made = true;
//               track3d = linetrackalgo.makeTrackCluster3D( img_v, badchimg_v, pts_a, pts_b, straight_track );
// 	      track3d.start_index = i;
// 	      track3d.end_index   = j;
//               pass0_endpts_connected.push_back( thiscombo );
//               //}

//             }

//             badcombo.track_made = track_made;
//             pass0_combos.insert( Pass0Pair(thiscombo,badcombo) );

//           }
//           else if ( pass==1 ) {
//             // A* star

//             // because this can be an expensive algorithm we use test heuristic to see if we should run astar
//             //bool shallwe = passTrackTest( pts_a, pts_b, img_v, badchimg_v );
//             // for debugging specific tracks

//             // We use information from the line tests to determine to run A*
//             Combo_t thiscombo(i,j);

//             // We already connect the two?
//             auto it_combo = pass0_combos.find( thiscombo );
//             if ( it_combo!=pass0_combos.end() ) {

//               if ( (*it_combo).second.track_made )
//               continue; // we have connected these points,  move on

//               // end points much have way to travel to it
//               if ( (*it_combo).second.goodend<0.10 || (*it_combo).second.goodstart<0.10 ) {
//                 if ( _config.verbosity>1 ) {
//                   std::cout << "combo (" << i << "," << j << ") failed pass0 end heuristic "
//                     << " fracstart=" << (*it_combo).second.goodstart << " "
//                     << " fracend=" << (*it_combo).second.goodend
//                     << std::endl;
//                 }
//                 continue;
//               }

//               // use pass0 calculation to determine if we run A*
//               if ( (*it_combo).second.fracGood<_config.astar3d_min_goodfrac || (*it_combo).second.fracMajCharge<_config.astar3d_min_majfrac ) {
//                 if ( _config.verbosity>1)
//                   std::cout << "failed pass0 heuristic [" << i << "," << j << "]: "
//                             << " fracgood=" << (*it_combo).second.fracGood << " "
//                             << " fracmajcharge=" << (*it_combo).second.fracMajCharge << " "                            
//                             << std::endl;
//                 continue;
//               }

//             }
//             else {
//               // combo not found
//               std::cout << "Pass0 combo not found." << std::endl;
//               std::cin.get();
//               continue;
//             }

//             if ( _config.verbosity>1)
//                 std::cout << "passed pass0 heuristic [" << i << "," << j << "]: "
//                           << " fracgood=" << (*it_combo).second.fracGood << " "
//                           << " fracmajcharge=" << (*it_combo).second.fracMajCharge << " "                            
//                           << std::endl;
            
//             // 3D A* path-finding
//             bool goal_reached = false;
//             track3d = runAstar3D( pts_a, pts_b, img_v, badchimg_v, tagged_v, img_compressed_v, badch_compressed_v, past_tagged_compressed_v, goal_reached );
// 	    track3d.start_index = i;
// 	    track3d.end_index   = j;
//             std::vector< BMTrackCluster2D >& planetracks = track3d.plane_paths;

//             if ( goal_reached ) {
//               track_made = true;
//             }
//           } //end of pass 2 (a*)
//           // ====================================================================================

//           ntotsearched++;          

//           if ( track_made ) {
//             // if we made a good track, we mark the end points as used. we also tag the path through the image
//             track3d.markImageWithTrack( img_v, badchimg_v, _config.thresholds, _config.tag_neighborhood, pass_tagged_v );
//             pass_tracks.emplace_back( std::move(track3d) );
//           }
//         }//loop over pt b
//       }//loop over pt a

//       // for combos where tracks were made, we remove them from the search if for some conditions
//       for ( auto &combo : pass0_endpts_connected ) {

// 	// we tried to remove those tracks in the middle, but they were not easily discernable from those tagged some distance from the end
// 	// could explore a more conservative cut here in the future...
	
//         // if we connected them. we no longer need the points
//         space_point_used.at( combo.first )  = true;
//         space_point_used.at( combo.second ) = true;        
//       }

//       // merge tagged images from this pass to final tagged images
//       for (size_t p=0; p<tagged_v.size(); p++ ) {
//         larcv::Image2D& final_tagged = tagged_v.at(p);
//         larcv::Image2D& pass_tagged  = pass_tagged_v.at(p);
//         final_tagged += pass_tagged;
//         final_tagged.binary_threshold( 1, 0, 255 );
//       }

//       // POST-PROCESS TRACKS
//       if ( pass==0 ) {
// 	// linear post-processor
//         std::vector< BMTrackCluster3D > filtered_tracks = linear_postprocess.process( pass_tracks );
//         // fill the output tracks
//         for ( auto &trk : filtered_tracks ) {
//           trackclusters.emplace_back( std::move(trk) );
//         }
//       }
//       else if ( pass==1 ) {
// 	// astar post-processor
//         for ( auto &trk : pass_tracks )
//           trackclusters.emplace_back( std::move(trk) );
//       }

//     }//end of loop over passes

//     // boring book-keeping stuff ...
//     for ( auto& track : trackclusters ) {
//       used_endpoints_indices.at( track.start_index ) = 1;
//       used_endpoints_indices.at( track.end_index ) = 1;      
//     }
    
//     float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
//     std::cout << "total paths searched: " << ntotsearched << " of " << npossible << " possible combinations. time=" << elapsed_secs << " secs" << std::endl;
//     std::cout << "number of tracks created: " << trackclusters.size() << std::endl;
    
//     return 0;
//   }

//   int BoundaryMuonTaggerAlgo::markImageWithTrackClusters( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
//                                                           const std::vector< larlitecv::BMTrackCluster3D >& tracks,
//                                                           std::vector<int>& goodlist, std::vector<larcv::Image2D>& markedimgs ) {
//     for ( int itrack=0; itrack<(int)tracks.size(); itrack++ ) {
//       if ( goodlist.at(itrack)==false ) continue;
//       tracks.at(itrack).markImageWithTrack( imgs, badchimgs, _config.thresholds, _config.tag_neighborhood, markedimgs );
//     }
    
//     return 0;
//   }
 
//   BMTrackCluster3D BoundaryMuonTaggerAlgo::runAstar3D( const BoundarySpacePoint& start_pt, const BoundarySpacePoint& end_pt, 
//         const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
//         const std::vector<larcv::Image2D>& img_compressed_v, const std::vector<larcv::Image2D>& badch_compressed_v, const std::vector<larcv::Image2D>& tagged_compressed_v,
//         bool& goal_reached ) {

//     // This wraps/interfaces with the AStar algorithm on our image. runs one start to end check.
//     // inputs
//     // ------
//     //  start: starting boundary point
//     //  end: ending boundary point
//     //  img: image with charge deposition we will follow
//     //  start_pad: radius around starting point where we will include all hits in case not start point directly on track
//     //  end_pad:   radius around ending point where we will include all hits in case end point not directly on track

//     BMTrackCluster3D track3d;

//     // Run AStar3DAlgo

//     // configure
//     const larlitecv::AStar3DAlgoConfig& astar_config = _config.getAStarConfig();
//     /*
//     astar_config.astar_threshold.resize(3,0);
//     astar_config.astar_threshold[0] = 100.0;
//     astar_config.astar_threshold[1] = 100.0;
//     astar_config.astar_threshold[2] = 150.0;
//     astar_config.astar_neighborhood.resize(3,6);
//     astar_config.astar_start_padding = 3;
//     astar_config.astar_end_padding   = 3;
//     astar_config.lattice_padding = 10;
//     astar_config.accept_badch_nodes = true;
//     astar_config.min_nplanes_w_hitpixel = 3;
//     */

//     // collect meta/translate start/goal tick/wires to row/col in compressed image
//     std::vector< const larcv::ImageMeta* > meta_compressed_v;
//     std::vector<int> start_cols(img_compressed_v.size(),0);
//     std::vector<int> start_rows(img_compressed_v.size(),0);
//     std::vector<int> goal_cols(img_compressed_v.size(),0);
//     std::vector<int> goal_rows(img_compressed_v.size(),0);    
//     std::vector< const larcv::Pixel2D* > start_pix;
//     std::vector< const larcv::Pixel2D* > goal_pix;

//     for (size_t p=0; p<img_compressed_v.size(); p++) {

//       // get start/end point informatio in compressed image
//       const larcv::ImageMeta* ptr_meta = &(img_compressed_v.at(p).meta());

//       const larlitecv::BoundaryEndPt& start_endpt = start_pt.at(p);
//       start_rows[p] =  ptr_meta->row( img_v.at(p).meta().pos_y( start_endpt.row ) );
//       start_cols[p] =  ptr_meta->col( img_v.at(p).meta().pos_x( start_endpt.col ) );
//       start_pix.push_back( new larcv::Pixel2D( start_endpt.getcol(), start_endpt.getrow() ) );

//       const larlitecv::BoundaryEndPt& goal_endpt  = end_pt.at(p);
//       goal_rows[p]  =  ptr_meta->row( img_v.at(p).meta().pos_y( goal_endpt.row ) );
//       goal_cols[p]  =  ptr_meta->col( img_v.at(p).meta().pos_x( goal_endpt.col ) );
//       goal_pix.push_back( new larcv::Pixel2D( goal_endpt.getcol(), goal_endpt.getrow() ) );      

//       meta_compressed_v.push_back( ptr_meta );
//     }


//     larlitecv::AStar3DAlgo algo( astar_config );
//     std::vector<larlitecv::AStar3DNode> path;
//     goal_reached = false;
//     int goalhit = 0;
//     try {
//       path = algo.findpath( img_compressed_v, badch_compressed_v, tagged_compressed_v,
//         start_rows.front(), goal_rows.front(), start_cols, goal_cols, goalhit );
//       if ( goalhit==1 )
//         goal_reached = true;
//     }
//     catch (const std::exception& e) {
//       std::cout << "exception running astar3dalgo::findpath: " << e.what() << std::endl;
//       BMTrackCluster3D track3d;
//       for (int p=0; p<3; p++ ) {
//         delete start_pix.at(p);
//         delete goal_pix.at(p);
//       }
//       return track3d; // return empty track
//     }

//     float nbad_nodes = 0;
//     float total_nodes = 0;
//     for ( auto& node : path ) {
//       if ( node.badchnode )
//         nbad_nodes+=1.0;
//       total_nodes+=1.0;
//     }

//     // if majority are bad ch nodes, reject this track
//     if ( nbad_nodes/total_nodes>0.5 )
//       goal_reached = false;

//     track3d = AStarNodes2BMTrackCluster3D( path, img_v, badch_v, start_pix, goal_pix, 
//       _config.tag_neighborhood[0], _config.track_step_size, _config.thresholds );
//     for (int p=0; p<3; p++ ) {
//       delete start_pix.at(p);
//       delete goal_pix.at(p);
//     }    

//     return track3d;
//   }

//   void BoundaryMuonTaggerAlgo::process2Dtracks( std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters2D,
//                                                 const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
//                                                 std::vector< BMTrackCluster3D >& tracks, std::vector<int>& goodlist ) {
//     // DEPRECATED 
//     std::vector< std::vector<float> > valid_range(3);
//     for (int p=0; p<3; p++) {
//       valid_range[p].resize(2);
//       valid_range[p][0] =   -20.;
//       valid_range[p][1] = 1100.0;
//     }

//     goodlist.resize( trackclusters2D.size(), 0 );

//     //float trig_tick = 2400;
//     //float drift_v   = 0.106865;

//     // loop over 
//     int itrack = -1;
//     for ( auto &track2d : trackclusters2D ) {
//       itrack++;
//       std::cout << "=============================================" << std::endl;
//       std::cout << "PROCESSING TRACK 2D #" << itrack << std::endl;
//       for (int p=0; p<3; p++)
//         std::cout << "  plane " << p << " start=(" << track2d.at(p).start.col << "," << track2d.at(p).start.row << ") ";
//       std::cout << std::endl;
//       for (int p=0; p<3; p++)
//         std::cout << "  plane " << p << " end=(" << track2d.at(p).end.col << "," << track2d.at(p).end.row << ") ";
//       std::cout << std::endl;
          

//       bool issame = false;
//       if ( tracks.size()>0 ) {
//         // we check to see if the current track is basically on the same path as a previous track
//         for ( auto const& track3d : tracks ) {
//           issame = compare2Dtrack( track2d, track3d, img_v.at(0).meta(), 5.0, 5.0 );
//           if ( issame ) {
//             break;
//           }
//         }
//       }
//       if ( issame ) continue; // to next 2d track

//       // if different, we process track into 3D
//       BMTrackCluster3D track3d = process2Dtrack( track2d, img_v, badchimg_v );
//       if ( track3d.path3d.size()==0 ) {
//         std::cout << "error producing this track" << std::endl;
//         continue; // onto next 2d track
//       }

//       // save track
//       track3d.track2d_index = itrack;
//       tracks.emplace_back( std::move(track3d) );
//       goodlist.at(itrack) = 1; // mark as good

//     }//end of loop over 2D tracks
//   }
  

//   BMTrackCluster3D BoundaryMuonTaggerAlgo::process2Dtrack( std::vector< larlitecv::BMTrackCluster2D >& track2d, 
//                                                            const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v ) {
    
//     // parameters to move to config file at some point
//     std::vector< std::vector<float> > valid_range(2);
//     valid_range[0].resize(2);
//     valid_range[1].resize(2);
//     valid_range[0][0] = -100;
//     valid_range[0][1] = 1200;
//     valid_range[1][0] = -150.0;
//     valid_range[1][1] =  150.0;
    
//     float trig_tick = 2400;
//     float drift_v   = 0.106865;
    
//     // we find path length, number of nodes, path, breakdown for each 2d track
//     // we also look for dumb tracks that have crazy kinks
//     float pathlength[3];
//     std::vector< float > edgelength[3];
//     std::vector< std::vector<float> > edgedir[3];
//     std::vector< std::vector<float> > nodepos[3];
//     bool isgood[3] = { true, true, true };
    
//     BMTrackCluster3D track3d;

//     for (size_t p=0; p<3; p++) {
//       const BMTrackCluster2D& planetrack = track2d.at(p);
//       pathlength[p] = 0.0;
//       if ( planetrack.pixelpath.size()<2 ) {
//         isgood[p] = false;
//         //std::cout << "p=" << p << " does not have enough pixels" << std::endl;
//         continue; // skip this plane
//       }
//       // follow the 2d track paths and file path info variables above
//       for (int i=1; i<(int)planetrack.pixelpath.size(); i++) {

//         std::vector<float> pos1(2,0.0);
//         pos1[0] = planetrack.pixelpath.at(i).X();
//         pos1[1] = planetrack.pixelpath.at(i).Y();

//         std::vector<float> pos0(2,0.0);
//         pos0[0] = planetrack.pixelpath.at(i-1).X();
//         pos0[1] = planetrack.pixelpath.at(i-1).Y();

//         float dx = pos1[0]-pos0[0];
//         float dy = pos1[1]-pos0[1];
//         float dist = sqrt( dx*dx + dy*dy );

//         std::vector<float> ndir(2,0.0);
//         ndir[0] = dx/dist;
//         ndir[1] = dy/dist;
        
// //      std::cout << "pix " << i << "(" << pos0[0] << "," << pos0[1] << ") -> (" << pos1[0] << "," << pos1[1] << ")"
// //                << ". (dx,dy)=(" << dx << "," << dy << ") dist=" << dist << std::endl;
        
//         if ( i>=2 ) {
//           std::vector<float> lastdir  = edgedir[p].at(i-2);
//           float lastcos = lastdir[0]*ndir[0] + ndir[1]*lastdir[1];
//           if ( lastcos<0.70 ) { // around 45 degrees
//             // kink. probably bad
//             // std::cout << "kink found: p=" << p << " lastcos=" << lastcos << " step=" << i << ":"
//             //        << "(" << planetrack.pixelpath.at(i-1).X() << "," << planetrack.pixelpath.at(i-1).Y() << ") -> "
//             //        << "(" << planetrack.pixelpath.at(i).X() << ","<< planetrack.pixelpath.at(i).Y() << "). "
//             //        << " dir: "
//             //        << "(" << lastdir[0] << ", " << lastdir[1] << ") -> "
//             //        << "(" << ndir[0] << "," << ndir[1] << ")" << std::endl;
//             //isgood[p] = false;
//           }
//         }
//         edgelength[p].push_back( dist );
//         pathlength[p] += dist;
//         edgedir[p].emplace_back( std::move(ndir) );
//         nodepos[p].emplace_back( std::move(pos1) );
//       }//end of loop over path
//       // add end point
//       std::vector<float> pos(2,0.0);
//       pos[0] = planetrack.pixelpath.back().X();
//       pos[1] = planetrack.pixelpath.back().Y();
//       nodepos[p].emplace_back( std::move(pos) );
//       edgelength[p].push_back(0);
//       edgedir[p].emplace_back( std::vector<float>(2,0.0) );
//     }//end of loop over planes

    
//     // use path variables to build 3D model
//     int ngood = 0;
//     float longest_pathlength = 0;
//     for (int p=0; p<3; p++) {
//       if ( isgood[p] ) { 
//         ngood++;
//         if ( pathlength[p]>longest_pathlength ) {
//           longest_pathlength = pathlength[p];
//         }
//       }
//     }

//     if ( ngood<2 ) 
//       return track3d; // empty track

//     int badplane = -1;
//     if ( ngood==2 ) {
//       // clear out path of BMTrackCluster2D from bad plane. we are going to remake it
//       for (int p=0; p<3; p++) {
//         if (!isgood[p]) {
//           BMTrackCluster2D& badtrack = track2d.at(p);
//           badtrack.pixelpath.clear();
//           larcv::Pixel2D pixel( badtrack.start.col, badtrack.start.row ); // (X,Y)
//           //badtrack.pixelpath.emplace_back( std::move(pixel) );
//           badplane = p;
//           break;
//         }
//       }
//     }
    
//     // make 3d path. we want to take cm steps
//     int nsteps = (int)longest_pathlength;
//     int current_node[3] = { 0, 0, 0 };
//     float node_ds[3] = {0, 0, 0 };
          
//     // std::cout << "Track with " << nsteps << " steps. number of nodes="
//     //        << "(" << nodepos[0].size() << "," << nodepos[1].size() << "," << nodepos[2].size() << ")"
//     //        << std::endl;
//     for (int istep=0; istep<nsteps; istep++) {
        
//       // get wire at this step
//       int imgcol[3] = { -1, -1, -1 };
//       int wireid[3] = { -1, -1, -1 };
//       int pid[3] = { -1, -1, -1 };
//       int ip=0;
//       float avetick = 0.0;
//       for (int p=0; p<3; p++) {
//         if ( !isgood[p] ) {
//           continue;
//         }
//         float imgpos[2] = { 0., 0. };
//         for (int v=0; v<2; v++) {
//           imgpos[v] = nodepos[p].at( current_node[p] )[v] + edgedir[p].at( current_node[p] )[v]*node_ds[v];
//         }
//         imgcol[ip] = (int)imgpos[0];
//         wireid[ip] = (int)imgcol[ip]*img_v.at(p).meta().pixel_width();
//         if ( wireid[ip]<0 ) wireid[ip]=0;
//         if ( wireid[ip]>=(int)larutil::Geometry::GetME()->Nwires(p) ) wireid[ip] = (int)larutil::Geometry::GetME()->Nwires(p)-1;
//         pid[ip] = p;
//         int imgpos_row = (int)imgpos[1];
//         if ( imgpos_row<0 ) imgpos_row = 0;
//         if ( imgpos_row>=(int)img_v.at(p).meta().rows() ) imgpos_row = (int)img_v.at(p).meta().rows()-1;
//         avetick += img_v.at(p).meta().pos_y( imgpos_row );
//         //std::cout << " imgpos[1] p=" << p << ", currentnode=" << current_node[p] << ": " << imgpos[1] << " tick=" << imgpos_row << std::endl;
//         ip++;
//       }
//       if ( ip==0 ) {
//         //std::cout << "no good plane?" << std::endl;
//         continue;
//       }
//       avetick /= float(ip);
//       //std::cout << "istep=" << istep << " avetick=" << avetick << " ip=" << ip << std::endl;
      
//       // now intersect depending on number of good wires
//       int crosses = 0;
//       std::vector<float> intersection;
//       if ( ngood==2 ) {
//         crosses = 0;
//         larcv::UBWireTool::wireIntersection( pid[0], wireid[0], pid[1], wireid[1], intersection, crosses );
//         if ( pid[0]==0 && pid[1]==1 ) pid[2] = 2;
//         else if ( pid[0]==0 && pid[1]==2 ) pid[2] = 1;
//         else if ( pid[0]==1 && pid[1]==2 ) pid[2] = 0;
//         else if ( pid[0]==1 && pid[1]==0 ) pid[2] = 2;
//         else if ( pid[0]==2 && pid[1]==0 ) pid[2] = 1;
//         else if ( pid[0]==2 && pid[1]==1 ) pid[2] = 0;
//         double worldpos[3] = { 0, (double)intersection[1], (double)intersection[0] };
//         wireid[2] = larutil::Geometry::GetME()->WireCoordinate( worldpos, pid[2] );
//         if ( wireid[2]<0 ) wireid[2]=0;
// 	wireid[2] = ( wireid[2]>=img_v.at(0).meta().max_x() ) ? img_v.at(0).meta().max_x()-1.0 : wireid[2];
//         imgcol[2] = wireid[2]/img_v.at(0).meta().pixel_width();
//         // we can backfill the missing bad plane
//         larcv::Pixel2D pix( imgcol[2], (int)img_v.at(pid[2]).meta().row( avetick ) );
//         track2d.at( badplane ).pixelpath.emplace_back( std::move(pix) );
//       }
//       else if ( ngood==3 ) {
//         crosses = 1;
// //      std::vector< std::vector<int> > wirelists(3);
// //      for (int p=0; p<3; p++)
// //        wirelists[p].push_back( wireid[p] );
// //      std::vector< std::vector<int> > intersections3plane;
// //      std::vector< std::vector<float> > vertex3plane;
// //      std::vector<float> areas3plane;
// //      std::vector< std::vector<int> > intersections2plane;
// //      std::vector< std::vector<float> > vertex2plane; 
// //      larcv::UBWireTool::findWireIntersections( wirelists, valid_range, intersections3plane, vertex3plane, areas3plane, intersections2plane, vertex2plane );
//         std::vector<int> wirelists(3,0);
//         for (int p=0; p<3; p++) {
//           wirelists[p] = wireid[p];
//           if ( wirelists[p]<0 ) wirelists[p] = 0;
//         }
//         std::vector<float> vertex3plane;
//         double tri_area = 0.0;
//         larcv::UBWireTool::wireIntersection( wirelists, intersection, tri_area, crosses );
//         //if ( tri_area>3.0 ) {
//         //std::cout << "warning, large intersection area for wires u=" << wireid[0]  << " v=" << wireid[1] << " y=" << wireid[2] << " area=" << tri_area << std::endl;
//         //}
//       }
      
//       std::vector<double> point3d(3,0.0);
//       point3d[0] = (avetick-trig_tick)*0.5*drift_v;
//       if ( intersection.size()>=2 ) {
//         point3d[1] = intersection[1];
//         point3d[2] = intersection[0];
//       }
      
//       // now update path variables
//       for (int p=0; p<3; p++) {
//         if ( !isgood[p] ) continue;
//         float stepsize = pathlength[p]/float(nsteps);
//         float next_ds = node_ds[p]+stepsize;
//         if ( current_node[p]<(int)edgelength[p].size() ) {
//           if ( next_ds>edgelength[p].at( current_node[p] ) ) {
//             node_ds[p] = next_ds-edgelength[p].at( current_node[p] );
//             current_node[p]++;
//           }
//           else {
//             node_ds[p] = next_ds;
//           }
//         }
//       }
      
//       if ( istep==0 ) {
//         // if first, log data
//         track3d.tick_start = avetick;
//         track3d.row_start  = img_v.at(0).meta().row( track3d.tick_start );
//         track3d.start_wire.resize(3,0);
//         track3d.start3D = point3d;
//         if ( ngood==3 ) {
//           for (int p=0; p<3; p++) {
//             track3d.start_wire[p] = wireid[p];
//           }
//         }
//         else if ( ngood==2 ) {
//           track3d.start_wire[pid[0]] = wireid[0];
//           track3d.start_wire[pid[1]] = wireid[1];
//           track3d.start_wire[pid[2]] = wireid[2];
//         }
//       }
//       else if ( istep+1==nsteps ) {
//         // set end info
//         track3d.tick_end = avetick;
//         track3d.row_end  = img_v.at(0).meta().row( track3d.tick_end );
//         track3d.end_wire.resize(3,0);
//         track3d.end3D = point3d;
//         if ( ngood==3 ) {
//           for (int p=0; p<3; p++) {
//             track3d.end_wire[p] = wireid[p];
//           }
//         }
//         else if ( ngood==2 ) {
//           track3d.end_wire[pid[0]] = wireid[0];
//           track3d.end_wire[pid[1]] = wireid[1];
//           track3d.end_wire[pid[2]] = wireid[2];
//         }
//       }
      
//       if ( track3d.path3d.size()==0 || track3d.path3d.back()!=point3d ) {
//         // std::cout << "step=" << istep << " node=(" << current_node[0] << "," << current_node[1] << "," << current_node[2] << ") "
//         //        << "tick=" << avetick << " "
//         //        << "pos=(" << point3d[0] << "," << point3d[1] << "," << point3d[2] << ") " 
//         //        << "(p=" << pid[0] << ", wire=" << wireid[0] << ") + "
//         //        << "(p=" << pid[1] << ", wire=" << wireid[1] << ") "
//         //        << "(p=" << pid[2] << ", wire=" << wireid[2] << ") "
//         //        << "crosses=" << crosses << " (isec=" << intersection.size() << ")" << std::endl;
//         track3d.path3d.emplace_back( std::move( point3d ) );
//       }
      
//     }//end of loop over steps

//     if ( ngood==2 ) {
//       larcv::Pixel2D pix( track2d.at(badplane).end.col, track2d.at(badplane).end.row );
//       //track2d.at(badplane).pixelpath.emplace_back( std::move(pix) );
//     }
    
//     // finishing touches:
//     // copy boundaryendpts
//     for (int p=0; p<3; p++) {
//       track3d.start_endpts.push_back( track2d.at(p).start );
//       track3d.end_endpts.push_back( track2d.at(p).end );
//     }
    
//     track3d.start_type = track3d.start_endpts.at(0).type;
//     track3d.end_type = track3d.end_endpts.at(0).type;
    

//     return track3d;
//   }
                                               
//   bool BoundaryMuonTaggerAlgo::compare2Dtrack( const std::vector< BMTrackCluster2D >& track2d, const BMTrackCluster3D& track3d, const larcv::ImageMeta& meta,
//                                                float path_radius_cm, float endpt_radius_cm ) {
//     // we get the start and end points of the 2D track
//     // we see if they are both close to the trajectory of the 3D track
//     // either point must be fruther away than the path_radius value
//     // end points must also be further away than the endpt_radius_cm value

//     float trig_tick = 2400;
//     float drift_v   = 0.106865; // cm/usec
//     float pix_rows_to_cm = meta.pixel_height()*0.5*drift_v; // cm/row
//     float pix_cols_to_cm = meta.pixel_width()*0.3; // cm

//     // std::cout << " 3D track start points (w,t): " << std::endl;
//     // for (int p=0; p<3; p++ )
//     //   std::cout << "  p" << p << ":: " << track3d.start_endpts[p].w << ", " << track3d.start_endpts[p].t << std::endl;
//     // std::cout << " 3D track end points (w,t): " << std::endl;
//     // for (int p=0; p<3; p++ )
//     //   std::cout << "  p" << p << ":: " << track3d.end_endpts[p].w << ", " << track3d.end_endpts[p].t << std::endl;
      

//     // check end points first
//     float start_dists[2] = { 0.};
//     float end_dists[2]   = { 0.};
//     float distances[3][2];
//     for ( int p=0; p<3; p++) {
//       if ( track2d.at(p).pixelpath.size()==0 ) {
//         distances[p][0] = -1;
//         distances[p][1] = -1;
//         continue;
//       }
//       //std::cout << "testing 2d endpoint p=" << p << "(" << track2d[p].start.w << ", " << track2d[p].start.t << ")" << std::endl;
//       float dt = (track2d[p].start.row-track3d.start_endpts[p].row)*pix_rows_to_cm;
//       float dw = (track2d[p].start.col-track3d.start_endpts[p].col)*pix_cols_to_cm;
//       start_dists[0] = sqrt( dt*dt + dw*dw );
//       // in case track reversed
//       dt = (track2d[p].end.row-track3d.start_endpts[p].row)*pix_rows_to_cm;
//       dw = (track2d[p].end.col-track3d.start_endpts[p].col)*pix_cols_to_cm;
//       start_dists[1] = sqrt( dt*dt + dw*dw );
      

//       dt = (track2d[p].end.row-track3d.end_endpts[p].row)*pix_rows_to_cm;
//       dw = (track2d[p].end.col-track3d.end_endpts[p].col)*pix_cols_to_cm;
//       end_dists[0] = sqrt( dt*dt + dw*dw );
//       dt = (track2d[p].start.row-track3d.end_endpts[p].row)*pix_rows_to_cm;
//       dw = (track2d[p].start.col-track3d.end_endpts[p].col)*pix_cols_to_cm;
//       end_dists[1] = sqrt( dt*dt + dw*dw );
      
//       int use_idx = 0;
//       if ( start_dists[0] > start_dists[1] ) use_idx = 1;
//       distances[p][0] = start_dists[use_idx];
//       distances[p][1] = end_dists[use_idx];
      
//       //std::cout << "distance between 2d endpoints on plane=" << p << " start=" << distances[p][0] << " end=" << distances[p][1] << std::endl;
//     }
    
//     bool allclose = true;
//     for (int p=0; p<3; p++) {
//       if ( distances[p][0]<0 ) continue;
//       if ( distances[p][0]>endpt_radius_cm || distances[p][1]>endpt_radius_cm )  {
//         allclose = false;
//         break;
//       }
//     }
//     if ( allclose ) {
//       //std::cout << "track too similar in end points" << std::endl;
//       return true;
//     }

    
//     // ok now check along 3D path
//     // need 3d position of end points
//     double start_area;
//     double end_area;
//     std::vector< float > start_poszy;
//     std::vector< float > end_poszy;
//     std::vector<int> start_wids;
//     std::vector<int> end_wids;
//     std::vector<int> goodplanes;
//     int crosses[2] = {0};
//     int ngoodplanes = 0;
//     for (int p=0; p<3; p++) {
//       if ( track2d[p].pixelpath.size()==0 ) {
//         continue;
//       }
//       ngoodplanes++;
//       goodplanes.push_back(p);
//       start_wids.push_back( (int)track2d[p].start.col*meta.pixel_width() );
//       end_wids.push_back(   (int)track2d[p].end.col*meta.pixel_width() );
//     }
    
//     if ( ngoodplanes==3 ) {
//       larcv::UBWireTool::wireIntersection( start_wids, start_poszy, start_area, crosses[0] );
//       larcv::UBWireTool::wireIntersection( end_wids, end_poszy, end_area, crosses[1] );
//     }
//     else {
//       larcv::UBWireTool::wireIntersection( goodplanes[0], start_wids[0], goodplanes[1], start_wids[1], start_poszy, crosses[0] );
//       larcv::UBWireTool::wireIntersection( goodplanes[0], end_wids[0],   goodplanes[1], end_wids[1],   end_poszy,   crosses[1] );
//     }

//     // if ( crosses[0]==0 ) {
//     //   std::cout << "Start point doesn't make a valid 3D point! start wid=(" << start_wids[0] << "," << start_wids[1];
//     //   if ( ngoodplanes==3 )
//     //  std::cout << "," << start_wids[2];
//     //   std::cout << ")" << std::endl;
//     // }
//     // if ( crosses[1]==0 ) {
//     //   std::cout << "End point doesn't make a valid 3D point! end wid=(" << end_wids[0] << "," << end_wids[1];
//     //   if ( ngoodplanes==3 )
//     //  std::cout << "," << end_wids[2];
//     //   std::cout << ")" << std::endl;
//     // }
//     if ( crosses[0]==0 || crosses[1]==0 )
//       return false;
    
//     std::vector< double > start3d(3,0.0);
//     std::vector< double > end3d(3,0.0);

//     for (int p=0; p<3; p++) {
//       start3d[0] += (track2d[p].start.row*meta.pixel_height()-trig_tick)/3.0;
//       end3d[0]   += (track2d[p].end.row*meta.pixel_height()-trig_tick)/3.0;
//     }
//     start3d[0] *= (0.5*drift_v);
//     start3d[1] = start_poszy[1];
//     start3d[2] = start_poszy[0];

//     end3d[0] *= (0.5*drift_v);
//     end3d[1] = end_poszy[1];
//     end3d[2] = end_poszy[0];

//     // std::cout << "track 2d 3D endpoints: "
//     //        << " start=(" << start3d[0] << "," << start3d[1] << "," << start3d[2] << ") "
//     //        << " end=(" << end3d[0] << "," << end3d[1] << "," << end3d[2] << ") "
//     //        << std::endl;
//     // std::cout << "track3d 3D endpoints: "
//     //        << " start=(" << track3d.path3d.front()[0] << "," << track3d.path3d.front()[1] << "," << track3d.path3d.front()[2] << ") "
//     //        << " end=(" << track3d.path3d.back()[0] << "," << track3d.path3d.back()[1] << "," << track3d.path3d.back()[2] << ") "
//     //        << std::endl;
    
//     // ok now that 3d start and end position found, check how close they are to the 3d path
//     // we use the geoalgo tools from larlite
//     // note that vector<double> has been typedefd as geoalgo::Point_t 
//     ::geoalgo::GeoAlgo algo;
//     ::geoalgo::Trajectory_t traj( track3d.path3d );
//     ::geoalgo::Point_t  start_pt( start3d );
//     ::geoalgo::Point_t    end_pt( end3d );
//     double closest_dist_start     = sqrt( algo.SqDist( start_pt, traj ) );
//     double closest_dist_end       = sqrt( algo.SqDist( end_pt, traj ) );
//     //std::cout << "closest distance of start/end points to trajectory: start=" << closest_dist_start << " end=" << closest_dist_end << std::endl;
//     if ( closest_dist_start<path_radius_cm && closest_dist_end<closest_dist_end ) {
//       return true;
//     }

//     return false;
//   }

  void BoundaryMuonTaggerAlgo::CollectCandidateBoundaryPixels( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchs,
    std::vector< dbscan::dbPoints >& combo_points, std::vector< std::vector< std::vector<int> > >& combo_cols,
    std::vector< larcv::Image2D >& matchedspacepts ) {
    // goal of this method is to search the images for pixels consistent with boundary crossings
    // (U,V,Y) combinations that meet at the edges of the detector are stored in class attributes matchalgo_loose, matchalgo_tight
    // for each row in the image, we send in information about which columns have charge above it
    // we then get in return whether a set of hits across each plane match to a boundary crossing

    const int nplanes = (int)imgs.size();
    if ( nplanes!=3 ) {
      throw std::runtime_error("CollectCandidateBoundaryPixels requires 3 planes. Sorry. Feel free to improve this.");
    }

    const int ncrossings = (int)BoundaryMuonTaggerAlgo::kNumCrossings;
    const larcv::ImageMeta& meta = imgs.at(0).meta();
    TRandom3 rand( time(NULL) );
    //TRandom3 rand( 12345 );        

    // we need the wire downsampling factor
    int dsfactor = int( meta.pixel_width()+0.1 ); // add 0.1 to round properly
    
    // we save col number of pixels above threshold in each plane with this vector.
    // it is the size of the neighborhood we'll look at. _config.neighborhoods defines the radius.
    std::vector< std::vector<int> > abovethresh(nplanes); 
    for (int p=0; p<nplanes; p++) {
      abovethresh[p].resize( 2*_config.neighborhoods.at(p)+1, 0 );
    }

    // reserve space for hit vector. marks the clumns where hits occur for a given a row.
    // we have one for each plane;
    std::vector< std::vector< int > > hits(nplanes);
    for (int p=0; p<nplanes; p++ ) {
      hits[p].resize( meta.cols() );
    }

    // storage for boundary combinations we'll cluster
    if ( (int)combo_points.size()!=ncrossings ) {
      combo_points.clear();
      combo_points.resize(ncrossings);
    }
    if ( (int)combo_cols.size()!=ncrossings ) {
      combo_cols.clear();
      combo_points.resize(ncrossings);
    }

    // misc. trackers
    int total_combos = 0;

    // now loop over over the time of the images
    for (size_t r=0; r<meta.rows(); r++) {

      // Find boundary combos using search algo
      //std::vector< int > hits[3];
      int nhits[nplanes];
      //std::cout << "start r="<< r << std::endl;
      for (int p=0; p<nplanes; p++) {
        // for given row, we mark which columns have pixels above threshold.
        // we mark that column but also -neighborhoods and +neighboorhood columns around the central pixel
        // the hit markers go into the hits vector
        nhits[p] = 0;
        memset( hits[p].data(), 0, sizeof(int)*hits[p].size() ); // clear hit marker with an old-fashion C-method
        const larcv::Image2D& img = imgs.at(p);
        //const larcv::Image2D& badchimg = badchs.at(p);
        for (int c=0; c<(int)meta.cols(); c++) {
          int wid = dsfactor*c;
          float val = img.pixel( r, c );
          //int lastcol = -1;
          if ( val > _config.thresholds.at(p) ) {
            for (int n=-_config.neighborhoods[p]; n<=_config.neighborhoods[p]; n++) {
              // bound check the requested column
              if ( wid+n>=(int)hits[p].size() ) continue;
              if ( wid+n<0 ) continue;
              hits[p][wid+n] = 1;
              //lastcol = wid+n;
            }
            nhits[p]++;
          }
          else if ( _config.hitsearch_uses_badchs && badchs.at(p).pixel( r, c )>0 ) {
            // use badchs. can toggle off/on as this increases false positive rate bigly
            hits[p][wid] = 1;
            //lastcol = c;
          }
          // because we mark columns that are +neighboorhood from current column, we can skip
          // to the last column-1 (end of loop will apply +1) to speed up slightly
          //if ( lastcol>=0 && lastcol<(int)meta.cols() && lastcol+1>c )
          //  c = lastcol-1;
        }
      }//end of plane loop
      // time this hit collecting
      //std::cout << "[row=" << r << ",t=" << meta.pos_y(r) << "] hits=" << nhits[0] << "," << nhits[1] << "," << nhits[2] << std::endl;
      
      // get boundary combos consistent which charge hits
      // two passes, a tight and loose pass
      std::vector< std::vector<BoundaryCombo> > matched_combos(ncrossings);
      matchalgo_tight->findCombos( hits[0], hits[1], hits[2], 
                                   badchs, true,
                                   matched_combos );
      matchalgo_loose->findCombos( hits[0], hits[1], hits[2], 
                                   badchs, true,
                                   matched_combos );
      
      // mark up image, filter out combinations for clustering
      for ( int pt=0; pt<(int)matched_combos.size(); pt++ ) {
        const std::vector<BoundaryCombo>& combos = matched_combos.at(pt);
        //std::cout << "  combos: type=" << pt << " number=" << combos.size() << std::endl;
        int idx_combo = -1;
        for ( auto &combo : combos ) {
          idx_combo++; // we index the vector combos, so others can refer to it
          int wirecols[3] = { combo.u()/dsfactor, combo.v()/dsfactor, combo.y()/dsfactor };
          // above boundary combo search will count bad channels as "hits"          
          // we require that at least 2 planes have good wire hits.
          int nbadchs = 0;
          for ( int p=0; p<nplanes; p++) {
            if ( badchs.at(p).pixel(r,wirecols[p])>0 ) 
              nbadchs++;
          }
          if (nbadchs>=2) {
            continue; // combo due two badchs
          }

          // otherwise set match

          // get position to mark up images
          // we plot top/bottom points in (z,t) (as y fixed by boundary)
          // we plot upstream/downstream in (y,t) (as z fixed by boundary)
          float x = 0;
          if ( pt==BoundaryMuonTaggerAlgo::top || pt==BoundaryMuonTaggerAlgo::bot ) {
            // top and bottom use the z value
            x = combo.pos[2];
          }
          else {
            // up stream and downstream use the y value
            x = combo.pos[1]+117.0; // y-values come back between [-117,117], we offset to have all positive values
          }

          // get average charge
          float charge = 0.0;
          int nhascharge = 0;
          for (int p=0; p<3; p++) {
            if ( imgs.at(p).pixel( r, wirecols[p] )>_config.thresholds.at(p) 
                 || badchs.at(p).pixel( r, wirecols[p] )>0 ) {
              nhascharge++;
              charge += imgs.at(p).pixel( r, wirecols[p] );
            }
          }
          charge /= float( 3.0-nbadchs );
          if ( charge>_config.thresholds.at(0) ) {
             if ( _config.save_endpt_images ) {
              // for optional visualization
              float prev_val = matchedspacepts.at( pt ).pixel( r, x ); // detector-space image
              matchedspacepts.at(pt).set_pixel( r, x, prev_val+charge );
              // we use matchedpixels image for later, when clusters are transferred back to image-space
              //for (int p=0; p<3; p++) {
              //  prev_val = matchedpixels.at(3*pt+p).pixel( r, wirecols[p] );
              //  matchedpixels.at(3*pt+p).set_pixel( r, wirecols[p], prev_val+charge );
              //}
            }
            // save (x,y) point for clustering
            std::vector<double> pt_combo(2,0.0); // this is (z,y)
            //pt_combo[0] = x;
            //pt_combo[1] = r + 0.1*float(idx_combo%10);
            pt_combo[0] = x + 1.0e-6*(1.0+rand.Uniform()); // prevents exact sample point from messing up spatial search tree
	    pt_combo[1] = r;
            combo_points[pt].emplace_back( std::move(pt_combo) );
            // save (u,v,y) column
            std::vector<int> combo_col(3);
            for (int p=0; p<3; p++) combo_col[p] = wirecols[p];
            combo_cols[pt].emplace_back( combo_col );
            total_combos++;
          }//if passes charge threshold
        }//end of loop over combos
      }//end of end point types
    }//end of row loop
  }//end of function Collect...


  void BoundaryMuonTaggerAlgo::ClusterBoundaryHitsIntoEndpointCandidates( const std::vector<larcv::Image2D>& imgs, 
									  const std::vector<larcv::Image2D>& badchs, 
									  const std::vector< dbscan::dbPoints >& combo_points, const std::vector< std::vector< std::vector<int> > >& combo_cols, 
									  std::vector< BoundarySpacePoint >& end_points,
									  std::vector< larcv::Image2D >& matchedpixels ) {
    // we cluster the combo_points (in detector space) and produce a list of BoundaryEndPt's (a set for each boundary type)
    // inputs
    // -------
    // imgs: ADC image
    // badchs: badch images
    // combo_points: realspace point list
    // combo_cols: image-coordinates. in sync with combo_points
    // 
    // outputs
    // -------
    // end_points: output spacepoints
    // matchedpixels: 

    const int nplanes = (int)imgs.size();
    const int ncrossings = BoundaryMuonTaggerAlgo::kNumCrossings;

    std::vector< larcv::Image2D > workspace;
    for (int p=0; p<nplanes; p++) {
      workspace.push_back( larcv::Image2D( imgs.at(p).meta() ) );
    }
    for (int pt=0; pt<ncrossings; pt++) {

      // We cluster the matched points
      dbscan::DBSCANAlgo dbalgo;
      //std::cout << "  starting dbscan for boundary type: " << pt << ". number of points: " << combo_points[pt].size() << std::endl;
      if ( (int)combo_points[pt].size()<_config.boundary_cluster_minpixels.at(0)) {
        //for (int i=0; i<combo_points[pt].size(); i++) {
        //std::cout << " (" << i << "): " << combo_points[pt].at(i)[0] << ", " << combo_points[pt].at(i)[1] << std::endl;
        //}
        //std::cout << "   not enough points to make a cluster: skipping clustering" << std::endl;
        continue;
      }
      dbscan::dbscanOutput clout = dbalgo.scan( combo_points[pt], _config.boundary_cluster_minpixels.at(0), _config.boundary_cluster_radius.at(0), false, 0.0 );
      //ann::ANNAlgo::cleanup();

      //std::cout << "  number of clusters: " << clout.clusters.size() << std::endl;
      // for each cluster we want to place an endpoint in image space
      // requirements
      //   1) place each point end point on a pixel with charge. Otherwise cause problems.
      //   2) each cluster will have an endpoint for each plane
      //   3) we should pick end points that sit on a cluster of charge in image-space that are similar
      //   4) 
      for (size_t ic=0; ic<clout.clusters.size(); ic++) {
        // loop over clusters in the real space points
        
        if ( (int)clout.clusters.at(ic).size() >= _config.boundary_cluster_minpixels.at(0) ) {
          //std::cout << "Find the endpoints for cluster pt=" << pt << " id=" << ic << " size=" << clout.clusters.at(ic).size() << std::endl;

          const dbscan::dbCluster& detspace_cluster = clout.clusters.at(ic);
          BoundarySpacePoint sppt = DefineEndpointFromBoundaryCluster( (BoundaryMuonTaggerAlgo::Crossings_t)pt, detspace_cluster, imgs, badchs, 
								       combo_points.at(pt), combo_cols.at(pt), matchedpixels );
          if (nplanes==(int)sppt.size()) {
            end_points.emplace_back( std::move(sppt) );
          }
          
        }//end of if cluster size is large enough
      }//end of cluster loop
    }//end of boundary point type

  }//end of clustering function

  BoundarySpacePoint BoundaryMuonTaggerAlgo::DefineEndpointFromBoundaryCluster( const BoundaryMuonTaggerAlgo::Crossings_t crossing_type, 
										const dbscan::dbCluster& detspace_cluster, 
										const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchs,
										const dbscan::dbPoints& combo_points, const std::vector< std::vector<int> >& combo_cols,
										std::vector<larcv::Image2D>& matchedpixels ) {

    // inputs
    // ------
    // crossing_type: boundary crossing type. will be given to output spacepoint
    // detspace_cluster: dbscan cluster of realspace hits. collection of indexes
    // imgs: ADC image
    // badch: badchimg
    // combo_points: all realspace hits
    // combo_cols: corresponding image-space hits in sync with combo_points
    //
    // output
    // ------
    // matchedpixels: 
    
    const int nplanes = (int)imgs.size();

    // we transfer information from this cluster into image space.  we mark pixels in image space with a hit
    std::vector< larcv::Image2D > workspace;
    for (int p=0; p<nplanes; p++) {
      workspace.push_back( larcv::Image2D( imgs.at(p).meta() ) );
      workspace.back().paint(0.0);
    }    

    // loop through hit in real space cluster. collect pixels in image space to cluster once more.
    TRandom3 rand( time(NULL) );
    //TRandom3 rand( 12345 );    

    // we organize image-space pixels in a vector that is sorted by (row,charge)
    // this is how we separate them
    struct PixelPt_t {
      std::vector<int> cols;
      int row;
      int nplanes_w_charge;
      int idxhit;
      float q;
      float boundary_dim;
    };
    struct PixelSorter_t {
      // we sort by number:
      // 1) planes with charge
      // 2) charge
      // 3) closest to boundary (can be stocastic)
      bool operator()(PixelPt_t lhs, PixelPt_t rhs) {
	if ( lhs.nplanes_w_charge>rhs.nplanes_w_charge ) {
	  return true;
	}
	else if ( lhs.nplanes_w_charge==rhs.nplanes_w_charge ) {
	  // next comparison: prox to edge
	  if ( lhs.boundary_dim < rhs.boundary_dim ) {
	    return true;
	  }
	  else if ( lhs.boundary_dim==rhs.boundary_dim ) {
	    // next comparison: q
	    if ( lhs.q > rhs.q ) {
	      return true;
	    }
	  }
	}
        return false;
      };
    } my_sorter;
    std::vector< PixelPt_t > sorted_imagespace_pixels(detspace_cluster.size());
    
    // loop overcluster. collect pixel information
    for (size_t ihit=0; ihit<detspace_cluster.size(); ihit++) {
      PixelPt_t& pixel = sorted_imagespace_pixels[ihit];
      int idxhit = detspace_cluster.at(ihit);
      pixel.idxhit = idxhit;
      pixel.cols.resize(3,0);      
      int row = combo_points[idxhit][1];
      pixel.row = row;
      pixel.q = 0;
      pixel.nplanes_w_charge = 0;

      if ( crossing_type==0 ) {
	// top: get close to 117
	pixel.boundary_dim = fabs(117-combo_points[idxhit][0]);
      }
      else if ( crossing_type==1 ) {
	// bottom: get close to 0
	pixel.boundary_dim = fabs(combo_points[idxhit][0]);
      }
      else if ( crossing_type==2 ) {
	// upstream: get close to 0
	pixel.boundary_dim = fabs(combo_points[idxhit][0]);
      }
      else if ( crossing_type==3 ) {
	// downstream: get close to 1036
	pixel.boundary_dim = fabs(combo_points[idxhit][0]-1036.0);
      }
      else {
	// never get here
	assert(false);
      }
      
      for (int p=0; p<nplanes; p++) {
        int col = combo_cols[idxhit][p]; // get the position in the image for this point
	pixel.cols[p] = col;
        // look for charge in neighborhood of this point
        int neighborhood = 0; // should not use a neighborhood
        for (int n=-neighborhood; n<=neighborhood; n++) {
          if ( col+n<0 || col+n>=(int)imgs.at(p).meta().cols() ) continue;
          float badchq = badchs.at(p).pixel( (int)combo_points[idxhit][1], col+n );
          float q      = imgs.at(p).pixel( (int)combo_points[idxhit][1], col+n );          
          if ( ( q > _config.thresholds.at(p) || (badchq>0 && n==0) ) && workspace[p].pixel(row,col+n)==0 ) {

            workspace[p].set_pixel( row, col+n, 10.0 );
            if ( _config.save_endpt_images ) {
              matchedpixels.at(nplanes*crossing_type+p).set_pixel( (int)combo_points[idxhit][1], col+n, 255.0 );
            }

	    pixel.nplanes_w_charge++;
	    pixel.q += q;
	    
          }//if pixel above thresh
        }// loop over neighborhood
      }//end of loop over plane
    }//end of loop over hits in real space
    
    // we sort
    std::sort( sorted_imagespace_pixels.begin(), sorted_imagespace_pixels.end(), my_sorter );

    // // debug: sort check
    // std::cout << "pixel sort check: " << std::endl;
    // for ( auto& pixel : sorted_imagespace_pixels ) {
    //   std::cout << " nplanes=" << pixel.nplanes_w_charge << " x=" << pixel.boundary_dim << " q=" << pixel.q << std::endl;
    // }
    // std::cout << std::endl;
    // std::cin.get();

    
    // we pick the best pixel to make the space point
    PixelPt_t& best = sorted_imagespace_pixels.front();

    if ( best.nplanes_w_charge<=1 ) {
      // bad best point
      BoundarySpacePoint emptysp;
      return emptysp; // return empty container
    }

    std::vector< BoundaryEndPt > endpt_v;
    for ( int p=0; p<nplanes; p++)  {
      BoundaryEndPt endpt(best.row,best.cols[p],CrossingToBoundaryEnd(crossing_type) );
      endpt_v.emplace_back( std::move(endpt) );
    }
    float x = (imgs.front().meta().pos_y( best.row ) - 3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);
    BoundarySpacePoint spacepoint( CrossingToBoundaryEnd(crossing_type), std::move(endpt_v), imgs.front().meta() );
    
    return spacepoint;
  }//end of end point definition

  void BoundaryMuonTaggerAlgo::GenerateEndPointMetaData( const std::vector< std::vector< BoundaryEndPt > >& endpts, const std::vector<larcv::Image2D>& img_v,
    const int rmax_window, const int rmin_window, const int col_window, 
    std::vector< std::vector<dbscan::ClusterExtrema> >& candidate_metadata ) {

    const int nplanes = (int)img_v.size();

    for ( auto const& endpt_v : endpts ) {
      // we want to grab some meta data about the end point
      // we want to cluster the pixels inside a window around the end point

      int row = endpt_v.front().row; // clusters on each plane should have same row
      // define the window
      int row_start = row - rmin_window;
      int row_end   = row + rmax_window;

      dbscan::DBSCANAlgo algo;

      // goals on each plane:
      // 1) collect hits in neighborhood
      // 2) see if cluster connects end point and window boundaries
      // 2) somehow need to track if end point is surrounded by badchs      

      // we collect the following for each plane (later can be part of struct)
      std::vector<dbscan::dbscanOutput> clout_v;
      std::vector<dbscan::dbPoints> hits_v;
      std::vector<dbscan::ClusterExtrema> extrema_v;

      for (int p=0; p<nplanes;p++) {
        dbscan::dbPoints planepts;
        int col = endpt_v.at(p).col;
        int col_start = col-col_window;
        int col_end   = col+col_window;
        for ( int r=row_start; r<=row_end; r++ ) {
          if ( r<0 || r>=(int)img_v.at(p).meta().rows() ) continue;
          for (int c=col_start; c<=col_end; c++ ) {
            if ( c<0 || c>=(int)img_v.at(p).meta().cols() ) continue;

            if ( img_v.at(p).pixel(r,c)>_config.thresholds.at(p) ) {
              std::vector<double> hit(2);
              hit[0] = c;
              hit[1] = r;
              planepts.emplace_back( std::move(hit) );
            }

          }
        }
        dbscan::dbscanOutput clout = algo.scan( planepts, 3, 5.0, false, 0.0 );
        std::vector<double> centerpt(2);
        centerpt[0] = col;
        centerpt[1] = row;
        int matching_cluster = clout.findMatchingCluster( centerpt, planepts, 3.0 );

        // initialize extrema search
        std::vector<dbscan::ClusterExtrema> cl_extrema;
        for (int ic=0; ic<(int)clout.clusters.size();ic++) {
          dbscan::ClusterExtrema extrema = dbscan::ClusterExtrema::FindClusterExtrema( clout.clusters.at(ic), planepts );
          cl_extrema.emplace_back( std::move(extrema) );
        }

        // fill the output
        if ( matching_cluster>=0 ) {
          extrema_v.emplace_back( std::move(cl_extrema.at(matching_cluster)) );
        }
        else {
          // create empty extrema object
          dbscan::ClusterExtrema empty = dbscan::ClusterExtrema::MakeEmptyExtrema();
          extrema_v.emplace_back( std::move(empty) );
        }

        clout_v.emplace_back( std::move(clout) );
        hits_v.emplace_back( std::move(planepts) );
      }//end of plane loop to gather information

      candidate_metadata.emplace_back( std::move(extrema_v) );
    }//end of loop over end points
  }

  void BoundaryMuonTaggerAlgo::CheckClusterExtrema(  const std::vector< BoundarySpacePoint >& endpts, 
    const std::vector< std::vector<dbscan::ClusterExtrema> >& candidate_metadata, const std::vector<larcv::Image2D>& img_v,
    std::vector<int>& passes_check ) {
    // we we see if we can match cluster extrema for 2 of 3 planes.
    // this is to ID clusters that have been put on completely different tracks

    if ( passes_check.size()!=endpts.size() )
      passes_check.resize(endpts.size(),1);

    const larcv::ImageMeta& meta = img_v.at(0).meta();

    for ( int idx_endpt=0; idx_endpt<(int)endpts.size(); idx_endpt++ ) {
      const BoundarySpacePoint& endpt_v = endpts.at(idx_endpt);
      const std::vector< dbscan::ClusterExtrema >& endpt_metadata_v = candidate_metadata.at(idx_endpt);
      std::vector<int> wids_top( endpts.size() );
      std::vector<int> wids_bot( endpts.size() );
      int nfilled = 0;
      for (int p=0; p<(int)endpt_v.size(); p++) {
        if ( !endpt_metadata_v.at(p).isempty() ) {
          wids_top[p] = img_v.at(p).meta().pos_x( endpt_metadata_v.at(p).topmost()[0] );
          wids_bot[p] = img_v.at(p).meta().pos_x( endpt_metadata_v.at(p).bottommost()[0] );        
          nfilled++;
        }
      }
      if ( nfilled==3 ) {
        int crosses_t = 0;
        std::vector<float> poszy_t(2);
        double tri_t = 0;
        larcv::UBWireTool::wireIntersection( wids_top, poszy_t, tri_t, crosses_t );

        int crosses_b = 0;
        std::vector<float> poszy_b(2);
        double tri_b = 0;
        larcv::UBWireTool::wireIntersection( wids_bot, poszy_b, tri_b, crosses_b );

        std::cout << "EndPt #" << idx_endpt << " tick=" << meta.pos_y( endpt_v.at(0).row )  << " cols=(" << endpt_v.at(0).col << "," << endpt_v.at(1).col << "," << endpt_v.at(2).col << ")"
          << " top-intersection area=" << tri_t << " bot-intersection area=" << tri_b << std::endl;
      }
      else {
        std::cout << "EndPt #" << idx_endpt << " tick=" << meta.pos_y( endpt_v.at(0).row )  << " cols=(" << endpt_v.at(0).col << "," << endpt_v.at(1).col << "," << endpt_v.at(2).col << ")"
          << " only has charge on " << nfilled << " planes" << std::endl;
      }
    }
  }

  void BoundaryMuonTaggerAlgo::SelectOnlyTrackEnds( const std::vector< BoundarySpacePoint >& endpts, 
    const std::vector<larcv::Image2D>& img_v, const int rmax_window, const int rmin_window, const int col_width,
    std::vector< int >& cluster_passed ) {
    // we want to remove false positive flash-end detections.
    // this means removing flash-tags that occur in the middle of a track.
    // such tags have the consequence of causing many multiple tracks, greatly inflatingn the number of aStar searchs
    // (though in principle, these repeate 3D tracks can be removed/merged -- but let's try to get rid of them here)

    // we store pass/fail tag here
    cluster_passed.resize( endpts.size(), 0 );

    int idx_cluster = 0;
    for ( auto const& endpt_v : endpts ) {
      // to determine if track end, we do something very simple:
      // we define a window around the cluster pt in all planes
      // within the window we collect charge and cluster.
      // we find the extrema of all the points
      // the extrema tells us whether to check the horizontally of vertically
      // does both or only one extrema reach the end?
      // if only one, it passes
      // if both, it fails as being midtrack

      const int nplanes = (int)img_v.size();
      int row = endpt_v.front().row; // clusters on each plane should have same row
      // define the window
      int row_start = row - rmin_window;
      int row_end   = row + rmax_window;

      dbscan::DBSCANAlgo algo;

      // goals on each plane:
      // 1) collect hits in neighborhood
      // 2) see if cluster connects end point and window boundaries
      // 2) somehow need to track if end point is surrounded by badchs      

      std::vector<dbscan::dbscanOutput> clout_v;
      std::vector<dbscan::dbPoints> hits_v;
      std::vector<int> matching_cluster_v;
      std::vector<std::vector<dbscan::ClusterExtrema> > extrema_v;

      for (int p=0; p<nplanes;p++) {
        dbscan::dbPoints planepts;
        int col = endpt_v.at(p).col;
        int col_start = col-col_width;
        int col_end   = col+col_width;
        for ( int r=row_start; r<=row_end; r++ ) {
          if ( r<0 || r>=(int)img_v.at(p).meta().rows() ) continue;
          for (int c=col_start; c<=col_end; c++ ) {
            if ( c<0 || c>=(int)img_v.at(p).meta().cols() ) continue;

            if ( img_v.at(p).pixel(r,c)>_config.thresholds.at(p) ) {
              std::vector<double> hit(2);
              hit[0] = c;
              hit[1] = r;
              planepts.emplace_back( std::move(hit) );
            }

          }
        }
        dbscan::dbscanOutput clout = algo.scan( planepts, 3, 5.0, false, 0.0 );
        std::vector<double> centerpt(2);
        centerpt[0] = col;
        centerpt[1] = row;
        int matching_cluster = clout.findMatchingCluster( centerpt, planepts, 3.0 );

        // initialize extrema search
        std::vector<dbscan::ClusterExtrema> cl_extrema;
        for (int ic=0; ic<(int)clout.clusters.size();ic++) {
          dbscan::ClusterExtrema extrema = dbscan::ClusterExtrema::FindClusterExtrema( clout.clusters.at(ic), planepts );
          cl_extrema.emplace_back( std::move(extrema) );
        }
        clout_v.emplace_back( std::move(clout) );
        hits_v.emplace_back( std::move(planepts) );
        matching_cluster_v.push_back(matching_cluster);
        extrema_v.emplace_back( std::move(cl_extrema) );
      }//end of plane loop to gather information

      // now that we have the info we want
      // we try to see if we can connect the central cluster to the ends of the window
      std::vector<int> boundaries_reached(nplanes,0);
      std::vector<int> reached_by_badchs(nplanes,0); 
      int num_no_cluster = 0;

      for (int p=0; p<nplanes; p++) {
        int num_boundaries_reached = 0;
        int col = endpt_v.at(p).col;
        int col_start = col-col_width;
        int col_end   = col+col_width;
        if ( col_start<0 ) col_start = 0;
        if ( col_end>=(int)img_v.at(p).meta().cols()) col_end = (int)img_v.at(p).meta().cols()-1;
        if ( matching_cluster_v.at(p)<0) {
          std::cout << "  cluster #" << idx_cluster << " no cluster on plane=" << p << std::endl;
          num_no_cluster++;
          continue;
        }

        // first check the simplest thing. Does the central cluster reach the ends of the window?
        dbscan::ClusterExtrema& extrema = extrema_v.at( p ).at( matching_cluster_v.at(p) );
        if ( extrema.leftmost()[0]<=col_start ) num_boundaries_reached++;
        if ( extrema.rightmost()[0]>=col_end   ) num_boundaries_reached++;
        if ( extrema.bottommost()[1]<=row_start ) num_boundaries_reached++;
        if ( extrema.topmost()[1]>=row_end   ) num_boundaries_reached++;

        boundaries_reached[p] = num_boundaries_reached;

        std::cout << "  cluster #" << idx_cluster << " extrema: "
          <<  " top=(" << img_v.at(p).meta().pos_y(extrema.topmost()[1]) << ") vs row_start=" << img_v.at(p).meta().pos_y(row_end)
          <<  " bot=(" << img_v.at(p).meta().pos_y(extrema.bottommost()[1]) << ") vs row_end=" << img_v.at(p).meta().pos_y(row_start)
          <<  " left=("  << img_v.at(p).meta().pos_x(extrema.leftmost()[0]) << ") vs col_start=" << img_v.at(p).meta().pos_x(col_start)
          <<  " right=(" << img_v.at(p).meta().pos_x(extrema.rightmost()[0]) << ") vs col_end=" << img_v.at(p).meta().pos_x(col_end)
          <<  " boundaries reached=" << num_boundaries_reached
          << std::endl;                  

        // otherwise, we try to follow cluster to the ends of the boundary
        // bool end_reached = false;
        // implement this later

      }//end of plane loop

      int nplanes_with_crossing_cluster = 0;
      for (int p=0; p<nplanes; p++) {
        if ( boundaries_reached[p]>=2 ) nplanes_with_crossing_cluster++;
      }

      if ( nplanes_with_crossing_cluster<=1 && num_no_cluster<2 ) {
        // only keep these
        cluster_passed[idx_cluster] = 1;
      }
      std::cout << "filter cluster#" << idx_cluster << " tick=" << img_v.at(0).meta().pos_y(row) 
        << " cols=(" << endpt_v[0].col << "," << endpt_v[1].col << "," << endpt_v[2].col << ")"
        << " planes with no cluster=" << num_no_cluster
        << " planes with crossing cluster=" << nplanes_with_crossing_cluster
        << std::endl;
      idx_cluster++;
    }//end of cluster loop

  }// end of function




}//end of namespace
