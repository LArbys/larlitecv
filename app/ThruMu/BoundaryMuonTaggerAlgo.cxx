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

#include "AStar3DAlgo.h"
#include "Linear3DChargeTagger.h"
#include "Linear3DPostProcessor.h"
#include "AStarNodes2BMTrackCluster3D.h"

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

    TRandom3 rand( time(NULL) );

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

    // create Image2D objects to store the result of boundary matching
    // we use meta to make the same sized image as the input
    if ( _config.save_endpt_images ) {
      for (int i=0; i<ncrossings; i++) { // top, bottom, upstream, downstream x 3 planes
        larcv::Image2D matchimage( meta );
        matchimage.paint(0.0);
        matchedspacepts.emplace_back( std::move(matchimage) );
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
    CollectCandidateBoundaryPixels( imgs, badchs, combo_points, combo_cols, matchedspacepts );
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

    // we need the wire downsampling factor
    int dsfactor = int( meta.pixel_width()+0.1 ); 
    
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
            pt_combo[0] = x;
            //pt_combo[1] = r + 0.1*float(idx_combo%10);
            pt_combo[1] = r + 0.1*rand.Uniform(); // prevents exact sample point from messing up spatial search tree
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


  void BoundaryMuonTaggerAlgo::ClusterBoundaryHitsIntoEndpointCandidates( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchs, 
    const std::vector< dbscan::dbPoints >& combo_points, const std::vector< std::vector< std::vector<int> > >& combo_cols, 
    std::vector< BoundarySpacePoint >& end_points,
    std::vector< larcv::Image2D >& matchedpixels ) {
    // we cluster the combo_points (in detector space) and produce a list of BoundaryEndPt's (a set for each boundary type)

    const int nplanes = (int)imgs.size();
    const int ncrossings = BoundaryMuonTaggerAlgo::kNumCrossings;

    std::vector< larcv::Image2D > workspace;
    for (int p=0; p<nplanes; p++) {
      workspace.push_back( larcv::Image2D( imgs.at(p).meta() ) );
    }
    for (int pt=0; pt<ncrossings; pt++) {
      
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
    const dbscan::dbPoints& combo_points, const std::vector< std::vector<int> >& combo_cols, std::vector<larcv::Image2D>& matchedpixels ) {

    const int nplanes = (int)imgs.size();

    // we transfer information from this cluster into image space.  we mark pixels in image space with a hit
    std::vector< larcv::Image2D > workspace;
    for (int p=0; p<nplanes; p++) {
      workspace.push_back( larcv::Image2D( imgs.at(p).meta() ) );
      workspace.back().paint(0.0);
    }    

    // loop through hit in real space cluster. collect pixels in image space to cluster once more.
    TRandom3 rand( time(NULL) );    
    dbscan::dbPoints chargepts[3]; // charge per plane    
    std::vector<int> min_row(nplanes,-1);
    std::vector<int> max_row(nplanes,-1);
    std::vector<float> cluster_q(nplanes,0);
    int abs_min_row = -1;
    int abs_max_row = -1;

    // we organize image-space pixels in a vector that is sorted by (row,charge)
    struct PixelPt_t {
       int col;
       int row;
       float q;
    };
    struct PixelSorter_t {
      bool operator()(PixelPt_t lhs, PixelPt_t rhs) {
        if ( lhs.row<rhs.row ) return true;
        else if ( lhs.row==rhs.row ) {
          if ( lhs.q<rhs.q) return true;
        }
        return false;
      };
    } my_sorter;
    std::vector< std::vector<PixelPt_t> > sorted_imagespace_pixels(nplanes);

    for (size_t ihit=0; ihit<detspace_cluster.size(); ihit++) {
      int idxhit = detspace_cluster.at(ihit);
      for (int p=0; p<nplanes; p++) {
        int row = combo_points[idxhit][1];
        if ( abs_min_row<0 || row<abs_min_row )
          abs_min_row = row;
        if ( abs_max_row<0 || row>abs_max_row )
          abs_max_row = row;
        int col = combo_cols[idxhit][p]; // get the position in the image for this point
        // look for charge in neighborhood of this point
        int neighborhood = _config.neighborhoods[p]*_config.type_modifier[(int)crossing_type];
        for (int n=-neighborhood; n<=neighborhood; n++) {
          if ( col+n<0 || col+n>=(int)imgs.at(p).meta().cols() ) continue;
          float badchq = badchs.at(p).pixel( (int)combo_points[idxhit][1], col+n );
          float q      = imgs.at(p).pixel( (int)combo_points[idxhit][1], col+n );          
          if ( ( q > _config.thresholds.at(p) || (badchq>0 && n==0) ) && workspace[p].pixel(row,col+n)==0 ) {
            // define charge point in image space. jiggle col,row because ANN fails if points on top one another
            std::vector<double> imagespacept(2);
            imagespacept[0] = double(col+n) + 0.1*rand.Uniform();
            imagespacept[1] = double(row) + 0.1*rand.Uniform();

            chargepts[p].push_back( imagespacept );
            // mark this pixel as used
            workspace[p].set_pixel( row, col+n, 10.0 );
            if ( _config.save_endpt_images ) {
              matchedpixels.at(nplanes*crossing_type+p).set_pixel( (int)combo_points[idxhit][1], col+n, 255.0 );
            }

            // set max or min row
            if ( min_row[p]==-1 || min_row[p]>row)
              min_row[p] = row;
            if ( max_row[p]==-1 || max_row[p]<row )
              max_row[p] = row;

            cluster_q[p] += q;

            // fill PixelPt_t container
            PixelPt_t pix;
            pix.col = col+n;
            pix.row = row;
            pix.q = q;
            sorted_imagespace_pixels.at(p).emplace_back( std::move(pix) );

          }//if pixel above thresh
        }// loop over neighborhood
      }//end of loop over plane
    }//end of loop over hits in real space

    // we sort
    for ( int p=0; p<nplanes; p++ ) {
      std::sort( sorted_imagespace_pixels.at(p).begin(), sorted_imagespace_pixels.at(p).end(), my_sorter );
    }

    // with sorted list of image-space points, we now step through time, pairing same-time triples, to find consistent positions
    // we keep the pixel that is closest to the wall
    int current_row = abs_min_row;
    std::vector<int> plane_idx(nplanes,0);
    int best_row = 0;
    std::vector<int> best_cols(3,0);
    float best_dwall = -1;
    std::vector<float> best_poszy(2,0.0);
    bool more_combos = true;

    //std::cout << "Start combo search" << std::endl;

    while ( more_combos ) {
      // go to pixel in each plane that is at least this row
      std::vector<int> wids(nplanes);
      std::vector<int> plane_row(nplanes,0);
      for (int p=0; p<nplanes; p++) {
        // keep moving up pixel list until we find a pixel that is same row or further
        //std::cout << "in combo with plane_idx=" << plane_idx[p] << " number of image-space pixels=" << sorted_imagespace_pixels.at(p).size() << std::endl;

        while ( plane_idx[p]+1<(int)sorted_imagespace_pixels.at(p).size() && current_row > sorted_imagespace_pixels[p][plane_idx[p]].row  ) {
          plane_idx[p]++;
        }
        if ( plane_idx[p]<(int)sorted_imagespace_pixels.at(p).size() ) {
          wids[p]      = sorted_imagespace_pixels[p][plane_idx[p]].col;
          plane_row[p] = sorted_imagespace_pixels[p][plane_idx[p]].row;
        }
      }
      if ( wids[0]==0 && wids[1]==0 && wids[2]==0 )
        break; // empty sorted pixels. this is a garbage cluster

      // std::cout << "current_row=" << current_row 
      //  << " plane_idx's=[" 
      //  << " " << plane_idx[0] << "/" << sorted_imagespace_pixels[0].size() << "r=" << sorted_imagespace_pixels[0][plane_idx[0]].row << ","
      //  << " " << plane_idx[1] << "/" << sorted_imagespace_pixels[1].size() << "r=" << sorted_imagespace_pixels[1][plane_idx[1]].row << ","
      //  << " " << plane_idx[2] << "/" << sorted_imagespace_pixels[2].size() << "r=" << sorted_imagespace_pixels[2][plane_idx[2]].row << "]" 
      //  << std::endl;                         

      // are they the same row?
      int smallest_row = -1;
      for ( int p=0; p<nplanes; p++) {
        if ( smallest_row<0 || plane_row[p]<smallest_row ) {
          smallest_row = plane_row[p];
        }
      }
      if (smallest_row>current_row) {
        current_row = smallest_row;
        //std::cout << " -- out of sync. move to next." << std::endl;
        continue; // look again for an alignment
      }

      // otherwise we have same-tick pixel to test intersection
      std::vector<float> poszy;
      double triarea = 0.;
      int crosses = 0;
      larcv::UBWireTool::wireIntersection( wids, poszy, triarea, crosses );

      if ( crosses==1 && triarea<1 ) {
        // good intersection. find dwall
        float thisdwall = 0.;
        if ( crossing_type==top ) thisdwall = 117.0-poszy[1];
        else if ( crossing_type==bot ) thisdwall = poszy[1] + 117.0;
        else if ( crossing_type==upstream ) thisdwall = poszy[0];
        else if ( crossing_type==downstream ) thisdwall = 1040-poszy[0];
        //std::cout << "  -- valid point dwall= " <<    thisdwall << std::endl;
        if ( best_dwall<0 || thisdwall<best_dwall ) {
          best_dwall = thisdwall;
          best_row = current_row;
          best_cols = wids;
          best_poszy = poszy;
          //std::cout << "  -- update." << std::endl;
        }
      }

      // we increment to move forward. but which plane's list?
      // the one with the most remaining pixels at this current row i guess
      std::vector<int> num_remaining(nplanes,0);
      for ( int p=0; p<nplanes; p++) {
        for (int idx=plane_idx[p]+1; idx<(int)sorted_imagespace_pixels[p].size(); idx++) {
          if ( sorted_imagespace_pixels[p][idx].row==current_row ) num_remaining[p]++;
          else if ( sorted_imagespace_pixels[p][idx].row>current_row ) break;
        }
      }

      int increment_p = 0;
      int max_remaining = num_remaining[0];
      for (int p=1; p<nplanes; p++) {
        if ( num_remaining[p]>max_remaining) {
          increment_p = p;
          max_remaining = num_remaining[p];
        }
      }

      plane_idx[increment_p]++;
      for (int p=0; p<nplanes; p++){
        if ( plane_idx[p]<(int)sorted_imagespace_pixels[p].size()  ) {
          if ( sorted_imagespace_pixels[p][plane_idx[p]].row>current_row )
            current_row = sorted_imagespace_pixels[p][plane_idx[p]].row;
        }
        else {
          more_combos = false;
        }
      }
    }
    
    // we have defined an image-space cluster of pixels on each plane in chargepts 
    // we have the min and max row and charge for this cluster on all three planes
    // std::cout << "Candidate Endpt(type=" << crossing_type << ")"
    //         << " min ticks=(" << meta.pos_y(max_row[0]) << "," << meta.pos_y(max_row[1]) << "," << meta.pos_y(max_row[2]) << ") "
    //         << " max_rows=(" <<  meta.pos_y(min_row[0]) << "," << meta.pos_y(min_row[1]) << "," << meta.pos_y(min_row[2]) << ") "
    //         << " endpt-dwall=" << best_dwall
    //         << " endpt-tick=" << meta.pos_y(best_row)
    //         << " endpt-cols=" << best_cols[0] << "," << best_cols[1] << "," << best_cols[2] << ") "
    //         << " best-tri=" << best_triarea
    //         << std::endl;


    // create the end point
    std::vector< BoundaryEndPt > endpt_v; 

    if ( best_row==0 && best_cols[0]==0 && best_cols[1]==0 && best_cols[2]==0 ) {
      BoundarySpacePoint emptysp;
      return emptysp; // return empty container
    }

    for ( int p=0; p<nplanes; p++)  {
      BoundaryEndPt endpt(best_row,best_cols[p], CrossingToBoundaryEnd(crossing_type) );
      endpt_v.emplace_back( std::move(endpt) );
    }
    
    float x = (imgs.front().meta().pos_y( best_row ) - 3200.0)*(larutil::LArProperties::GetME()->DriftVelocity()*0.5);    
    BoundarySpacePoint spacepoint( CrossingToBoundaryEnd(crossing_type), std::move(endpt_v), x, best_poszy[1], best_poszy[0] );

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
