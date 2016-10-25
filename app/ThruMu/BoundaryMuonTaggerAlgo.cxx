#include "BoundaryMuonTaggerAlgo.h"

// std
#include <vector>
#include <cmath>
#include <assert.h>
#include <ctime>

// larcv
#include "UBWireTool/UBWireTool.h"

#include "AStarGridAlgo.h"
#include "BoundaryEndPt.h"
#include "BoundaryIntersectionAlgo.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

namespace larlitecv {

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
							 std::vector< std::vector<BoundaryEndPt> >& end_points,
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

    int ncrossings = 4;
    //int nplanes = imgs.size();

    if ( !_config.checkOK() )  {
      std::cout << "[BOUNDARY MUON TAGGER ALGO ERROR] Not configured." << std::endl;
      return kErr_NotConfigured;
    }
    
    if ( imgs.size()<3 ) {
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
    enum { top=0, bot, upstream, downstream };

    for (int i=0; i<ncrossings; i++) { // top, bottom, upstream, downstream x 3 planes
      larcv::Image2D matchimage( meta );
      matchimage.paint(0.0);
      matchedspacepts.emplace_back( std::move(matchimage) );
      for (int p=0; p<3; p++) {
	larcv::Image2D matchimage2( imgs.at(p).meta() );
	matchimage2.paint(0.0);
	matchedpixels.emplace_back( std::move(matchimage2) );
      }
    }

    // we need the wire downsampling factor
    int dsfactor = int( meta.pixel_width()+0.1 ); 
    const clock_t begin_time = clock();
    
    // we save col number of pixels above threshold in each plane with this
    std::vector<int> abovethresh[3]; 
    for (int p=0; p<3; p++) {
      abovethresh[p].resize( 2*_config.neighborhoods.at(p)+1, 0 );
    }

    // reserve space for hit vector. marks where hits occur.
    std::vector< int > hits[3];
    for (int p=0; p<3; p++ ) {
      hits[p].resize( meta.cols() );
    }

    // storage for boundary combinations we'll cluster
    dbscan::dbPoints combo_points[4];
    std::vector< std::vector<int> > combo_cols[4];

    // misc. trackers
    float tot_hit_collecting = 0;
    int total_combos = 0;
    
    // now loop over over the time of the images
    for (int r=0; r<meta.rows(); r++) {

      // Find boundary combos using search algo
      //std::vector< int > hits[3];
      const clock_t begin_hit_collecting = clock();
      int nhits[3] = {0};
      for (int p=0; p<3; p++) {
	memset( hits[p].data(), 0, sizeof(int)*hits[p].size() );
	const larcv::Image2D& img = imgs.at(p);                                                                                                                                
	//const larcv::Image2D& badchimg = badchs.at(p);
	for (int c=0; c<meta.cols(); c++) {
	  int wid = dsfactor*c;
	  float val = img.pixel( r, c );
	  int lastcol = -1;
	  if ( val > _config.thresholds.at(p) ) {
	    for (int n=-_config.neighborhoods[p]; n<=_config.neighborhoods[p]; n++) {
	      if ( wid+n>=hits[p].size() ) continue;
	      if ( wid+n<0 ) continue;
	      hits[p][wid+n] = 1;
	      lastcol = wid+n;
	    }
	    nhits[p]++;
	  }
// 	  if ( badchimg.pixel( r, c )>0 ) {
// 	    hits[p][wid] = 1;
// 	  }
	  if ( lastcol>=0 && lastcol<meta.cols() )
	    c = lastcol;
	}
      }
      tot_hit_collecting += float( clock()-begin_hit_collecting )/CLOCKS_PER_SEC;
      //std::cout << "[row=" << r << ",t=" << meta.pos_y(r) << "] hits=" << nhits[0] << "," << nhits[1] << "," << nhits[2] << std::endl;
      
      // get boundary combos consistent which charge hits
      std::vector< std::vector<BoundaryCombo> > matched_combos;
      matchalgo.findCombos( hits[0], hits[1], hits[2], 
			    _config.neighborhoods.at(0), _config.neighborhoods.at(1), _config.neighborhoods.at(2), 
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
	  int nbadchs = 0;
	  for ( int p=0; p<3; p++) {
	    if ( badchs.at(p).pixel(r,wirecols[p])>0 ) 
	      nbadchs++;
	  }
	  if (nbadchs>=2) {
	    continue; // combo due two badchs
	  }

	  // otherwise set match

	  // get pos
	  float x = 0;
	  if ( pt==0 || pt==1 ) {
	    // top and bottom use the z value
	    x = combo.pos[2];
	  }
	  else {
	    // up stream and downstream use the y value
	    x = combo.pos[1]+116.0;
	  }
	  int x_i = (int)x;

	  // get total charge
	  float charge = 0.0;
	  bool has_charge[3] = { false, false, false };
	  int nhascharge = 0;
	  for (int p=0; p<3; p++) {
	    if ( imgs.at(p).pixel( r, wirecols[p] )>_config.thresholds.at(p) 
		 || badchs.at(p).pixel( r, wirecols[p] )>0 ) {
	      has_charge[p] = true;
	      nhascharge++;
	      charge += imgs.at(p).pixel( r, wirecols[p] );
	    }
	  }
	  charge /= float( 3.0-nbadchs );
	  if ( charge>_config.thresholds.at(0) ) {
	    float prev_val = matchedspacepts.at( pt ).pixel( r, x_i );
	    matchedspacepts.at(pt).set_pixel( r, x_i, prev_val+charge );
	    for (int p=0; p<3; p++) {
	      prev_val = matchedpixels.at(3*pt+p).pixel( r, wirecols[p] );
	      matchedpixels.at(3*pt+p).set_pixel( r, wirecols[p], prev_val+charge );
	    }
	    // save (x,y) point for clustering
	    std::vector<double> pt_combo(2,0.0); // this is (z,y)
	    pt_combo[0] = x;
	    pt_combo[1] = r + 0.1*float(idx_combo%10);
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
    
    // cluster each boundary type
    std::cout << "  begin clustering" << std::endl;
    clock_t begin_clustering = clock();
    for (int pt=0; pt<4; pt++) {
      
      dbscan::DBSCANAlgo dbalgo;
      //std::cout << "  starting dbscan for boundary type: " << pt << ". number of points: " << combo_points[pt].size() << std::endl;
      if ( combo_points[pt].size()<_config.boundary_cluster_minpixels.at(0)) {
	//for (int i=0; i<combo_points[pt].size(); i++) {
	//std::cout << " (" << i << "): " << combo_points[pt].at(i)[0] << ", " << combo_points[pt].at(i)[1] << std::endl;
	//}
	//std::cout << "   not enough points to make a cluster: skipping clustering" << std::endl;
	continue;
      }
      dbscan::dbscanOutput clout = dbalgo.scan( combo_points[pt], _config.boundary_cluster_minpixels.at(0), _config.boundary_cluster_radius.at(0), false, 0.0 );
      //ann::ANNAlgo::cleanup();

      //std::cout << "  number of clusters: " << clout.clusters.size() << std::endl;
      for (int ic=0; ic<clout.clusters.size(); ic++) {

	if ( clout.clusters.at(ic).size() > 2 ) {
	  //std::cout << "Find the endpoints for cluster pt=" << pt << " id=" << ic << " size=" << clout.clusters.at(ic).size() << std::endl;

	  dbscan::dbPoints chargepts[3]; // charge per plane
	  int largest_qpt[3] = {0,0,0};
	  float larget_q[3]  = {0,0,0};

	  // loop through hit in real space cluster. collect pixels in image space to cluster once more.
	  for (int ihit=0; ihit<clout.clusters.at(ic).size(); ihit++) {
	    int idxhit = clout.clusters.at(ic).at(ihit);
	    for (int p=0; p<3; p++) {
	      int col = combo_cols[pt][idxhit][p];
	      for (int n=-_config.neighborhoods[p]; n<=_config.neighborhoods[p]; n++) {
		if ( col+n<0 || col+n>=imgs.at(p).meta().cols() ) continue;
		float q = imgs.at(p).pixel( (int)combo_points[pt][idxhit][1], col+n );
		if ( q > _config.thresholds.at(p) ) {
		  std::vector<double> qpt(2,0.0);
		  qpt[0] = col+n + 0.1*rand.Uniform();
		  qpt[1] = combo_points[pt][idxhit][1]+0.1*rand.Uniform();
		  chargepts[p].emplace_back( std::move( qpt ) );
		  if ( q>larget_q[p] ) {
		    largest_qpt[p] = chargepts[p].size()-1;
		    larget_q[p] = q;
		  } 
		}//if pixel above thresh
	      }// loop over neighborhood
	    }//end of loop over plane
	  }//end of loop over hits in real space
	  
	  // define the endpoint
	  std::vector< BoundaryEndPt > endpt_v_min; // at min time of cluster
	  std::vector< BoundaryEndPt > endpt_v_max; // at max time of cluster
	  int end_vote[3] = { 0, 0, 0 }; // votes by each plane for min or max time

	  for (int p=0; p<3; p++) {
	    // cluster the charge hits on the plane
	    //std::cout << "  " << chargepts[p].size() << " charge points on plane=" << p << std::endl;
	    if ( chargepts[p].size()< _config.boundary_cluster_minpixels.at(p)+1 ) {
	      // don't have enough here
	      //std::cout << "  not enough to cluster. just make a point from largest charge" << std::endl;
	      BoundaryEndPt endpt( chargepts[p].at( largest_qpt[p] )[1],  chargepts[p].at( largest_qpt[p] )[0], (BoundaryEndPt::BoundaryEnd_t)pt );
	      endpt_v_min.push_back( endpt );
	      endpt_v_max.push_back( endpt );
	      continue;
	    }
 	    //for (int iq=0; iq<chargepts[p].size(); iq++)
	    //  std::cout << "  (" << iq << ") " << chargepts[p].at(iq)[0] << ", " << chargepts[p].at(iq)[1] << std::endl;
	    
	    dbscan::dbscanOutput q_clout = dbalgo.scan( chargepts[p], _config.boundary_cluster_minpixels.at(p), _config.boundary_cluster_radius.at(p), false, 0.0 );
	    //ann::ANNAlgo::cleanup();
	    int largest_cluster = -1;
	    int largest_size = 0;
	    for ( int icq=0; icq<q_clout.clusters.size(); icq++ ) {
	      if ( largest_cluster<0 || q_clout.clusters.at(icq).size()>largest_size ) {
		largest_cluster = icq;
		largest_size = q_clout.clusters.at(icq).size();
	      }
	    }

	    int idxhit_tmin, idxhit_tmax, idxhit_wmin, idxhit_wmax;
	    getClusterEdges( chargepts[p], imgs, q_clout, largest_cluster, idxhit_tmin, idxhit_tmax, idxhit_wmin, idxhit_wmax );

	    // the end points form a line, we follow outward on those lines to see which one is more likely the end of the track
	    float pos[2][2] = { { (float)chargepts[p].at( idxhit_tmin )[0], (float)chargepts[p].at( idxhit_tmin )[1] },
				{ (float)chargepts[p].at( idxhit_tmax )[0], (float)chargepts[p].at( idxhit_tmax )[1] } };
	    float dir[2] = { pos[1][0]-pos[0][0], pos[1][1]-pos[0][1] };
	    float norm = 0;
	    for (int i=0; i<2; i++) {
	      norm += dir[i]*dir[i];
	    }
	    norm = sqrt(norm);
	    for (int i=0; i<2; i++)
	      dir[i] /= norm;
	    int nsteps[2] = { 0, 0 }; // steps before no more charge
	    int nzeros[2] = { 0, 0 };
	    for (int j=0; j<2; j++) {
	      for (int i=0; i<40; i++) {
		if ( nzeros[j]>3 ) break; // cut it off
		int r = pos[j][1] + (2*j-1)*dir[1]*i;
		int c = pos[j][0] + (2*j-1)*dir[0]*i;
		if ( r<2 || r>=imgs.at(p).meta().rows()-2 || c<2 || c>=imgs.at(p).meta().cols()-2 )
		  break;
		int npixs = 0;
		for ( int dr=-2; dr<=2; dr++) {
		  for (int dc=-2; dc<=2; dc++) {
		    if ( imgs.at(p).pixel( r+dr, c+dc )> _config.thresholds.at(p) 
			 || badchs.at(p).pixel( r+dr, c+dc ) > 0 )
		      npixs++;
		  }
		}
		if (npixs>0) {
		  nsteps[j]++;
		  nzeros[j] = 0;
		}
		else {
		  nzeros[j]++;
		}
	      }//end of loop over steps
	    }//end of loop over min/max
	    
	    std::cout << "endpt_type=" << pt << " plane=" << p
		      << " (" << imgs.at(p).meta().pos_x( (int)pos[0][0] ) << "," << imgs.at(p).meta().pos_y( (int)pos[0][1] ) << ") "
		      << "vs (" << imgs.at(p).meta().pos_x( (int)pos[1][0] ) << "," << imgs.at(p).meta().pos_y( (int)pos[1][1] ) << "): "
		      << "nsteps[0]=" << nsteps[0] << " vs. nsteps[1]=" << nsteps[1]
		      << std::endl;

// 	    if ( nsteps[0]<nsteps[1] ) {
// 	      BoundaryEndPt endpt( (int)pos[0][1], (int)pos[0][0], (BoundaryEndPt::BoundaryEnd_t)pt );
// 	      endpt_v.emplace_back( endpt );
// 	    }
// 	    else {
// 	      BoundaryEndPt endpt( (int)pos[1][1], (int)pos[1][0], (BoundaryEndPt::BoundaryEnd_t)pt );
// 	      endpt_v.emplace_back( endpt );
// 	    }


	    BoundaryEndPt endpt_min( (int)pos[0][1], (int)pos[0][0], (BoundaryEndPt::BoundaryEnd_t)pt );
	    endpt_min.dir[0] = -dir[0];
	    endpt_min.dir[1] = -dir[1];
	    endpt_v_min.emplace_back( endpt_min );
	    BoundaryEndPt endpt_max( (int)pos[1][1], (int)pos[1][0], (BoundaryEndPt::BoundaryEnd_t)pt );
	    endpt_max.dir[0] = dir[0];
	    endpt_max.dir[1] = dir[1];
	    endpt_v_max.emplace_back( endpt_max );
	    
	    // vote for the end point
	    end_vote[p] = nsteps[0]-nsteps[1];
	  }//end of loop over plane
	  // store it
	  int totvote = 0;
	  for (int p=0; p<3; p++) totvote += end_vote[p];
	  if ( totvote<0 )
	    end_points.emplace_back( std::move(endpt_v_min) );
	  else
	    end_points.emplace_back( std::move(endpt_v_max) );
	  
	}//end of if cluster size is large enough
      }//end of cluster loop
    }//end of boundary point type
    float elapsed_clustering = float( clock()-begin_clustering )/CLOCKS_PER_SEC;

    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "boundary pixel search took " << elapsed_secs << " secs" << std::endl;
    std::cout << "  hit collecting time: " << tot_hit_collecting << " secs" << std::endl;
    std::cout << "  clustering time: " << elapsed_clustering << " secs" << std::endl;
    std::cout << "  total number of combos found: " << total_combos << std::endl;

    
    return kOK;
  }

  int BoundaryMuonTaggerAlgo::makeTrackClusters3D( std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
						   const std::vector< const std::vector< BoundaryEndPt >* >& spacepts,
						   std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters ) {
    
    // wrap references into a struct so i can treat them more like an array
    // the order matches the order in BoundaryEnd_t enum
    int nendpts = (int)spacepts.size();

    // pair up containers
    int ntotsearched = 0;
    int npossible = 0;
    // poor-man's profiling
    const clock_t begin_time = clock();
    std::vector<int> modimg(nendpts,0);
    for (int i=0; i<nendpts; i++) {
      for (int j=i+1; j<nendpts; j++) {
	const std::vector< BoundaryEndPt >& pts_a = *(spacepts[i]);
	const std::vector< BoundaryEndPt >& pts_b = *(spacepts[j]);

	if ( pts_a.at(0).type==pts_b.at(0).type ) continue; // don't connect same type
	npossible++;
	std::cout << "[ path-finding for endpoints (" << i << "," << j << ") "
		  << "of type (" << pts_a.at(0).type << ") -> (" << pts_b.at(0).type << ") ]" << std::endl;

	// we modify the image in order to mitigate against space-charge effects
	/*
	const std::vector< BoundaryEndPt >* pts[2] = { &pts_a, &pts_b };
	int pt_idx[2] = { i, j };
	for (int ipt=0; ipt<2; ipt++) {
	  if ( modimg.at( pt_idx[ipt] )==1 ) continue; // don't remod for this point
	  modimg.at( pt_idx[ipt] ) = 1; // mark that mod made
	  BoundaryEndPt::BoundaryEnd_t pt_type = (*pts[ipt])[0].type;
	  if ( pt_type==BoundaryEndPt::kTop ) {
	    // we add to the right for plane=0, to the left in plane=1
	    for (int n=-5; n<=50; n++) {
	      float prev = img_v.at(0).pixel( (*pts[ipt])[0].t, (*pts[ipt])[0].w+n );
	      img_v.at(0).set_pixel( (*pts[ipt])[0].t, (*pts[ipt])[0].w+n, prev+20.0 );
	    }
	    for (int n=5; n>=-50; n-- ) {
	      float prev = img_v.at(1).pixel( (*pts[ipt])[1].t, (*pts[ipt])[1].w+n );
	      img_v.at(1).set_pixel( (*pts[ipt])[1].t, (*pts[ipt])[1].w+n, prev+20.0 );
	    }
	  }
	  else if ( pt_type==BoundaryEndPt::kBottom ) {
	    // we add to the right for plane=0, to the left in plane=1
	    for (int n=-5; n<=50; n++) {
	      float prev = img_v.at(1).pixel( (*pts[ipt])[1].t, (*pts[ipt])[1].w+n );
	      img_v.at(1).set_pixel( (*pts[ipt])[1].t, (*pts[ipt])[1].w+n, prev+20.0 );
	    }
	    for (int n=5; n>=-50; n-- ) {
	      float prev = img_v.at(0).pixel( (*pts[ipt])[0].t, (*pts[ipt])[0].w+n );
	      img_v.at(0).set_pixel( (*pts[ipt])[0].t, (*pts[ipt])[0].w+n, prev+20.0 );
	    }
	  }
	}//end of loop over pts
	*/

	for (int p=0; p<3; p++) {
	  int col_a = pts_a.at(p).w;
	  int row_a = pts_a.at(p).t;
	  int col_b = pts_b.at(p).w;
	  int row_b = pts_b.at(p).t;
	  std::cout << "  plane=" << p << ": "
		    << " (w,t): (" << img_v.at(p).meta().pos_x( col_a ) << ", " << img_v.at(p).meta().pos_y( row_a ) << ") ->"
		    << " (" << img_v.at(p).meta().pos_x( col_b ) << "," << img_v.at(p).meta().pos_y( row_b ) << ")"
		    << std::endl;	    
	}
	
	//use test heuristic to see if we should run astar
	bool shallwe = passTrackTest( pts_a, pts_b, img_v, badchimg_v );
	
	if ( !shallwe ) { 
	  std::cout << "failed heuristic." << std::endl;
	  continue; // we shant
	}
	
	std::vector< BMTrackCluster2D > planetracks;
	int ncompleted = 0;
	for (int p=0; p<3; p++) {
// 	  std::cout << "  plane=" << p << ": "
// 		    << " (c,r): (" << pts_a.at(p).w << "," << pts_a.at(p).t << ") ->"
// 		    << " (" << pts_b.at(p).w << "," << pts_b.at(p).t << "), "
// 		    << " (w,t): (" << img_v.at(p).meta().pos_x( col_a ) << ", " << img_v.at(p).meta().pos_y( row_a ) << ") ->"
// 		    << " (" << img_v.at(p).meta().pos_x( col_b ) << "," << img_v.at(p).meta().pos_y( row_b ) << ")"
// 		    << std::endl;
	  
	  BMTrackCluster2D track = runAstar( pts_a.at(p), pts_b.at(p), img_v.at(p), badchimg_v.at(p), 5, 5, 2, true );
	  std::cout << "  p=" << p << ": pathsize=" << track.pixelpath.size() << std::endl;
	  if ( track.pixelpath.size()>3 ) ncompleted++;
	  planetracks.emplace_back( std::move( track ) );
	}
	
	if ( ncompleted>=2 ) {
	  trackclusters.emplace_back( std::move(planetracks) );
	}
	ntotsearched++;
      }//loop over pt b
    }//loop over pt a
    
    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "total paths searched: " << ntotsearched << " of " << npossible << " possible combinrations. time=" << elapsed_secs << " secs" << std::endl;
    
    return 0;
  }

  int BoundaryMuonTaggerAlgo::markImageWithTrackClusters( const std::vector<larcv::Image2D>& imgs, 
							  const std::vector< std::vector< larlitecv::BMTrackCluster2D > >& trackclusters,
							  std::vector<larcv::Image2D>& markedimgs ) {
    for ( auto &img : imgs ) {
      const larcv::ImageMeta& meta = img.meta();
      larcv::Image2D markedimg( meta );
      markedimg.paint(0.0);

      int plane = (int)meta.plane();
      
      const std::vector< larlitecv::BMTrackCluster2D >& planetracks = trackclusters.at(plane);

      for ( auto &track : planetracks ) {
	for ( auto &pixel : track.pixelpath ) {
	  int col = pixel.X();
	  int row = pixel.Y();

	  for ( int dc=-_config.astar_neighborhood.at(plane); dc<=_config.astar_neighborhood.at(plane); dc++ ) {
	    for ( int dr=-_config.astar_neighborhood.at(plane); dr<=_config.astar_neighborhood.at(plane); dr++ ) {
	      int r = row+dr;
	      int c = col+dc;
	      if ( r<0 || r>=meta.rows() ) continue;
	      if ( c<0 || c>=meta.cols() ) continue;
	      float val = img.pixel( row, col );
	      if ( val>_config.thresholds.at(plane) ) {
		markedimg.set_pixel( r, c, 100.0 );
	      }
	    }
	  }
	}//end of pixel list
      }//end of planetracks
      markedimgs.emplace_back( markedimg );
    }
    
    return 0;
  }

  BMTrackCluster2D BoundaryMuonTaggerAlgo::runAstar( const BoundaryEndPt& start_pt, const BoundaryEndPt& end_pt, 
						     const larcv::Image2D& img, const larcv::Image2D& badchimg,
						     int start_pad, int end_pad, int verbose, bool use_badchs ) {
    // This wraps/interfaces with the AStar algorithm on our image. runs one start to end check.
    // inputs
    // ------
    //  start: starting boundary point
    //  end: ending boundary point
    //  img: image with charge deposition we will follow
    //  start_pad: radius around starting point where we will include all hits in case not start point directly on track
    //  end_pad:   radius around ending point where we will include all hits in case end point not directly on track

    
    larlitecv::AStarAlgoConfig astar_config;
    astar_config.astar_threshold    = _config.astar_thresholds;
    astar_config.astar_neighborhood = _config.astar_neighborhood;
    astar_config.astar_start_padding = start_pad;
    astar_config.astar_end_padding   = end_pad;
    larlitecv::AStarGridAlgo algo( astar_config );
    algo.setVerbose(verbose);
    algo.setBadChImage(badchimg);

    larlitecv::AStarNode start( start_pt.w, start_pt.t );
    larlitecv::AStarNode goal( end_pt.w, end_pt.t );
    
    BMTrackCluster2D track2d;
    track2d.start = start_pt;
    track2d.end   = end_pt;
    track2d.pixelpath.clear();

    if ( abs(start.row-goal.row)>820 ) { // make this a parameter 
      return track2d; // empty track
    }
    
    std::vector< larlitecv::AStarNode > path = algo.findpath( img, start.row, start.col, goal.row, goal.col, 5.0, use_badchs );
    for ( auto& node : path ) {
      larcv::Pixel2D pixel( node.col, node.row );
      pixel.Intensity( img.pixel( node.row, node.col ) );
      track2d.pixelpath += pixel; // uniary operator!
    }
    return track2d;
  }

  std::vector<BMTrackCluster2D> BoundaryMuonTaggerAlgo::runAstar3planes( const std::vector< BoundaryEndPt >& start_pts, const std::vector< BoundaryEndPt >& end_pts,
									 const std::vector< larcv::Image2D >& img, int start_pad, int end_pad ) {
    std::vector<BMTrackCluster2D> planetracks;
    return planetracks;
  }


  /*
  void BoundaryMuonTaggerAlgo::filterShortAndKinkedTracks( int mintracklength, std::vector< larlitecv::BMTrackCluster2D >& tracks, 
							   std::vector< larlitecv::BMTrackCluster2D >& filtered, bool savefiltered ) {
    // this removes tracks with large kinks
    // input:
    //   mintracklength: filter out tracks shorter than this
    //   input: a list of bmtrackclusters.  we will remove elements
    //   output: filtered tracks. will be filled based on 'savefiltered' values
    //   savefiltered (optional): if true, passes filtered tracks into output container

    std::vector<int> remove_index_list;
    
    // loop over input tracks
    for ( auto &track : tracks ) {
      continue;
    }
    
  }
  */
  
  void BoundaryMuonTaggerAlgo::matchTracksStage1( const std::vector< larcv::Image2D >& imgs,
						  const std::vector< std::vector< larlitecv::BMTrackCluster2D >* >& plane2dtracks, 
						  std::vector< larlitecv::BMTrackCluster3D >& output  ) {
    // we do dumb O(N3) matching search.  
    // in the future we could structure start and end points into a tree for more efficient search: to make it Nlog(N)
    // also data products are starting to become heavy, might want to start using storage container and refer to indicices
    //     we'll see when things start to get bogged down.

    std::set< larlitecv::BMTrackCluster3D > track_set;

    for (int p1=0; p1<1; p1++) {
      for (int p2=p1+1; p2<2; p2++) {
	
	const std::vector< larlitecv::BMTrackCluster2D >& p1tracks = *(plane2dtracks.at(p1));
	const std::vector< larlitecv::BMTrackCluster2D >& p2tracks = *(plane2dtracks.at(p2));

	for ( int idx1=0; idx1<p1tracks.size(); idx1++) {
	  const larlitecv::BMTrackCluster2D& track1 = p1tracks.at(idx1);
	  std::vector< int > match_indices;
	  std::vector< float > tot_start_diff;
	  std::vector< float > tot_end_diff;
	  for ( int idx2=0; idx2<p2tracks.size(); idx2++) {
	    const larlitecv::BMTrackCluster2D& track2 = p2tracks.at(idx2);
	    float sdiff,ediff;
	    bool start2start;
	    if ( !doTracksMatch( track1, track2, sdiff, ediff, start2start ) ) 
	      continue; // move on
	    
	    // check the third plane. what is it?
	    int p3 = -1;
	    if ( p1==0 && p2==1 ) p3 = 2;
	    else if ( p1==0 && p2==2 ) p3 = 1;
	    else if ( p1==1 && p2==2 ) p3 = 0;
	    if ( p3==-1 ) {
	      std::cout << "[BoundaryMuonTagger::matchTracksStage2] Invalid plane combination." << std::endl;
	      assert( false );
	    }

	    // look for third plane match
	    const std::vector< larlitecv::BMTrackCluster2D >& p3tracks = *(plane2dtracks.at(p3));
	    for ( int idx3=0; idx3<(int)p3tracks.size(); idx3++ ) {
	      const larlitecv::BMTrackCluster2D& track3 = p3tracks.at(idx3);
	      
	      float p1sdiff, p2sdiff, p1ediff, p2ediff;
	      bool p1start2start, p2start2start;
	      
	      if ( doTracksMatch( track1, track3, p1sdiff, p1ediff, p1start2start ) && doTracksMatch( track2, track3, p2sdiff, p2ediff, p2start2start ) ) {
		
		std::cout << " match on planes=(" << p1 << "," << p2 << "," << p3 << ") tracks=(" << idx1 << "," << idx2 << "," << idx3 << ")" << std::endl;

		// wow, made it all this way
		larlitecv::BMTrackCluster3D track3d;
		//
		track3d.trackidx.resize(3);
		track3d.trackidx.at(0) = idx1;
		track3d.trackidx.at(1) = idx2;
		track3d.trackidx.at(2) = idx3;

		// start/end time is average of all three
		int t_start  = track1.start.t;
		t_start += ( (start2start) ?   track2.start.t : track2.end.t );
		t_start += ( (p1start2start) ? track3.start.t : track3.end.t );
		t_start /= 3;

		int t_end = track1.end.t;
		t_end += ( (start2start) ?   track2.end.t : track2.start.t );
		t_end += ( (p1start2start) ? track3.end.t : track3.start.t );
		t_end /= 3;
		
		track3d.row_start = t_start;
		track3d.row_end   = t_end;
		track3d.tick_start = imgs.at(0).meta().pos_y( track3d.row_start );
		track3d.tick_end   = imgs.at(0).meta().pos_y( track3d.row_end );
		
		// set type
		track3d.start_type = track1.start.type;
		track3d.end_type   = track1.end.type;

		// set start wires
		track3d.start_wire.resize(3);
		track3d.start_wire.at(0) = track1.start.w;
		track3d.start_wire.at(1) = ( (start2start) ?   track2.start.w : track2.end.w );
		track3d.start_wire.at(2) = ( (p1start2start) ? track3.start.w : track3.end.w );

		// set end wires
		track3d.end_wire.resize(3);
		track3d.end_wire.at(0) = track1.end.w;
		track3d.end_wire.at(1) = ( (start2start) ?   track2.end.w : track2.start.w );
		track3d.end_wire.at(2) = ( (p1start2start) ? track3.end.w : track3.start.w );
		
		// boundary intersection algo
		larlitecv::BoundaryIntersectionAlgo algo;
		algo.determine3Dpoint( track3d.start_wire, track3d.start3D, track3d.start_type );
		algo.determine3Dpoint( track3d.end_wire, track3d.end3D, track3d.end_type );
		
		track_set.insert( track3d );
		output.emplace_back( track3d );
		
	      }
	    }
	    
	  }//end of loop over plane 2
	}
      }//end of loop over plane2
    }//end of loop over plane3

    std::cout << "Number of tracks in the set: " << track_set.size() << std::endl;

  }


  bool BoundaryMuonTaggerAlgo::doTracksMatch( const larlitecv::BMTrackCluster2D& track1, const larlitecv::BMTrackCluster2D& track2, 
					      float& start_t_diff, float& end_t_diff, bool& start2start ) {
    // inputs:
    //  track1: a TrackCluster2D
    //  track2: a TrackCluster2D
    // outputs:
    //  start_t_diff: difference in start point times
    //  end_t_diff: difference in end point times
    //  start2start: true if start point to start point matches. false if start to end points match.

    // does the first end point match in type?
    start2start = true; // 0=start; 1=end
    if ( track1.start.type == track2.start.type )
      start2start = true;
    else if ( track1.start.type!=track2.end.type )
      start2start = false;
    else
      return false; // does not match either
    
    // check if start_matches in time
    start_t_diff = 1e9;
    if ( start2start )
      start_t_diff = abs( track1.start.t - track2.start.t );	      
    else if ( start2start==1 )
      start_t_diff = abs( track1.start.t - track2.end.t );
    
    if ( start_t_diff > 10 )
      return false; // move on to the next one
    
    // does the end point match in type and time?
    end_t_diff = 1e9;
    bool end_matches_in_type = false;
    if ( start2start ) {
      end_t_diff = abs( track1.end.t - track2.end.t );
      if ( track1.end.type==track2.end.type ) end_matches_in_type = true;
    }
    else {
      end_t_diff = abs( track1.end.t - track2.start.t );
      if ( track1.end.type==track2.start.type ) end_matches_in_type = true;
    }
    
    if ( end_t_diff > 10 || !end_matches_in_type ) 
      return false;
    
    // profit!!!
    return true;
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
    
  }

  void BoundaryMuonTaggerAlgo::getClusterEdges( const dbscan::dbPoints& points,  const std::vector< larcv::Image2D >& imgs,
						const dbscan::dbscanOutput& clout, int idx_cluster,
						int& idxhit_tmin, int& idxhit_tmax, int& idxhit_wmin, int& idxhit_wmax ) {
    idxhit_tmin = -1;
    idxhit_tmax = -1;
    idxhit_wmin = -1;
    idxhit_wmax = -1;

    double hit_tmin[2];
    double hit_tmax[2];
    double hit_wmin[2];
    double hit_wmax[2];
		      
    // we find the extrema hits in a cluster
    for (int ichit=0; ichit<(int)clout.clusters.at(idx_cluster).size(); ichit++) {
      int hitidx = clout.clusters.at(idx_cluster).at(ichit);
      const std::vector<double>& hit = points.at(hitidx);
      //float hit_z = hit.at(0);
      //int ycol    = cols.at(2);
      //int row     = (int)hit.at(1);

      //if ( imgs.at(2).pixel( row, ycol )<10.0 ) continue;

      // tmin hit
      if ( idxhit_tmin==-1 || hit[1]<hit_tmin[1] ) {
	idxhit_tmin = hitidx;
	hit_tmin[0] = hit[0];
	hit_tmin[1] = hit[1];
      }
      // tmax hit
      if ( idxhit_tmax==-1 || hit[1]>hit_tmax[1] ) {
	idxhit_tmax = hitidx;
	hit_tmax[0] = hit[0];
	hit_tmax[1] = hit[1];
      }
//       // wmin hit
//       if ( idxhit_wmin==-1 || hit[0]<hit_wmin[0] ) {
// 	hit_wmin[0] = hit[0];
// 	hit_wmin[1] = hit[1];
//       }
//       //  wmax hit
//       if ( idxhit_wmax==-1 || hit[0]>hit_wmin[0] ) {
// 	hit_wmax[0] = hit[0];
// 	hit_wmax[1] = hit[1];
//       }
    }//end of loop over hit indices of cluster
  }//end of getClusterEdges

  void BoundaryMuonTaggerAlgo::calcTrackTest( const BoundaryEndPt& start, const BoundaryEndPt& end, 
					      const larcv::Image2D& img,  const larcv::Image2D& badchimg, 
					      float angle, float pix_thresh, int time_win, int wire_win,
					      std::vector<float>& q_in_angle, std::vector<int>& pixels_in_angle, std::vector<int>& badpixs_in_angle ) {
    // we test the pixels in the neighborhood of the end points
    // if they have enough pixels going in the direction to the other end point, the test passes

    q_in_angle.resize(2,0.);
    pixels_in_angle.resize(2,0);
    badpixs_in_angle.resize(2,0);


    float dir[2] = { (float)(end.w-start.w), (float)(end.t-start.t) };
    float dist = 0.;
    for (int i=0; i<2; i++)
      dist  += dir[i]*dir[i];
    dist = sqrt(dist);
    if (dist==0) // wtf
      return;
    for (int i=0; i<2; i++)
      dir[i] /= dist;
      
    const larcv::ImageMeta& meta = img.meta();
    
    const BoundaryEndPt* pts[2] = { &start, &end };
    
    float cos_angle = cos(angle*3.14159/180.0);

    for (int pt=0; pt<2; pt++) {
      const BoundaryEndPt& endpt = *(pts[pt]);
      for (int dt=-time_win; dt<=time_win; dt++) {
	if ( endpt.t+dt<0 || endpt.t+dt>=meta.rows() ) continue;
	for (int dw=-wire_win; dw<=wire_win; dw++) {
	  if ( endpt.w+dw<0 || endpt.w+dw>=meta.cols() ) continue;
	  if ( dt==0 && dw==0 ) continue;
	  int t=endpt.t+dt;
	  int w=endpt.w+dw;
	  bool isbad = false;
	  if ( badchimg.pixel(t,w)>0 ) isbad = true;
	  if ( img.pixel( t, w )<pix_thresh || isbad ) continue;
	  // direction to pixel
	  float pixnorm = sqrt( dw*dw+dt*dt );
	  float pixdir[2] = { dw/pixnorm, dt/pixnorm };
	  float pixcos = pixdir[0]*dir[0] + pixdir[1]*dir[1];
	  if ( pixcos*(1.0-pt*2.0) > 1.0-cos_angle  ) {
	    if ( !isbad ) {
	      q_in_angle[pt] += img.pixel( t, w );
	      pixels_in_angle[pt]++;
	    }
	    else {
	      badpixs_in_angle[pt]++;
	    }
	  }//end of within target cos
	}//end of wire region loop
      }//end of time region loop
    }//end of pt loop
  }

  bool BoundaryMuonTaggerAlgo::passTrackTest( const std::vector<BoundaryEndPt>& start_v, const std::vector<BoundaryEndPt>& end_v,
                                              const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v ) {
    
    // we first get track test results for each plane
    bool passes_test[3] =   { false, false, false };
    bool passes_wbadch[3] = { false, false, false };
    float angle = 10.0;
    float pix_thresh = 10.0;
    int time_win = 10;
    int wire_win = 10;
    std::cout << "test heuristic: ";
    for (int p=0; p<3; p++) {
      const BoundaryEndPt& start = start_v.at(p);
      const BoundaryEndPt& end   = end_v.at(p);
      const larcv::Image2D& img  = img_v.at(p);
      const larcv::Image2D& bad  = badchimg_v.at(p);
      std::vector<float> q_in_angle;
      std::vector<int> pixels_in_angle;
      std::vector<int> badpixs_in_angle;
      calcTrackTest( start, end, img, bad, angle, pix_thresh, time_win, wire_win, q_in_angle, pixels_in_angle, badpixs_in_angle );
      std::cout << "  p=" << p 
		<< "  q=(" << q_in_angle[0] << "," << q_in_angle[1] << ")"
		<< "  pix=(" << pixels_in_angle[0] << "," << pixels_in_angle[1] << ")"
		<< "  bad=(" << badpixs_in_angle[0] << "," << badpixs_in_angle[1] << ")"
		<< std::endl;
      if ( pixels_in_angle[0]>8 && pixels_in_angle[1]>8 ) {
	passes_test[p] = true;
	if ( badpixs_in_angle[0]>0 || badpixs_in_angle[1]>0 ) 
	  passes_wbadch[p] = true;
      }
    }
    
    int npasses = 0;
    for ( int p=0; p<3; p++) {
      if ( passes_test[p] )
	npasses++;
    }

    if ( npasses>=2 )
      return true;
    return false;
  }
  
}
