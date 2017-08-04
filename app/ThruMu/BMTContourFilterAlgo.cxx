#include "BMTContourFilterAlgo.h"
#include <assert.h>
#include <sstream>
#include <exception>

#include "UBWireTool/UBWireTool.h"

#ifdef USE_OPENCV
#include "CVUtil/CVUtil.h"
#endif


namespace larlitecv {

  // =======================================================================================================
  BMTContourFilterAlgo::BMTContourFilterAlgo() {};

  BMTContourFilterAlgo::~BMTContourFilterAlgo() {};

  bool BMTContourFilterAlgo::buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
					   const std::vector<float>& pos3d, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					   const float max_dist2contour ) {
    // we rely on contours produced by BMTCV::analyzeImages + BMTCV::splitcontours
    
    // a blank image to build merged cluster
    if ( clusterpix_v.size()==0 ) {
      for (int p=0; p<(int)img_v.size(); p++) {
	larcv::Image2D img( img_v[p].meta() );
	img.paint(0.0);
	clusterpix_v.emplace_back( std::move(img) );
      }
    }
    
    std::vector<int> imgcoords;
    try {
      imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos3d, img_v.front().meta(), img_v.size() );
    }
    catch (...) {
      std::cout << __FILE__ << ":" << __LINE__ << " Spacepoint could not project into the image." << std::endl;
      return false;
    }
    
    // first check if boundary point is in a contour past?
    // ----------------------------------------------    
    // std::cout << __FILE__ << ":" << __LINE__ << " "
    // 	      << " pos=" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2]
    // 	      << "  imgcoords=(" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")"
    // 	      << std::endl;

    // int nplanes_in_contour = 0;
    // for (size_t p=0; p<img_v.size(); p++) {
    //   int row = imgcoords[0];
    //   int col = imgcoords[p+1];
    //   if ( clusterpix_v[p].pixel(row,col)>0 )
    // 	nplanes_in_contour++;
    // }


    // if ( nplanes_in_contour==(int)img_v.size() ) {
    //   std::cout << __FILE__ << ":" << __LINE__ << " already in a contour. skipping." << std::endl;
    //   return;
    // }

    // Convert Image2D position into a cv::Point
    std::vector<cv::Point> imgpt;
    for (size_t p=0; p<img_v.size(); p++) {
      cv::Point pt( imgcoords[p+1], imgcoords[0] );
      imgpt.emplace_back( std::move(pt) );
    }


    cuts.resize(2,false);
    
    // ok now we need the seed cluster
    ContourCluster foundcluster;
    cuts[0] = isPointInContour( imgpt, img_v, plane_contours_v, max_dist2contour, foundcluster );

    if ( !cuts[0] )
      return false;
    
    // seed cluster analysis. for the lazy
    cuts[1]  = analyzeSeedContours( imgpt, img_v, 10.0, foundcluster );
    
    // now we extend as far as we can
    /*
    std::vector<cv::Mat> cvimg_v;
    for ( auto const& img : img_v ) {
      cv::Mat cvimg = larcv::as_gray_mat( img, 8.0, 256.0, 1.0 );
      cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img, 10.0, 100.0 );
      cvimg_v.emplace_back( std::move(cvrgb) );
    }

    // perform the expansion loop    
    int maxnsteps = 20;
    int iter = 0;
    bool extended = true;
    while ( extended && iter<maxnsteps ) {
      std::cout << "Expansion " << iter << " ------------------------------------" << std::endl;
      bool extendbyflailing = false;
      extended = ratchetCluster( foundcluster, plane_contours_v, img_v, badch_v, clusterpix_v );
      if ( !extended  ) {
	extendbyflailing = true;
	extended = extendClusterGroup( foundcluster, plane_contours_v, img_v, badch_v, clusterpix_v );
      }
 
      iter++;
    }//end of expansion loop

    std::cin.get();
    */

    /*
    // for debug/visualize
    for (int p=0; p<3; p++) {

      for (int p=0; p<3; p++) {
	cv::Mat& cvrgb = cvimg_v[p];
	if ( cluster.earlyContours[p].size()>=2 ) {
	  int last     = cluster.earlyContours[p].size()-1;
	  int nextlast = cluster.earlyContours[p].size()-2;
	  cv::line( cvrgb, cluster.earlyContours[p][last]->getFitSegmentEnd(), cluster.earlyContours[p][nextlast]->getFitSegmentStart(), cv::Scalar(255,0,0,255), 1 );
	}	

	if ( cluster.earlyContours[p].size()>0 ) {
	  auto const& contour = *(cluster.earlyContours[p].back());
	  std::vector< std::vector<cv::Point> >  contour_v;
	  contour_v.push_back( contour );
	  cv::Scalar contourcolor(0,0,255,255);
	  if ( extendbyflailing )
	    contourcolor = cv::Scalar(255,0,255,255);
	  cv::drawContours( cvrgb, contour_v, 0, contourcolor, -1 );
	  cv::circle( cvrgb, cluster.lateEnd[p], 3, cv::Scalar(0,255,0,255), 1 );	
	}
	auto& earlypt = cluster.earlyEnd[p].back();
	cv::circle( cvrgb, earlypt, 3, cv::Scalar(0,255,0,255), 1 );	  	
      }

      
      std::stringstream imgname;
      imgname << "contourclusterdev_plane" << p << ".png";
      cv::imwrite( imgname.str(), cvimg_v[p] );
    }
    */
    
    // graph attempt
    //buildContourGraph( cluster, plane_contours_v, img_v, badch_v, clusterpix_v );

    
    // fill out cluster, make solid contour

    return cuts[0];
  }

  bool BMTContourFilterAlgo::isPointInContour( const std::vector<cv::Point>& imgpt, const std::vector<larcv::Image2D>& img_v,
					       const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					       const float max_dist2contour,
					       ContourCluster& outcluster ) {

    // ok now we need the seed cluster
    std::vector< const ContourShapeMeta* > plane_seedcontours;
    int nplanes_found = 0;
    std::vector<int> seed_idx(3,-1);
    for (size_t p=0; p<img_v.size(); p++) {

      bool contains = false;
      int containing_idx = -1;
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {
	// test imgpt
	bool bboxcontains = plane_contours_v[p][idx].getBBox().contains( imgpt[p] );
	if ( !bboxcontains )
	  continue;

	// more detailed test
	double dist = cv::pointPolygonTest( plane_contours_v[p][idx], imgpt[p], true );
	if (dist>=-max_dist2contour ) {
	  contains = true;
	  containing_idx = idx;
	}

	if ( contains )
	  break;
      }
	
      if ( contains ) {
	plane_seedcontours.push_back( &(plane_contours_v[p][containing_idx]) );
	nplanes_found++;
      }
      else {
	plane_seedcontours.push_back( NULL );
      }
      seed_idx[p] = containing_idx;
    }//end of loop over planes
    
    // we need at least 2 seed clusters
    if ( nplanes_found<2 ) { // 3 for now. for 2, we need an extension. but let's wait.
      std::cout << __FILE__ << ":" << __LINE__ << " Space point matches to contours in " << nplanes_found << " planes only" << std::endl;      
      return false;
    }

    // make the seed cluster
    // ContourCluster cluster( plane_seedcontours );
    // for (int p=0; p<3; p++)
    //   cluster.indices[p].insert( seed_idx[p] );

    // std::swap( outcluster, cluster );
    outcluster.addEarlyContours( plane_seedcontours );
    
    return true;
  }

  int BMTContourFilterAlgo::getIndexOfContainingContour( const int row, const int col, const std::vector<ContourShapeMeta>& contours_v, int min_cluster_size, float dist_tolerance ) {
    // first test bbox (faster test)
    // then do point poly test
    
    bool contains = false;
    int containing_idx = -1;
    cv::Point imgpt(col,row);
    for (int idx=0; idx<(int)contours_v.size(); idx++) {
      
      if ( contours_v[idx].size()<min_cluster_size )
	continue;
      
      // test imgpt
      bool bboxcontains = contours_v[idx].getBBox().contains( imgpt );
      if ( !bboxcontains )
	continue;
      
      // more detailed test
      double dist = cv::pointPolygonTest( contours_v[idx], imgpt, false );
      if ( dist>=-fabs(dist_tolerance) )
	return idx;
    }
    return -1;
  }  

  bool BMTContourFilterAlgo::ratchetCluster( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v ) {
    // what we do is we find the tmax and tmin of the clusters across planes
    // we then extend the direction of the current contours to that time until we run into another contour.
    // we build until we reach. we keep going, ratcheting until no more clusters.
    
    bool cluster_extended = false;

    float row_max = 0;
    float row_min = (float)img_v.front().meta().rows();

    // we get the time bounds of the current contour cluster
    for ( size_t p=0; p<img_v.size(); p++ ) {

      // get the current ends of the cluster
      const cv::Point& earlypt = cluster.earlyEnd[p].back(); // early end on the plane
      const cv::Point& latept  = cluster.lateEnd[p];         // late end on the plane

      if ( earlypt.y>row_max ) {
	row_max = earlypt.y;
      }
      if ( earlypt.y<row_min ) {
	row_min = earlypt.y;
      }
      if ( latept.y>row_max ) {
	row_max = latept.y;
      }
      if ( latept.y<row_min ) {
	row_min = latept.y;
      }
    }

    std::cout << __FILE__ << ":Ratchet Cluster bounds [" << row_min << "," << row_max << "]" << std::endl;
    
    // now for each plane, we extrapolate out from the ends, to the end extremes
    // we step until we run in the a new contour
    for ( size_t p=0; p<img_v.size(); p++ ) {

      // early extension
      const cv::Point& earlypt = cluster.earlyEnd[p].back();
      std::vector<float> earlydir = cluster.earlyDir[p].back();

      if ( earlypt.y==row_min )
	continue; // not extending in this plane

      if ( earlydir[1]==0 )
	continue; // outside the use case unfortunately... need a separate horizontal extension check
      
      float drow = row_min-earlypt.y;
      int nsteps = drow/(earlydir[1]/10.0);

      int lastrow = earlypt.y;
      int lastcol = earlypt.x;
      std::set<int> onpath_indices_set; // for dupilicate search
      std::vector<int> onpath_indices_v; // for listing in order of intersection
      for (int istep=0; istep<=nsteps; istep++) {
	cv::Point testpt( earlypt.x+(float(istep)*earlydir[0]/10.0), earlypt.y+(float(istep)*earlydir[1]/10.0) );
	if ( (int)testpt.x==lastcol && (int)testpt.y==lastrow )
	  continue;

	// update pt
	lastrow = testpt.y;
	lastcol = testpt.x;

	int idx = getIndexOfContainingContour( testpt.y, testpt.x, plane_contours_v[p], 10, 5.0 );
	if ( idx<0 || onpath_indices_set.find(idx)!=onpath_indices_set.end() )
	  continue;

	onpath_indices_set.insert(idx);
	onpath_indices_v.push_back(idx);
      }

      std::cout << __FILE__ << ":RatchetCluster number of on-path candidates " << onpath_indices_v.size() << std::endl;
      
      // evaluate quality of candidates

      for ( auto& idx : onpath_indices_v ) {

	auto const& early_contour = plane_contours_v[p][idx];

	float cosine = 0;
	for (int i=0; i<2; i++)
	  cosine += early_contour.getStartDir()[i]*earlydir[i];

	std::cout << __FILE__ << ":   candidate [" << idx << "] cosine=" << cosine << std::endl;
	
	if ( cosine>0.8 ) {
	  cluster.earlyContours[p].push_back( &early_contour );
	  cluster.earlyDir[p].push_back( early_contour.getStartDir() );
	  cluster.earlyEnd[p].push_back( early_contour.getFitSegmentStart() );
	  cluster.indices[p].insert( idx );
	  cluster_extended = true;
	  break;
	}
      }
    }//end of p loop

    return cluster_extended;
  }//end of extension method
  
  bool BMTContourFilterAlgo::extendClusterGroup( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
						 const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v ) {
    // for each of the active ends of the cluster, we look for the best match on each plane
    // appending requires that we can build a 3D consistent line segment
    // then we update the cluster
    
    bool cluster_extended = false;
    
    class ClusterMatchScore {
    public:
      ClusterMatchScore( int iidx, double fdist, double fcosine, double ftriscore )
      {
	dist   = fdist;
	cosine = fcosine;
	triscore = ftriscore;
	idx = iidx;
      };
      virtual ~ClusterMatchScore() {};
      double dist;
      double cosine;
      double triscore;
      int idx;
      bool operator<( const ClusterMatchScore& rhs ) const {
	// order by distance
	float nearfardivide = 50.0;
	if ( dist< nearfardivide && rhs.dist>nearfardivide ) {
	  return true;
	}
	else if ( dist>nearfardivide && rhs.dist<nearfardivide )
	  return false;
	else if ( dist<nearfardivide && rhs.dist<nearfardivide ) {
	  if ( cosine < rhs.cosine )
	    return true;
	  else
	    return false;
	}
	else if ( dist>=nearfardivide && rhs.dist>=nearfardivide ) {
	  if ( dist < rhs.dist )
	    return true;
	  else
	    return false;
	}
	
	return false;
      };
      
    };//end of clustermatchscore

    std::vector< std::vector<ClusterMatchScore> > cluster_scores;

    // score is based on triangle formed from end-to-end and intersection of end directions
    for ( size_t p=0; p<img_v.size(); p++ ) {

      const cv::Point& earlypt = cluster.earlyEnd[p].back(); // early end on the plane
      const cv::Point& latept = cluster.lateEnd[p]; // late end on the plane

      std::vector< ClusterMatchScore > early_extension_candidates_v;
      std::vector< ClusterMatchScore >  late_extension_candidates_v;
      
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {

	// we already used this cluster
	if ( cluster.indices[p].find(idx)!=cluster.indices[p].end() )
	  continue;

	// get contour
	const ContourShapeMeta& contour = plane_contours_v[p][idx];

	if ( contour.size()<10 )
	  continue;

	// get closet end to early end of current cluster
	float early2start = sqrt( (contour.getFitSegmentStart().x-earlypt.x)*(contour.getFitSegmentStart().x-earlypt.x)
				  +(contour.getFitSegmentStart().y-earlypt.y)*(contour.getFitSegmentStart().y-earlypt.y) );
	float early2end = sqrt( (contour.getFitSegmentEnd().x-earlypt.x)*(contour.getFitSegmentEnd().x-earlypt.x)
				+(contour.getFitSegmentEnd().y-earlypt.y)*(contour.getFitSegmentEnd().y-earlypt.y) );
	float early2startcos = 0;
	float early2endcos   = 0;
	for (int i=0; i<2; i++) {
	  early2startcos += contour.getStartDir()[i]*cluster.earlyDir[p].back()[i];
	  early2endcos   += contour.getEndDir()[i]*cluster.earlyDir[p].back()[i];
	}
		
	std::vector< cv::Point > earlypair;
	float mindist = 0;
	float mincos  = 0;
	std::vector<float> earlydir(2,0);
	if ( early2start<early2end ) {
	  earlypair.push_back( earlypt );
	  earlypair.push_back( contour.getFitSegmentStart() );
	  mindist = early2start;
	  mincos  = early2startcos;
	  earlydir = contour.getStartDir();
	}
	else {
	  earlypair.push_back( earlypt );
	  earlypair.push_back( contour.getFitSegmentEnd() );
	  mindist = early2end;
	  mincos  = early2endcos;
	  earlydir = contour.getEndDir();	  
	}

	// cos between connecting line and contour direction
	float econnectdir[2];
	econnectdir[0] = earlypair[1].x-earlypair[0].x;
	econnectdir[1] = earlypair[1].y-earlypair[0].y;
	float econnectnorm = 0;
	for (int i=0; i<2; i++)
	  econnectnorm += econnectdir[i]*econnectdir[i];
	econnectnorm = sqrt(econnectnorm);
	float econnectcos = 0;
	for (int i=0; i<2; i++) {
	  econnectdir[i] /= econnectnorm;
	  econnectcos += econnectdir[i]*earlydir[i];
	}

	// save the candidate info
	ClusterMatchScore earlymatchcand( idx, mindist, mincos, econnectcos );
	early_extension_candidates_v.emplace_back( earlymatchcand );


	// get closet end to late end of current cluster
	float late2start = sqrt( (contour.getFitSegmentStart().x-latept.x)*(contour.getFitSegmentStart().x-latept.x)
				  +(contour.getFitSegmentStart().y-latept.y)*(contour.getFitSegmentStart().y-latept.y) );
	float late2end = sqrt( (contour.getFitSegmentEnd().x-latept.x)*(contour.getFitSegmentEnd().x-latept.x)
				+(contour.getFitSegmentEnd().y-latept.y)*(contour.getFitSegmentEnd().y-latept.y) );
	float late2startcos = 0;
	float late2endcos   = 0;
	for (int i=0; i<2; i++) {
	  late2startcos += contour.getStartDir()[i]*cluster.lateDir[p][i];
	  late2endcos   += contour.getEndDir()[i]*cluster.lateDir[p][i];
	}
	
	std::vector< cv::Point > latepair;
	std::vector<float> latedir(2,0);
	float latemindist = 0;
	float latemincos  = 0;
	if ( late2start<late2end ) {
	  latepair.push_back( latept );
	  latepair.push_back( contour.getFitSegmentStart() );
	  latemindist = late2start;
	  latemincos  = late2startcos;
	  latedir = contour.getStartDir();
	  //ClusterMatchScore matchcand( idx, late2start, late2startcos, 0.0 );
	  //late_extension_candidates_v.emplace_back( matchcand );	  	  
	}
	else {
	  latepair.push_back( latept );
	  latepair.push_back( contour.getFitSegmentEnd() );
	  latemindist = late2end;
	  latemincos  = late2endcos;
	  latedir = contour.getEndDir();
	  //ClusterMatchScore matchcand( idx, late2end, late2endcos, 0.0 );
	  //late_extension_candidates_v.emplace_back( matchcand );	  	  
	}

	// cos between connecting line and contour direction
	float lconnectdir[2];
	lconnectdir[0] = latepair[1].x-latepair[0].x;
	lconnectdir[1] = latepair[1].y-latepair[0].y;
	float lconnectnorm = 0;
	for (int i=0; i<2; i++)
	  lconnectnorm += lconnectdir[i]*lconnectdir[i];
	lconnectnorm = sqrt(lconnectnorm);
	float lconnectcos = 0;
	for (int i=0; i<2; i++) {
	  lconnectdir[i] /= lconnectnorm;
	  lconnectcos += lconnectdir[i]*latedir[i];
	}

	// save the candidate info
	ClusterMatchScore latematchcand( idx, latemindist, latemincos, lconnectcos );
	late_extension_candidates_v.emplace_back( latematchcand );
	
      }//end of loop over plane conours

      std::sort( early_extension_candidates_v.begin(), early_extension_candidates_v.end() );
      std::sort( late_extension_candidates_v.begin(), late_extension_candidates_v.end() );      

      std::cout << "Plane " << p << " Early Extension matching" << std::endl;
      for (int iex=0; iex<(int)early_extension_candidates_v.size(); iex++) {
	auto& candidate = early_extension_candidates_v[iex];
	std::cout << " [" << candidate.idx << "] dist=" << candidate.dist << " cos=" << candidate.cosine << " triscore=" << candidate.triscore << std::endl;
	if ( candidate.cosine<-0.7 && candidate.dist<200.0 && (candidate.triscore<-0.5 || candidate.dist<10) ) {
	  std::cout << "  Match: " << candidate.cosine << " " << candidate.triscore << " " << candidate.dist << std::endl;
	  // append best to early match
	  auto const& early_contour = plane_contours_v[p][candidate.idx];
	  cluster.earlyContours[p].push_back( &early_contour );
	  cluster.earlyDir[p].push_back( early_contour.getStartDir() );
	  cluster.earlyEnd[p].push_back( early_contour.getFitSegmentStart() );
	  cluster.indices[p].insert( candidate.idx );
	  cluster_extended = true;
	  break;
	}
      }

      std::cout << "Plane " << p << " Late Extension matching" << std::endl;
      for (int iex=0; iex<(int)late_extension_candidates_v.size(); iex++) {
	auto& candidate = late_extension_candidates_v[iex];

	if ( candidate.cosine<-0.7 && candidate.dist<200.0 && (candidate.triscore<-0.5 || candidate.dist<10) ) {
	  std::cout << " [" << candidate.idx << "] dist=" << candidate.dist << " cos=" << candidate.cosine << std::endl;
	  // append best late match
	  auto& best_late = late_extension_candidates_v.front();
	  auto const& late_contour = plane_contours_v[p][best_late.idx];
	  cluster.lateContours[p] = &late_contour;
	  cluster.lateDir[p] = late_contour.getEndDir();
	  cluster.lateEnd[p] = late_contour.getFitSegmentEnd();
	  cluster.indices[p].insert( best_late.idx );
	  cluster_extended = true;
	  break;
	}
      }
      

      
    }//end of plane

    return cluster_extended;
  }//end of extension method


  bool BMTContourFilterAlgo::analyzeSeedContours( const std::vector<cv::Point>& imgpt, const std::vector<larcv::Image2D>& img_v,
						  const float max_dist2edge,
						  const ContourCluster& seedcluster ) {
    // for pairs of contours, we check that the edges of the contours in time are inside the image
    // if we have 3 planes, we check 3-plane consistency

    min_otherplane_v.clear();
    max_otherplane_v.clear(); 
    max_poszy_v.clear();
    min_poszy_v.clear();
    num_valid_pos = 0;
    num_combos = 0;

    int nplanes_w_contours = 0;
    for (int p1=0; p1<(int)seedcluster.size(); p1++) {
      if ( seedcluster[p1].size()==0 )
	continue;
      nplanes_w_contours++;
    }
    if ( nplanes_w_contours<2 )
      return false;
    
    num_combos = 1;
    for (int x=nplanes_w_contours; x>1; x--)
      num_combos *=x;
    num_combos /= 2;

    std::cout << __FILE__ << ":" << __FUNCTION__ << "--------------------------------" << std::endl;
    std::cout << " nplanes_w_contours=" << nplanes_w_contours << std::endl;
    std::cout << " number of combos=" << num_combos << std::endl;
    
    // We assume the seed cluster, so there should only be 1 or 0 contours per plane
    for (int p1=0; p1<(int)seedcluster.size(); p1++) {
      if ( seedcluster[p1].size()==0 )
	continue;

      for (int p2=p1+1; p2<(int)seedcluster.size(); p2++) {
	if ( seedcluster[p2].size()==0 )
	  continue;

	const ContourShapeMeta* contours[2];
	contours[0] = &(seedcluster[p1].front());
	contours[1] = &(seedcluster[p2].front());

	if ( contours[0]->getEndDir()[1]==0 || contours[1]->getEndDir()[1]==0 ) {
	  // not checking. so assume valid...
	  num_valid_pos++; 
	  continue;
	}

	// get min and max overlaps
	float max_common_row = -1;
	float min_common_row = -1;

	for (int i=0; i<2; i++) {
	  if ( max_common_row<0 || contours[i]->getFitSegmentStart().y > max_common_row ) {
	    max_common_row = contours[i]->getFitSegmentStart().y;	    
	  }
	  if ( min_common_row<0 || contours[i]->getFitSegmentEnd().y < min_common_row ) {
	    min_common_row = contours[i]->getFitSegmentEnd().y;
	  }

	  if ( min_common_row<0 || contours[i]->getFitSegmentStart().y < min_common_row ) {
	    min_common_row = contours[i]->getFitSegmentStart().y;	    
	  }
	  if ( max_common_row<0 || contours[i]->getFitSegmentEnd().y > max_common_row ) {
	    max_common_row = contours[i]->getFitSegmentEnd().y;
	  }
	  
	}

	// get interpolated positions
	std::vector< float > min_cols(2,0);
	for (int i=0; i<2; i++) {
	  std::vector<float> dir = contours[i]->getEndDir();
	  float dt = min_common_row - contours[i]->getFitSegmentStart().y;
	  min_cols[i]  = contours[i]->getFitSegmentStart().x + (dir[0]/dir[1])*dt;
	}

	std::vector< float > max_cols(2,0);
	for (int i=0; i<2; i++) {
	  std::vector<float> dir = contours[i]->getEndDir();
	  float dt = max_common_row - contours[i]->getFitSegmentStart().y;
	  max_cols[i]  = contours[i]->getFitSegmentStart().x + (dir[0]/dir[1])*dt;
	}
	
	// get intersection: minimum
	int min_otherplane = -1;
	int min_otherwire = -1;
	std::vector<float> min_intersection_zy;
	int min_crosses = 0;
	larcv::UBWireTool::getMissingWireAndPlane( p1, (int)min_cols[0], p2, (int)min_cols[1], min_otherplane, min_otherwire, min_intersection_zy, min_crosses ); 

	std::cout << "  "
		  << "MIN row=" << min_common_row << " "
		  << "seedcluster check planes=(" << p1 << "," << p2 << ") "
		  << "(z,y)=(" << min_intersection_zy[0] << "," << min_intersection_zy[1] << ") "
		  << "crosses=" << min_crosses << " "
		  << std::endl;

	// get intersection: maximum
	int max_otherplane = -1;
	int max_otherwire = -1;
	std::vector<float> max_intersection_zy;
	int max_crosses = 0;
	larcv::UBWireTool::getMissingWireAndPlane( p1, (int)max_cols[0], p2, (int)max_cols[1], max_otherplane, max_otherwire, max_intersection_zy, max_crosses ); 

	std::cout << "  "
		  << "MAX row=" << max_common_row << " "
		  << "seedcluster check planes=(" << p1 << "," << p2 << ") "
		  << "(z,y)=(" << max_intersection_zy[0] << "," << max_intersection_zy[1] << ") "
		  << "crosses=" << max_crosses << " "
		  << std::endl;
	
	// check positions
	bool max_valid = true;
	if ( max_intersection_zy[0]<-max_dist2edge || max_intersection_zy[0]>1036+max_dist2edge )
	  max_valid = false;
	if ( max_intersection_zy[1]<-117-max_dist2edge || max_intersection_zy[1]>117+max_dist2edge )
	  max_valid = false;
	if (max_valid) {
	  max_poszy_v.emplace_back( max_intersection_zy );
	  max_otherplane_v.push_back( max_otherplane );
	}
	
	bool min_valid = true;
	if ( min_intersection_zy[0]<-max_dist2edge || min_intersection_zy[0]>1036+max_dist2edge )
	  min_valid = false;
	if ( min_intersection_zy[1]<-117-max_dist2edge || min_intersection_zy[1]>117+max_dist2edge )
	  min_valid = false;
	if (min_valid) {
	  min_poszy_v.emplace_back( min_intersection_zy );
	  min_otherplane_v.push_back( min_otherplane );	  
	}

	if ( max_valid && min_valid) {
	  num_valid_pos++;
	}
	
      }//end of p2 loop
    }//end of p1 loop

    std::cout << "  Num valid planes: " << num_valid_pos << "/" << num_combos << std::endl;

    if ( num_valid_pos < num_combos )
      return false;
    
    return true;
  }
}
