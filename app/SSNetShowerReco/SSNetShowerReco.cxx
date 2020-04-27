#include "SSNetShowerReco.h"

#include "TMath.h"

#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPGraph.h"
#include "DataFormat/track.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {
namespace ssnetshowerreco {

  /**
   * constructor.
   *
   * sets some defaults for parameters
   */
  SSNetShowerReco::SSNetShowerReco(bool make_larlite, std::string outputname ) {

    _adc_tree_name = "wire";
    _calib_adc_tree_name = "calibrated";
    _ssnet_shower_image_stem = "ubspurn"; // sometimes ubspurn (when files made at FNAL)
    _vertex_tree_name = "test";
    _track_tree_name  = "trackReco";
    _partroi_tree_name = "segment";
    _mcshower_tree_name = "mcreco";
    _segment_tree_name = "segment";
    _instance_tree_name = "instance";
    _mctruth_name = "generator";
    _thrumu_tree_name = "wire";
    _ssnet_shower_tree_name = "ssnetshowerreco";

    _make_larlite = make_larlite;
    _Qcut = 10;
    _SSNET_SHOWER_THRESHOLD = 0.5;
    _use_calibrated_pixelsum2mev = false;
    _use_ncpi0 = false;
    _use_nueint = false;
    _use_bnb = false;
    _second_shower = false;
    _second_shower_adc_threshold = 5000;
    clear();
  }

  void SSNetShowerReco::initialize() {
    clear();
    setupAnaTree();
  }

  //set up output ana tree
  void SSNetShowerReco::setupAnaTree(){
    // OutFile->cd();
    std::string anatreename = _ssnet_shower_tree_name + "_anatree";
    _ana_tree = new TTree(anatreename.c_str(), "SSNet-based shower reco, analysis variables");
    _ana_tree->Branch("Run",&_run);
    _ana_tree->Branch("Subrun",&_subrun);
    _ana_tree->Branch("Event",&_event);
    _ana_tree->Branch("Vtxid",&_vtxid);
    _ana_tree->Branch("3dOpeningAngle",&_alpha);
    _ana_tree->Branch("Pi0Mass",&_pi0mass);
    _ana_tree->Branch("ImpactParameter1",&_impact1);
    _ana_tree->Branch("ImpactParameter2",&_impact2);
    _ana_tree->Branch("3dDirectionShower1",&_firstdirection);
    _ana_tree->Branch("3dDirectionShower2",&_seconddirection);
    _ana_tree->Branch("OverlapFractionY1",&_match_y1_vv);
    _ana_tree->Branch("OverlapFractionY2",&_match_y2_vv);
    _ana_tree->Branch("GapDistance1",&_shower_gap_vv);
    _ana_tree->Branch("GapDistance2",&_secondshower_gap_vv);
    _ana_tree->Branch("2DShowerStart1",&_shower_start_2d_vvv);
    _ana_tree->Branch("2DShowerStart2",&_secondshower_start_2d_vvv);
    _ana_tree->Branch("ShowerSumADC",&_shower_sumQ_vv);
    _ana_tree->Branch("ShowerLength",&_shower_shlength_vv);
    _ana_tree->Branch("ShowerEnergy",&_shower_energy_vv);
    _ana_tree->Branch("2dShowerDirection",&_shower_shangle_vv);
    _ana_tree->Branch("2dShowerOpeningAngle",&_shower_shopen_vv);
    _ana_tree->Branch("SecondShowerSumADC",&_secondshower_sumQ_vv);
    _ana_tree->Branch("SecondShowerLength",&_secondshower_shlength_vv);
    _ana_tree->Branch("SecondShowerEnergy",&_secondshower_energy_vv);
    _ana_tree->Branch("2dSecondShowerDirection",&_secondshower_shangle_vv);
    _ana_tree->Branch("2dSecondShowerOpeningAngle",&_secondshower_shopen_vv);
    _ana_tree->Branch("RecoVtxPos",&_vtx_pos_vv);
    _ana_tree->Branch("SmallQ",  &_smallQ);
    _ana_tree->Branch("SmallQ2", &_smallQ2);
    //trur variables for testing
    _ana_tree->Branch("truefid",&_truefid);
    _ana_tree->Branch("TrueVtxPos",&_true_vtx_3d_v);
    _ana_tree->Branch("TrueHasPi0",&_haspi0);
    _ana_tree->Branch("TrueCCNC",&_ccnc);
    _ana_tree->Branch("True3dDirectionShower2",&_seconddirection_true);
    _ana_tree->Branch("True3dDirectionShower1",&_firstdirection_true);
    _ana_tree->Branch("True2DShowerStart1",&_shower_start_2d_true_vvv);
    _ana_tree->Branch("True2DShowerStart2",&_secondshower_start_2d_true_vvv);
    _ana_tree->Branch("TrueShowerEnergy",&_shower_energy_true_vv);
    _ana_tree->Branch("TrueSecondShowerEnergy",&_secondshower_energy_true_vv);
    _ana_tree->Branch("TrueShowerRecoTrueDist",&_shower_recotrue_dist_v);
    _ana_tree->Branch("TrueSecondShowerRecoTrueDist",&_secondshower_recotrue_dist_v);

  }//end of setup ana tree funtion

  //write to outputfile
  void SSNetShowerReco::finalize(){
    // OutFile = TFile::Open("output_larlite.root","WRITE");
    // setupAnaTree();
    _ana_tree->Write();
  }//end of finalize function

  /**
   * clear result containers
   *
   */
  void SSNetShowerReco::clear() {
    _shower_energy_vv.clear();
    _shower_sumQ_vv.clear();
    _shower_shlength_vv.clear();
    _shower_shangle_vv.clear();
    _shower_shopen_vv.clear();
    _shower_gap_vv.clear();
    _shower_start_2d_vvv.clear();
    _secondshower_energy_vv.clear();
    _secondshower_sumQ_vv.clear();
    _secondshower_shlength_vv.clear();
    _secondshower_shangle_vv.clear();
    _secondshower_shopen_vv.clear();
    _secondshower_gap_vv.clear();
    _secondshower_start_2d_vvv.clear();
    _vtx_pos_vv.clear();
    _shower_ll_v.clear();
    _secondshower_ll_v.clear();
    _shower_pixcluster_v.clear();
    _true_energy_vv.clear();
    _true_shower_start_vv.clear();
    _true_shower_dir_vv.clear();
    _match_y1_vv.clear();
    _match_y2_vv.clear();
    _bestmatch_y1_vv.clear();
    _bestmatch_y2_vv.clear();
    _pi0mass.clear();
    _useformass.clear();
    _disttoint.clear();
    _impact1.clear();
    _impact2.clear();
    _firstdirection.clear();
    _seconddirection.clear();
    _alpha.clear();
    _vtxid.clear();
    _true_vtx_3d_v.clear();
    _seconddirection_true.clear();
    _firstdirection_true.clear();
    _shower_start_2d_true_vvv.clear();
    _secondshower_start_2d_true_vvv.clear();
    _shower_energy_true_vv.clear();
    _secondshower_energy_true_vv.clear();
    _shower_recotrue_dist_v.clear();
    _secondshower_recotrue_dist_v.clear();
    _smallQ.clear();
    _smallQ2.clear();
  }

  /**
   * use triple product to get area of triangle defined by (x_i,y_i)
   *
   */
  float SSNetShowerReco::_area( float x1, float y1,
                                float x2, float y2,
                                float x3, float y3 ) {
    return std::fabs((x1 * (y2 - y3) + x2 * (y3 - y1)  + x3 * (y1 - y2)) / 2.0);
  }

  /**
   * check is inside the triangle, using sum of partitions
   *
   * this was original test. deprecated.
   */
  bool SSNetShowerReco::_isInside(float x1, float y1,
                                  float x2, float y2,
                                  float x3, float y3,
                                  float x, float y ) {

    float A  = _area (x1 , y1 , x2 , y2 , x3 , y3);
    float A1 = _area (x  , y  , x2 , y2 , x3 , y3);
    float A2 = _area (x1 , y1 , x  , y  , x3 , y3);
    float A3 = _area (x1 , y1 , x2 , y2 , x  , y );

    if (std::fabs(A - (A1+A2+A3)) < 1e-6)
      return true;
    else
      return false;
  }

  /**
   * is test point (x3,y3) above(+),below(-1) the line defined by (x1,y1) and (x2,y2)
   */
  float SSNetShowerReco::_sign( float x1, float y1,
                                float x2, float y2,
                                float x3, float y3 ) {
    // (x1,y1) is the test point
    return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
  }

  /**
   * check if test point (x,y) is inside triangel defined by (x_i,y_i)
   *
   * verified to behave monotonically in contrast to _isInside1.
   * this is the default triangle containment test.
   *
   * @param[in] x test point x
   * @param[in] y test point y
   */
  bool SSNetShowerReco::_isInside2( float x1, float y1,
                                    float x2, float y2,
                                    float x3, float y3,
                                    float x, float y ) {

    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = _sign( x, y, x1, y1, x2, y2 );
    d2 = _sign( x, y, x2, y2, x3, y3 );
    d3 = _sign( x, y, x3, y3, x1, y1 );

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);

  }
  
  /**
   * get the enclosed charge in the whole image.
   *
   * @param[in]  chargeMap  Image with charge info
   * @param[in]  meta       Whole image meta
   * @param[out] cropCharge returned summed charge in crop
   * @param[in]  p          index of plane  
   * @param[in]  vtx_col    center column pixel (whole image coordinates)
   * @param[in]  vtx_row    center row pixel (whole image coordinates)
   *
   */
  void  SSNetShowerReco::_imageCharge( std::vector<std::vector<float>> img,
				       larcv::ImageMeta meta,
				       float& cropCharge, int p,
				       int vtx_col,
				       int vtx_row) {
	
    int tbounds[2] = { (int)(vtx_row-256*meta.pixel_height()),
		       (int)(vtx_row+256*meta.pixel_height()) };
    int wbounds[2] = { (int)(vtx_col-256*meta.pixel_width()),
		       (int)(vtx_col+256*meta.pixel_width()) };

    // correct bounds
    int wireMax = 2400;
    if(p == 2) wireMax = meta.max_x();
    
    if ( tbounds[0]<=meta.min_y() ) {
      tbounds[0] = meta.min_y()+meta.pixel_height();
      tbounds[1] = tbounds[0] + 512*meta.pixel_height();
    }
    if ( tbounds[1]>=meta.max_y() ) {
      tbounds[1] = meta.max_y() - meta.pixel_height();
      tbounds[0] = tbounds[1] - 512*meta.pixel_height();
    }
    if ( wbounds[0]<=meta.min_x() ) {
      wbounds[0] = meta.min_x()+meta.pixel_width();
      wbounds[1] = wbounds[0] + 512*meta.pixel_width();
    }
    if ( wbounds[1]>=wireMax ) {
      wbounds[1] = wireMax-meta.pixel_width();
      wbounds[0] = wbounds[1] - 512*meta.pixel_width();
    }
    
    std::cout << "tbounds: " << meta.row(tbounds[0]) << ", " << meta.row(tbounds[1]) << std::endl;
    std::cout << "wbounds: " << meta.col(wbounds[0]) << ", " << meta.col(wbounds[1]) << std::endl;
    std::cout << "First Row: " << img[0][0] << std::endl;
    std::cout << "First Col: " << img[0][1] << std::endl;
    std::cout << "Last Row: " << img[img.size()-1][0] << std::endl;
    std::cout << "Last Col: " << img[img.size()-1][1] << std::endl;
    float cropSum = 0;
    for (int ii =0;ii<(int)img.size();ii++){
      if(img[ii][0]>=meta.row(tbounds[1])&&img[ii][0]<meta.row(tbounds[0])&&img[ii][1]>=meta.col(wbounds[0])&&img[ii][1]<meta.col(wbounds[1])){
	cropSum+=img[ii][2];
      }
    }
    std::cout << "Crop sum: " << cropSum << std::endl;
    cropCharge = cropSum;

 } // end of imageCharge


  /**
   * get the enclosed charge inside a triangle.
   *
   * @param[in]  chargeMap Image with charge info
   * @param[in]  theta     Angle in radians in (col,row) coordinates that defines shower 2D direction
   * @param[out] sumIn     total pixel sum of pixels inside triangle
   * @param[out] triangle  set of 3 points that defines triangle used
   * @param[in]  vtx_col   start point of triangle
   * @param[in]  vtx_row   start point of triangle
   * @param[in]  shLen     length of triangle to use, radial line from vertex
   * @param[in]  shOpen    shower opening angle
   *
   */
  std::vector< std::vector<int> >  SSNetShowerReco::_enclosedCharge( std::vector<std::vector<float>> chargeMap,
                                                                     float theta,
                                                                     float& sumIn,
                                                                     std::vector< std::vector<float> >& triangle,
                                                                     bool calcE,
                                                                     int vtx_col,
                                                                     int vtx_row,
                                                                     float shLen,
                                                                     float shOpen) {

    std::vector< std::vector<int> > pix_v;

    sumIn = 0.;
    triangle.resize(3);
    for (int i=0; i<3; i++ ) {
      triangle[i].resize(2,0);
    }

    std::vector<int> vtx = { vtx_col, vtx_row };
    float t1X = vtx[0];
    float t1Y = vtx[1];

    float t2X = vtx[0] + shLen*cos(theta+shOpen);
    float t2Y = vtx[1] + shLen*sin(theta+shOpen);

    float t3X = vtx[0] + shLen*cos(theta-shOpen);
    float t3Y = vtx[1] + shLen*sin(theta-shOpen);

    triangle[0] = { t1X, t1Y };
    triangle[1] = { t2X, t2Y };
    triangle[2] = { t3X, t3Y };

    if ( shLen==0.0 ) {
      // weird things when we have zero area triangle
      bool hasCharge = false;
      // std::cout<<"------"<<vtx_row<<","<<vtx_col<<"------------"<<std::endl;
      for (int ii =0;ii<(int)chargeMap.size();ii++){
        // std::cout<<chargeMap[ii][0]<<","<<chargeMap[ii][1]<<std::endl;
        if(chargeMap[ii][0]==vtx_row&&chargeMap[ii][1]==vtx_col){
          sumIn = chargeMap[ii][2];
          hasCharge =true;
        }
      }
      //default to fall back on
      if (!hasCharge) sumIn = 10.0;
      return pix_v;
    }

    sumIn = 0;
    for (int ii =0;ii<(int)chargeMap.size();ii++){
      if ( _isInside2(t1X,t1Y,t2X,t2Y,t3X,t3Y,chargeMap[ii][1],chargeMap[ii][0])){
        if (calcE) sumIn+=chargeMap[ii][2];
        else sumIn+=1;
        std::vector<int> pix = { (int)chargeMap[ii][0], (int)chargeMap[ii][1] };
        pix_v.push_back( pix );
      }
    }
    return pix_v;
  }

  /*
   * vary direction (viz. angle) of shower at vertex to collect maximum charge
   *
   * @param[in]  chargeMap  Image with charge pixels
   * @param[in]  vtx_col    vertex in image
   * @param[in]  vtx_row    vertex in image
   * @param[in]  scanLen    fixed triangle length used to scan
   * @param[in]  scanOpen   fixed triangle opening angle used to scan
   * @return                best angle found; -9999 if no angle found charge
   */
  float SSNetShowerReco::_findDir( std::vector<std::vector<float>> chargeMap,
                                   int vtx_col,
                                   int vtx_row,
                                   float scanLen,
                                   float scanOpen ) {

    float bestDir   = -9999;
    float maxCharge = -9999;
    std::vector< std::vector<float> > bestTri;

    // coarse steps
    int coarseSteps = 72;
    std::vector<float> coarseAngs( coarseSteps, 0 );
    for ( int i=0; i<coarseSteps; i++ ) {

      coarseAngs[i]  = 2*TMath::Pi() / (float)coarseSteps * (float)i;
      float ang = coarseAngs[i];
      float sumQ;
      std::vector< std::vector<float> > triangle;
      _enclosedCharge( chargeMap, ang, sumQ, triangle, false,
                       vtx_col,vtx_row,scanLen,scanOpen);
      // std::cout << " find-dir ang=" << ang*180.0/TMath::Pi() << " deg; =" << ang << " rad;  sumQ=" << sumQ << std::endl;
      if ( sumQ > maxCharge ) {
        maxCharge = sumQ;
        bestDir   = ang;
        bestTri   = triangle;
      }
    }

    // Fine Tune
    float fineBestDir = -9999;
    maxCharge = -9999;

    for (int i=-5; i<=5; i++ ) {
      //angs = [bestDir + pi/180*i for i in range(-5,5)]
      float ang = bestDir + TMath::Pi()/180.0*(float)i;
      float sumQ;
      std::vector< std::vector<float> > triangle;
      _enclosedCharge( chargeMap, ang, sumQ, triangle, false,
                       vtx_col, vtx_row, scanLen, scanOpen);
      if ( sumQ > maxCharge ) {
        maxCharge = sumQ;
        fineBestDir   = ang;
        bestTri   = triangle;
      }
    }

    return fineBestDir;
  }


  /*
   * vary length of shower at vertex to collect maximum charge
   *
   * @param[in]  chargeMap  Image with charge pixels
   * @param[in]  theta      fixed triangle direction used in scan
   * @param[in]  vtx_col    vertex in image
   * @param[in]  vtx_row    vertex in image
   * @param[in]  scanOpen   fixed triangle opening angle used to scan
   * @return                best angle found
   */
  float SSNetShowerReco::_findLen( std::vector<std::vector<float>> chargeMap,
                                   float theta,
                                   int vtx_col, int vtx_row,
                                   float scanOpen ) {

    const int step     = 10;
    std::vector<float> lengths( (int)400/step, 0.0);
    for (int i=0; i<(int)400/step; i++ ) {
      lengths[i]  = step*i;
    }
    std::vector<float> charge;
    charge.reserve( lengths.size() );

    for ( auto const& length : lengths ) {
      float sumQ;
      triangle_t triangle;
      _enclosedCharge(chargeMap,theta,
                      sumQ, triangle, false,
                      vtx_col, vtx_row,
                      length,scanOpen);

      // std::cout << " find-len " << length << "; "
      //           << " area=" << _area( triangle[0][0], triangle[0][1],
      //                                 triangle[1][0], triangle[1][1],
      //                                 triangle[2][0], triangle[2][1] )
      //           << " sumQ=" << sumQ << std::endl;

      charge.push_back( sumQ );
    }

    float maxCharge = *std::max_element( charge.begin(), charge.end() );
    //std::cout << "    maxcharge=" << maxCharge << std::endl;
    int shStop = 0;
    for (size_t i=0; i<charge.size(); i++ ) {
      float x = charge[i];
      if ( x > 0.95*maxCharge ) {
        shStop   = i;
        break;
      }
    }

    return (shStop+1)*step;
  }

  /*
   * vary opening angle of shower at vertex to collect maximum charge
   *
   * @param[in]  chargeMap  Image with charge pixels
   * @param[in]  theta      fixed triangle direction used in scan
   * @param[in]  length     fixed triangle length used to scan
   * @param[in]  vtx_col    vertex in image
   * @param[in]  vtx_row    vertex in image
   * @return                best angle found
   */
  float SSNetShowerReco::_findOpen( std::vector<std::vector<float>> chargeMap,
                                    float theta,
                                    float length,
                                    int vtx_col, int vtx_row ) {

    float step    = 0.02;
    std::vector<float> opens( 14, 0 );
    for (int i=1; i<15; i++ ) {
      opens[i-1] = step*i;
    }
    std::vector<float> charge;
    charge.reserve( opens.size() );

    for (auto const& op : opens ) {
      float sumQ;
      std::vector< std::vector<float> > triangle;
      _enclosedCharge(chargeMap,theta,sumQ,triangle, false,
                      vtx_col, vtx_row, length, op);
      charge.push_back(sumQ);
    }

    float maxCharge = *std::max_element(charge.begin(), charge.end());
    int opStop = 0;
    for (size_t i=0; i<charge.size(); i++ ) {
      float x  = charge[i];
      if (x > 0.98*maxCharge) {
        opStop   = i;
        break;
      }
    }

    return (float)(opStop+1)*step;
  }


  /*
   * process one event using data from larcv and larlite IO managers
   *
   * @param[in]  iolcv  larcv input:   needs to contain ssnet images and pgraph with vertex information
   * @param[in]  ioll   larlite input: not used
   * @param[in]  ioimgs larcv input:   needs to contain charge images
   * @return            true if processing ok, false if not
   */
   std::vector<std::vector<float>> SSNetShowerReco::MakeImage2dSparse(const larcv::Image2D& input_img,
								      float threshold ) {
     //function to make sparse objects from the larcv Image2d
     //output: vector of (row,col, value)
     std::vector<std::vector<float>> output_vv;
     for (int r = 0; r<(int)input_img.meta().rows();r++){
       for (int c= 0; c<(int)input_img.meta().cols();c++){
         float pix_value = input_img.pixel(r,c);
         if (pix_value > threshold){
           std::vector<float> tmp_v = {(float) r, (float) c, pix_value};
           output_vv.push_back(tmp_v);
         }
       }//end of col loop
     }//end of row loop
     return output_vv;
   }// end of function

  //------------------------------------------------------------------------------  
  /*
   * sparsify given image using one threshold 
   *   inside some distance from radius and another threshold outside
   *
   * @param[in]  input_img   larcv input image
   * @param[in]  vtx_pix_row Image row where vertex is
   * @param[in]  vtx_pix_col Image col where vertex is
   * @param[in]  radius      Radius within which we use the special threshold
   * @param[in]  threshold_inside  Threshold to use inside the vertex
   * @param[in]  threshold_outside Threshold to use outside the vertex
   * @return     Sparse image, a vector of (row,column,pixel value)
   */
   std::vector<std::vector<float>>
   SSNetShowerReco::MakeImage2dSparseWithVertexThreshold(const larcv::Image2D& input_img,
							 const std::vector< std::vector<int> >& imgcoord_vv,
							 int plane,
							 float radius,
							 float threshold_inside,
							 float threshold_outside )
   {
     //function to make sparse objects from the larcv Image2d
     //output: vector of (row,col, value)
     std::vector<std::vector<float>> output_vv;     
     int ninside=0;
     for (int r = 0; r<(int)input_img.meta().rows();r++){
       for (int c= 0; c<(int)input_img.meta().cols();c++){
         float pix_value = input_img.pixel(r,c);
	 bool near_vtx = false;
	 for ( auto const& imgcoord_v : imgcoord_vv ) {
	   int vtx_pix_row = imgcoord_v[3];
	   int vtx_pix_col = imgcoord_v[plane];
	   float pix_rad = sqrt( (vtx_pix_row-r)*(vtx_pix_row-r) + (vtx_pix_col-c)*(vtx_pix_col-c) );
	   if ( pix_rad<=radius ) near_vtx = true;
	   if ( near_vtx ) break;
	 }

	 if ( near_vtx ) ninside++;
	 
	 if ( (!near_vtx && pix_value > threshold_outside )
	      || (near_vtx && pix_value > threshold_inside) ) {
	   std::vector<float> tmp_v = {(float) r, (float) c, pix_value};
	   output_vv.push_back(tmp_v);
	 }
       }//end of col loop
     }//end of row loop

     std::cout << "[ SSNetShowerReco::MakeImage2dSparseWithVertexThreshold ] "
	       << " nvertices=" << imgcoord_vv.size()
	       << " threshold inside=" << threshold_inside << " outside=" << threshold_outside
	       << " npix-inside=" << ninside++
	       << std::endl;
     
     return output_vv;
   }// end of function
  
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------

  bool SSNetShowerReco::process( larcv::IOManager& iolcv, larlite::storage_manager& ioll, int entry ) {//, larcv::IOManager& ioimgs ) {

    // get adc image (larcv)
    // get ssnet image (larcv)
    // get vertex (larcv)
    // get tracks
    // store in root ana tree, store in json file


    // clear result container
    clear();
    // OutFile->cd();
    //get second shower functions
    SecondShower SecondShower;
    Utils Utils;
    bool allowGap = true;
    bool makeDisp = false;
    bool useTrueVtx = false;

    std::vector<int> vtx2d_true;

    _run    = ioll.run_id();
    _subrun = ioll.subrun_id();
    _event  = ioll.event_id();
    std::cout<<"(r,s,e)="<<_run<<","<<_subrun<<","<<_event<<std::endl;
    //load in inputs

    larcv::EventImage2D* ev_adc
      = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _adc_tree_name );
    const std::vector<larcv::Image2D>& adc_v = ev_adc->Image2DArray();
    larcv::EventImage2D* ev_thrumu
      = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _thrumu_tree_name );
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();


    std::vector<std::vector<float>> adc_u_sparse = MakeImage2dSparse(adc_v[0],10.0);
    std::vector<std::vector<float>> adc_v_sparse = MakeImage2dSparse(adc_v[1],10.0);
    std::vector<std::vector<float>> adc_y_sparse = MakeImage2dSparse(adc_v[2],10.0);
    std::vector<std::vector<std::vector<float>>> adc_sparse_vvv = {adc_u_sparse,adc_v_sparse,adc_y_sparse};
    std::vector<larcv::ImageMeta> wire_meta = {adc_v[0].meta(),adc_v[1].meta(),adc_v[2].meta()};

    larcv::EventImage2D* ev_shower_score[3] = { nullptr };
    for ( size_t p=0; p<3; p++ ) {
      char treename[50];
      sprintf( treename, "%s_plane%d", _ssnet_shower_image_stem.c_str(), (int)p );
      ev_shower_score[p] =
        (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, treename );
    }

    // ===================================================================
    // GET VERTICES AND MC INFO
    // -------------------------

    larcv::EventPGraph* ev_vtx
      = (larcv::EventPGraph*)iolcv.get_data( larcv::kProductPGraph, _vertex_tree_name );
    larcv::EventROI* ev_partroi = nullptr;
    larlite::event_mcshower* ev_mcshower = nullptr;
    larcv::EventImage2D* ev_segment = nullptr;
    larcv::EventImage2D* ev_instance = nullptr;
    larlite::event_mctruth* ev_mctruth = nullptr;

    // Vertices

    if (_use_mc){
      ev_partroi = (larcv::EventROI*)(iolcv.get_data( larcv::kProductROI, _partroi_tree_name));
      ev_mcshower = ((larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,  _mcshower_tree_name));
    }

    if (_use_mc){
    // if (_use_nueint || _use_ncpi0){
      ev_segment = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _segment_tree_name );
      ev_instance = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _instance_tree_name );
    }
    if (_use_mc){
    // if (_use_bnb){
      ev_mctruth        = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,  _mctruth_name );
    }

    // get candidate vertices, if in fiducial volume, keep
    _vtx_pos_vv.clear();
    
    if (!useTrueVtx){
      int vtxid =0;
      for ( auto const& pgraph : ev_vtx->PGraphArray() ) {
        if ( pgraph.ParticleArray().size()==0 ) continue; // dont expect this
        auto const& roi = pgraph.ParticleArray().front();
        std::vector<double> vtx3d = { roi.X(), roi.Y(), roi.Z() };
        std::cout << "Load Vertex[" << vtxid << "] Pos (" << vtx3d[0] << "," << vtx3d[1] << "," << vtx3d[2] << ")" << std::endl;
        if ((vtx3d[0]>    0.001) && (vtx3d[0] <  255.999) && (vtx3d[1]> -116.499) && (vtx3d[1] < 116.499)
            && (vtx3d[2]>    0.001) && (vtx3d[2] < 1036.999)){

          _vtx_pos_vv.push_back( vtx3d );
          _vtxid.push_back(vtxid);
          // std::cout<<"vertex id: "<<vtxid<<std::endl;

        }//end of if statement
        else std::cout<<"Outside fiducial volume!"<<std::endl;
        vtxid++;
      }//end of loop over vertices/pgraph
    }
    else{
      // MCC9 SCE correction
      TFile* newSCEFile = new TFile("/cluster/tufts/wongjiradlab/rshara01/ubdl/SCEoffsets_dataDriven_combined_fwd_Jan18.root");
      TH3F* sceDx = (TH3F*) newSCEFile->Get("hDx");
      TH3F* sceDy = (TH3F*) newSCEFile->Get("hDy");
      TH3F* sceDz = (TH3F*) newSCEFile->Get("hDz");

      std::vector<double> _scex(3,0);
      std::vector<double> _tx(3,0);
      std::vector<double> vtx3d_true;
			

      // get true variables---------------------------------------------------
      for(auto const& roi : ev_partroi->ROIArray()){
        if(std::abs(roi.PdgCode()) == 12 || std::abs(roi.PdgCode()) == 14) {
          _tx.resize(3,0);
          _tx[0] = roi.X();
          _tx[1] = roi.Y();
          _tx[2] = roi.Z();
          auto const offset = Utils.GetPosOffsets(_tx,sceDx,sceDy,sceDz);
          _scex = Utils.MakeSCECorrection(_scex,_tx,offset);

        }
      }
      newSCEFile->Close();
      vtx3d_true = _scex;
      vtx2d_true = Utils.getProjectedPixel(_scex, ev_adc->Image2DArray()[2].meta(), 3);

    }

    // ---------------------------------------------------------
    // GET MC INFO
    //-------------remove later----------------------------------
    if (_use_mc){
      // ev_partroi = (larcv::EventROI*)(iolcv.get_data( larcv::kProductROI, _partroi_tree_name));
      // TFile* newSCEFile = new TFile("/cluster/tufts/wongjiradlab/rshara01/bkp/SCEoffsets_dataDriven_combined_fwd_Jan18.root","read");
      // TH3F* sceDx = (TH3F*) newSCEFile->Get("hDx");
      // TH3F* sceDy = (TH3F*) newSCEFile->Get("hDy");
      // TH3F* sceDz = (TH3F*) newSCEFile->Get("hDz");
      std::vector<double> _scex(3,0);
      std::vector<double> _tx(3,0);
      for(auto const& roi : ev_partroi->ROIArray()){
        if(std::abs(roi.PdgCode()) == 12 || std::abs(roi.PdgCode()) == 14) {
          _tx.resize(3,0);
          _tx[0] = roi.X();
          _tx[1] = roi.Y();
          _tx[2] = roi.Z();
          // auto const offset = Utils.GetPosOffsets(_tx,sceDx,sceDy,sceDz);
          // _scex = Utils.MakeSCECorrection(_scex,_tx,offset);
        }
      }
      std::cout<<"TRUE VTX: "<<_tx[0]<<" "<<_tx[1]<<" "<<_tx[2]<<std::endl;
      bool vtxinfid = Utils.InsideFiducial(_tx[0],_tx[1],_tx[2]);
      if (vtxinfid) _truefid =1;
      else _truefid =0;
      _true_vtx_3d_v = _tx;
      // newSCEFile->Close();

      if (vtxinfid ){
        _true_energy_vv = SecondShower.SaveTrueEnergies(ev_mcshower);
        _true_shower_start_vv = SecondShower.SaveTrueStarts(ev_mcshower);
        _true_shower_dir_vv = SecondShower.SaveTrueDirections(ev_mcshower);
      }

      _haspi0 = 0;
      _ccnc = 0;
      for(int part =0;part<(int)ev_mctruth->at(0).GetParticles().size();part++){

        if (ev_mctruth->at(0).GetParticles().at(part).PdgCode() == 111) _haspi0 = 1;
        if (ev_mctruth->at(0).GetNeutrino().CCNC() == 1) _ccnc = 1;

      }
      if (_haspi0 ==1) std::cout<<"HAS PI0!"<<std::endl;
      std::cout<<"here"<<std::endl;

    }// end of if MC

    // ===================================================================================
    // GET PROJECTED IMAGE COORD OF VERTICES
    // --------------------------------------
    
    std::vector< std::vector<int>  > imgcoord_vv;

    for ( size_t ivtx=0; ivtx<_vtx_pos_vv.size(); ivtx++ ) {
      
      std::cout << "[SSNetShowerReco] Reconstruct vertex[" << ivtx << "] "
                << "pos=(" << _vtx_pos_vv[ivtx][0] << "," << _vtx_pos_vv[ivtx][1] << "," << _vtx_pos_vv[ivtx][2] << ")"
                << std::endl;

      //get 2d vertex position
      std::vector<int> imgcoord_v(4);  // (U,V,Y,tick)
      float tick = 3200 + _vtx_pos_vv[ivtx][0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
      for ( size_t p=0; p<3; p++ ) {
        imgcoord_v[p] = adc_v[p].meta().col( larutil::Geometry::GetME()->NearestWire( _vtx_pos_vv[ivtx], (int)p ), __FILE__, __LINE__ );
      }
      imgcoord_v[3] = adc_v[0].meta().row( tick, __FILE__, __LINE__ );
      imgcoord_vv.push_back( imgcoord_v );
    }
    

    // ============================================
    // MAKE SPARSE IMAGES FOR SSNET SCORES
    // OLD FIXED THRESHOLD
    // std::vector<std::vector<float>> ssnet_u_shower_sparse = MakeImage2dSparse(ev_shower_score[0]->Image2DArray()[0],_SSNET_SHOWER_THRESHOLD);
    // std::vector<std::vector<float>> ssnet_v_shower_sparse = MakeImage2dSparse(ev_shower_score[1]->Image2DArray()[0],_SSNET_SHOWER_THRESHOLD);
    // std::vector<std::vector<float>> ssnet_y_shower_sparse = MakeImage2dSparse(ev_shower_score[2]->Image2DArray()[0],_SSNET_SHOWER_THRESHOLD);

    // TWO THRESHOLD METHOD, for pixels near vertices
    std::vector<std::vector<float>> ssnet_u_shower_sparse
      = MakeImage2dSparseWithVertexThreshold(ev_shower_score[0]->Image2DArray()[0],
					     imgcoord_vv,
					     0,
					     10.0,
					     0.05,
					     _SSNET_SHOWER_THRESHOLD);
    std::vector<std::vector<float>> ssnet_v_shower_sparse
      = MakeImage2dSparseWithVertexThreshold(ev_shower_score[1]->Image2DArray()[0],
					     imgcoord_vv,
					     1,
					     10.0,
					     0.05,
					     _SSNET_SHOWER_THRESHOLD);
    std::vector<std::vector<float>> ssnet_y_shower_sparse
      = MakeImage2dSparseWithVertexThreshold(ev_shower_score[2]->Image2DArray()[0],
					     imgcoord_vv,
					     2,
					     10.0,
					     0.05,
					     _SSNET_SHOWER_THRESHOLD);
    std::vector<std::vector<std::vector<float>>> ssnet_sparse_vvv = {ssnet_u_shower_sparse,ssnet_v_shower_sparse,ssnet_y_shower_sparse};
    
    std::vector<std::vector<float>> ssnet_u_shower_sparse_high = MakeImage2dSparse(ev_shower_score[0]->Image2DArray()[0],.5);
    std::vector<std::vector<float>> ssnet_v_shower_sparse_high = MakeImage2dSparse(ev_shower_score[1]->Image2DArray()[0],.5);
    std::vector<std::vector<float>> ssnet_y_shower_sparse_high = MakeImage2dSparse(ev_shower_score[2]->Image2DArray()[0],.5);
    std::vector<std::vector<std::vector<float>>> ssnet_sparse_high_vvv = {ssnet_u_shower_sparse,ssnet_v_shower_sparse,ssnet_y_shower_sparse};
    
    std::vector<std::vector<float>> ssnet_u_shower_sparse_low = MakeImage2dSparse(ev_shower_score[0]->Image2DArray()[0],.05);
    std::vector<std::vector<float>> ssnet_v_shower_sparse_low = MakeImage2dSparse(ev_shower_score[1]->Image2DArray()[0],.05);
    std::vector<std::vector<float>> ssnet_y_shower_sparse_low = MakeImage2dSparse(ev_shower_score[2]->Image2DArray()[0],.05);
    std::vector<std::vector<std::vector<float>>> ssnet_sparse_low_vvv = {ssnet_u_shower_sparse_low,ssnet_v_shower_sparse_low,ssnet_y_shower_sparse_low};


    //if in fid, save to branch
    //------------------------------------------------------------------------

    //mask out non shower pixels from adc
    std::vector<std::vector<std::vector<float>>> masked_adc_vvv;
    std::vector<std::vector<std::vector<float>>> masked_adc_high_vvv;
    std::vector<std::vector<std::vector<float>>> masked_adc_low_vvv;


    //make object to save triangles of showers to fill later.
    //for each u1,u2,v1,v2,y1,y2: save (x1,y1,x2,y2,x3,y3)
    std::vector<std::vector<float>> shower_points_vv(6,std::vector<float> (6,-1));


    for ( size_t p=0; p<3; p++ ) {

      std::vector<int> vtx2d_true;


      // mask the ADC image using ssnet - low threshold
      // Used for shower clustering
      int npixels_adc = adc_sparse_vvv[p].size();
      int npixels_ssnet = ssnet_sparse_vvv[p].size();
      std::vector<std::vector<float>> masked_plane_v;
      for ( int adc_i =0;adc_i<npixels_adc;adc_i++ ) {
        for ( int ssnet_i =0;ssnet_i<npixels_ssnet;ssnet_i++){
          float thrumupix = thrumu_v[p].pixel(adc_sparse_vvv[p][adc_i][0],adc_sparse_vvv[p][adc_i][1]);
          if (ssnet_sparse_vvv[p][ssnet_i][0] == adc_sparse_vvv[p][adc_i][0] && ssnet_sparse_vvv[p][ssnet_i][1] == adc_sparse_vvv[p][adc_i][1]){
              masked_plane_v.push_back({(float) adc_sparse_vvv[p][adc_i][0],adc_sparse_vvv[p][adc_i][1],adc_sparse_vvv[p][adc_i][2]});
          }//end of if statements
        }//end of loop through ssnet
      }//end of loop throug adc
      masked_adc_vvv.push_back(masked_plane_v);


      // mask the ADC image using ssnet - high threshold
      int npixels_ssnet_high = ssnet_sparse_high_vvv[p].size();
      std::vector<std::vector<float>> masked_plane_high_v;
      for ( int adc_i =0;adc_i<npixels_adc;adc_i++ ) {
        for ( int ssnet_i =0;ssnet_i<npixels_ssnet_high;ssnet_i++){
          if (ssnet_sparse_high_vvv[p][ssnet_i][0] == adc_sparse_vvv[p][adc_i][0] && ssnet_sparse_high_vvv[p][ssnet_i][1] == adc_sparse_vvv[p][adc_i][1]){
              masked_plane_high_v.push_back({(float) adc_sparse_vvv[p][adc_i][0],adc_sparse_vvv[p][adc_i][1],adc_sparse_vvv[p][adc_i][2]});
          }//end of if statements
        }//end of loop through ssnet
      }//end of loop throug adc
      masked_adc_high_vvv.push_back(masked_plane_high_v);
      
			
      // mask the ADC image using ssnet - high threshold
      int npixels_ssnet_low = ssnet_sparse_low_vvv[p].size();
      std::vector<std::vector<float>> masked_plane_low_v;
      for ( int adc_i =0;adc_i<npixels_adc;adc_i++ ) {
        for ( int ssnet_i =0;ssnet_i<npixels_ssnet_low;ssnet_i++){
          if (ssnet_sparse_low_vvv[p][ssnet_i][0] == adc_sparse_vvv[p][adc_i][0] && ssnet_sparse_low_vvv[p][ssnet_i][1] == adc_sparse_vvv[p][adc_i][1]){
              masked_plane_low_v.push_back({(float) adc_sparse_vvv[p][adc_i][0],adc_sparse_vvv[p][adc_i][1],adc_sparse_vvv[p][adc_i][2]});
          }//end of if statements
        }//end of loop through ssnet
      }//end of loop throug adc
      masked_adc_low_vvv.push_back(masked_plane_low_v);
    
    }//end of plane loop

    if (_use_nueint){
      _uplane_profile_vv = SecondShower.SaveTrueProfile(0, ev_segment, ev_instance,
            ev_mcshower,masked_adc_vvv);
      _vplane_profile_vv = SecondShower.SaveTrueProfile(1, ev_segment, ev_instance,
            ev_mcshower,masked_adc_vvv);
      _yplane_profile_vv = SecondShower.SaveTrueProfile(2, ev_segment, ev_instance,
            ev_mcshower,masked_adc_vvv);
    }
    // std::cout<<"Size of masked adc sparse vectors: \n";
    // std::cout<<"Masked ADC: "<<masked_adc_vvv[0].size()<<","<<masked_adc_vvv[1].size()<<","<<masked_adc_vvv[2].size()<<"\n";





    // now get shower energy per vertex, per plane

    for ( size_t ivtx=0; ivtx<_vtx_pos_vv.size(); ivtx++ ) {

      float pi0mass = 0;
      float disttoint = 0;
      float impact1 = 0;
      float impact2 = 0;
      std::vector<float> firstdirection (3,-1);
      std::vector<float> seconddirection (3,-1);
      float alpha;

      std::cout << "[SSNetShowerReco] Reconstruct vertex[" << ivtx << "] "
                << "pos=(" << _vtx_pos_vv[ivtx][0] << "," << _vtx_pos_vv[ivtx][1] << "," << _vtx_pos_vv[ivtx][2] << ")"
                << std::endl;

      //get 2d vertex position
      std::vector<int> imgcoord_v(4);  // (U,V,Y,tick)
      imgcoord_v[3] = 3200 + _vtx_pos_vv[ivtx][0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
      for ( size_t p=0; p<3; p++ ) {
        imgcoord_v[p] = larutil::Geometry::GetME()->NearestWire( _vtx_pos_vv[ivtx], (int)p );
      }
      std::vector<std::vector<int>> shower_start_2d_vv;
      std::vector<int>   shower_gap_v(3,0);
      std::vector<float> shower_energy_v(3,0);
      std::vector<float> shower_sumQ_v(3,0);
      std::vector<float> shower_shlength_v(3,0);
      std::vector<float> shower_shangle_v(3,0);
      std::vector<float> shower_shopen_v(3,0);
      std::vector<std::vector<int>> secondshower_start_2d_vv;
      std::vector<int>   secondshower_gap_v(3,0);
      std::vector<float> secondshower_energy_v(3,0);
      std::vector<float> secondshower_sumQ_v(3,0);
      std::vector<float> secondshower_shlength_v(3,0);
      std::vector<float> secondshower_shangle_v(3,0);
      std::vector<float> secondshower_shopen_v(3,0);
      std::vector<float> smallq_v(3,0);
      std::vector<float> smallq2_v(3,0);
				

      for ( size_t p=0; p<3; p++ ) {
        std::vector<int> vtx_pix2 = { (int)wire_meta[p].col(imgcoord_v[p]) , (int)wire_meta[p].row(imgcoord_v[3])};
        //std::cout << "[SSNetShowerReco]   Plane [" << p << "]" << std::endl;
        std::vector<int> vtx_pix = { (int)wire_meta[p].col(imgcoord_v[p]) , (int)wire_meta[p].row(imgcoord_v[3])};
        int vtx_pix_original[2] = { (int)wire_meta[p].col(imgcoord_v[p]) , (int)wire_meta[p].row(imgcoord_v[3])};
        //std::cout << "[SSNetShowerReco]     vertex pixel: (" << vtx_pix[0] << "," << vtx_pix[1] << ")" << std::endl;
        float shangle  = _findDir( masked_adc_vvv[p], vtx_pix[0], vtx_pix[1] );
        float shangle2 = -1;
        int step = -1;
        int step2 =-1;
        std::vector<int> shower_start_2d_v(2,-1);
        std::vector<int> secondshower_start_2d_v(2,-1);
        if (allowGap){
          int step = SecondShower.ChooseGapSize(vtx_pix[0],vtx_pix[1],shangle,
            masked_adc_vvv[p], 1);
          vtx_pix[0] = vtx_pix[0] + (step * cos(shangle));
          vtx_pix[1] = vtx_pix[1] + (step * sin(shangle));
        }

        // std::cout << "[SSNetShowerReco]     shower angle: " << shangle << std::endl;
        float shlength = _findLen( masked_adc_vvv[p],  shangle, vtx_pix[0], vtx_pix[1] );
        float shlength2 = -1;
        // std::cout << "[SSNetShowerReco]     shower length: " << shlength << std::endl;
        float shopen   = _findOpen( masked_adc_vvv[p], shangle, shlength, vtx_pix[0], vtx_pix[1] );
        float shopen2 =-1;
        // std::cout << "[SSNetShowerReco]     open angle: " << shopen << std::endl;
        float sumQ;
        float sumQ2 = -1 ;
        float smallQ;
        float smallQ2 = -1;
        triangle_t tri;
        std::vector< std::vector<int> > pixlist_v =
          _enclosedCharge( masked_adc_vvv[p], shangle, sumQ, tri, true, vtx_pix[0], vtx_pix[1], shlength, shopen );
	_imageCharge( masked_adc_low_vvv[p], wire_meta[p], smallQ, p, imgcoord_v[p], imgcoord_v[3]);
	smallq_v[p] = smallQ;

        float reco_energy = 0;
        float reco_energy2 = -1;
        if ( _use_calibrated_pixelsum2mev )
          reco_energy = sumQ*0.01324 + 37.83337; // calibrated images
        else
          reco_energy = sumQ*0.013456 + 2.06955; // uncalibrated images

        std::cout << "[SSNetShowerReco] plane[" << p << "] final sumQ=" << sumQ << " reco=" << reco_energy << std::endl;

        //run second shower code if requested
        if (_second_shower){
          //first mask out first shower
          //triangle coord: [x1,y1],[x2,y2],[x3,y3]
          std::vector<std::vector<double>> first_triangle;
          first_triangle.push_back({(double)vtx_pix[0],(double)vtx_pix[1]});
          first_triangle.push_back({vtx_pix[0] + shlength*cos(shangle+shopen),vtx_pix[1] + shlength*sin(shangle+shopen)});
          first_triangle.push_back({vtx_pix[0] + shlength*cos(shangle-shopen),vtx_pix[1] + shlength*sin(shangle-shopen)});


          std::vector<std::vector<float>> secondshower_adc_vv;
          //loop through planes
          float tot_in = 0;
          //loop through sparse image
          std::vector<std::vector<float>> tmp_masked_v;
          for (int ii = 0; ii<int(masked_adc_vvv[p].size());ii++){
            float adc_pix = masked_adc_vvv[p][ii][2];
            if(!_isInside2(first_triangle[0][0],first_triangle[0][1],
                first_triangle[1][0],first_triangle[1][1],first_triangle[2][0],
                first_triangle[2][1],masked_adc_vvv[p][ii][1],masked_adc_vvv[p][ii][0])){
              tot_in++;
              secondshower_adc_vv.push_back({masked_adc_vvv[p][ii][0],masked_adc_vvv[p][ii][1],adc_pix});
            }
          }//end of loop through sparse
          // std::cout<<"Size of new masked image "<<p<<": "<<secondshower_adc_vv.size()<<std::endl;

          //check remaining charge in image - have a threshold
          float remaining_adc =0;
          for (int ii=0;ii<(int)secondshower_adc_vv.size();ii++){
            remaining_adc+=secondshower_adc_vv[ii][2];
          }
          // save value
          _remaining_adc = remaining_adc;

          //only run second shower search if above adc threshold
          if (remaining_adc > _second_shower_adc_threshold){
            //run code again with mask allowed
            //reset to vtx coordinates

            shangle2  = _findDir( secondshower_adc_vv, vtx_pix2[0], vtx_pix2[1] );
            //always allow gap for second shower
            int step2 = SecondShower.ChooseGapSize(vtx_pix2[0],vtx_pix2[1],shangle2,
              secondshower_adc_vv, 2);
            vtx_pix2[0] = vtx_pix2[0] + (step2 * cos(shangle2));
            vtx_pix2[1] = vtx_pix2[1] + (step2 * sin(shangle2));

            // std::cout << "[SSNetShowerReco]    second shower angle: " << shangle2 << std::endl;
            shlength2 = _findLen( secondshower_adc_vv,  shangle2, vtx_pix2[0], vtx_pix2[1] );
            // std::cout << "[SSNetShowerReco]    second shower length: " << shlength2 << std::endl;
            shopen2   = _findOpen( secondshower_adc_vv, shangle2, shlength2, vtx_pix2[0], vtx_pix2[1] );
            // std::cout << "[SSNetShowerReco]    second open angle: " << shopen2 << std::endl;
            triangle_t tri2;
            std::vector< std::vector<int> > pixlist2_v =
              _enclosedCharge( secondshower_adc_vv, shangle2, sumQ2, tri2, true, vtx_pix2[0], vtx_pix2[1], shlength2, shopen2 );
            _imageCharge( masked_adc_low_vvv[p], wire_meta[p], smallQ2, p, vtx_pix[0], vtx_pix[1]);
	    smallq2_v[p] = smallQ2;

            if ( _use_calibrated_pixelsum2mev )
              reco_energy2 = sumQ2*0.01324 + 37.83337; // calibrated images
            else
              reco_energy2 = sumQ2*0.013456 + 2.06955; // uncalibrated images

            if (shangle2 < shangle+(shopen/2.0)&&shangle2 > shangle-(shopen/2.0) ){
	      std::cout<<"MATCHING ANGLE!" << std::endl;

	      std::vector< std::vector<double> > matching_triangle;
	      matching_triangle.push_back({(double)vtx_pix2[0],(double)vtx_pix2[1]});
	      matching_triangle.push_back({vtx_pix2[0] + shlength2*cos(shangle2+shopen2),vtx_pix2[1] + shlength2*sin(shangle2+shopen2)});
	      matching_triangle.push_back({vtx_pix2[0] + shlength2*cos(shangle2-shopen2),vtx_pix2[1] + shlength2*sin(shangle2-shopen2)});

	      
              // remask with this same direction fragment off
              secondshower_adc_vv.clear();
              //loop through sparse image
              for (int ii = 0; ii<int(masked_adc_vvv[p].size());ii++){
                float adc_pix = masked_adc_vvv[p][ii][2];
                if( !_isInside2(first_triangle[0][0],first_triangle[0][1],
				first_triangle[1][0],first_triangle[1][1],first_triangle[2][0],
				first_triangle[2][1],masked_adc_vvv[p][ii][1],masked_adc_vvv[p][ii][0]) &&
		    !_isInside2(matching_triangle[0][0],matching_triangle[0][1],
				matching_triangle[1][0],matching_triangle[1][1],matching_triangle[2][0],
				matching_triangle[2][1],masked_adc_vvv[p][ii][1],masked_adc_vvv[p][ii][0])		    
		    ){
                  secondshower_adc_vv.push_back({masked_adc_vvv[p][ii][0],masked_adc_vvv[p][ii][1],adc_pix});
                }
              }//end of loop through sparse
	      
	      if ( sumQ2>sumQ ) {
		std::cout << "SAME DIR SECOND SHOWER BIGGER THAN THE FIRST! REPLACE FIRST SHOWER" << std::endl;
		// bigger than first shower! modifying first shower"<<std::endl;
		//make 1st one == second
		shangle = shangle2;
		shopen = shopen2;
		shlength = shlength2;
		sumQ = sumQ2;
		tri = tri2;
		pixlist_v = pixlist2_v;
		reco_energy = reco_energy2;
		first_triangle.clear();
		first_triangle.push_back({(double)vtx_pix2[0],(double)vtx_pix2[1]});
		first_triangle.push_back({vtx_pix2[0] + shlength*cos(shangle+shopen),vtx_pix2[1] + shlength*sin(shangle+shopen)});
		first_triangle.push_back({vtx_pix2[0] + shlength*cos(shangle-shopen),vtx_pix2[1] + shlength*sin(shangle-shopen)});
	      }
	      else {
		// tmw: add to energy [maybe not great move]?
		sumQ += sumQ2;
		reco_energy += reco_energy2;
	      }

              //run code again with mask allowed
              //reset to vtx coordinates
              shangle2  = _findDir( secondshower_adc_vv, vtx_pix2[0], vtx_pix2[1] );
              //always allow gap for second shower
              step2 = SecondShower.ChooseGapSize(vtx_pix_original[0],vtx_pix_original[1],shangle2,
                secondshower_adc_vv, 2);
              vtx_pix2[0] = vtx_pix_original[0] + (step2 * cos(shangle2));
              vtx_pix2[1] = vtx_pix_original[1] + (step2 * sin(shangle2));

              // std::cout << "[SSNetShowerReco]    redo second shower angle: " << shangle2 << std::endl;
              shlength2 = _findLen( secondshower_adc_vv,  shangle2, vtx_pix2[0], vtx_pix2[1] );
              // std::cout << "[SSNetShowerReco]    second shower length: " << shlength2 << std::endl;
              shopen2   = _findOpen( secondshower_adc_vv, shangle2, shlength2, vtx_pix2[0], vtx_pix2[1] );
              // std::cout << "[SSNetShowerReco]    second open angle: " << shopen2 << std::endl;
              pixlist2_v =  _enclosedCharge( secondshower_adc_vv, shangle2, sumQ2, tri2, true, vtx_pix2[0], vtx_pix2[1], shlength2, shopen2 );
              reco_energy2 = 0;
              if ( _use_calibrated_pixelsum2mev )
                reco_energy2 = sumQ2*0.01324 + 37.83337; // calibrated images
              else
                reco_energy2 = sumQ2*0.013456 + 2.06955; // uncalibrated images

            } //end of overlap check
	    
            std::cout << "[SSNetShowerReco] plane[" << p << "] final sumQ2=" << sumQ2 << " reco2=" << reco_energy2 << std::endl;
	    
          }//end of remaining adc check
	  
          secondshower_gap_v[p] = step2;
          secondshower_start_2d_vv.push_back(vtx_pix2);
          secondshower_energy_v[p] = reco_energy2;
          secondshower_sumQ_v[p] = sumQ2;
          secondshower_shlength_v[p] = shlength2;
          secondshower_shangle_v[p]=shangle2;
          secondshower_shopen_v[p]=shopen2;
	  
        }//end of second shower search
        shower_gap_v[p] = step;
        shower_start_2d_vv.push_back(vtx_pix);
        shower_energy_v[p] = reco_energy;
        shower_sumQ_v[p] = sumQ;
        shower_shlength_v[p] = shlength;
        shower_shangle_v[p]=shangle;
        shower_shopen_v[p]=shopen;

        // make larlite
        // -------------
        larlite::shower shr;  // stores shower parameters
        larlite::shower shr2;  // stores second shower parameters
        larlite::larflowcluster pixcluster; // stores pixels

        // store vertex index
        shr.set_id( ivtx );
        shr2.set_id( ivtx );

        // store vertex
        TVector3 vtx3 = { _vtx_pos_vv[ivtx][0], _vtx_pos_vv[ivtx][1], _vtx_pos_vv[ivtx][2] };
        shr.set_start_point( vtx3 );
        shr2.set_start_point( vtx3 );

        // store plane
        shr.set_total_best_plane( p );
        shr2.set_total_best_plane( p );

        // set energy
        shr.set_total_energy( reco_energy );
        shr2.set_total_energy( reco_energy2 );
        // set sumq
        shr.set_total_energy_err( sumQ );
        shr2.set_total_energy_err( sumQ2 );
        
	// set cropSumq
        shr.set_dqdx( smallQ );
        shr2.set_dqdx( smallQ2 );

        // set length
        shr.set_length( shlength );
        shr2.set_length( shlength2 );

        // opening angle
        shr.set_opening_angle( shopen );
        shr2.set_opening_angle( shopen2 );

        // for direction store 3 vector
        // (x,y) are the wire and tick angle. (z) is the angle itself
        TVector3 pseudo_dir = { cos( shangle ), sin( shangle ), shangle };
        TVector3 pseudo_dir2 = { cos( shangle2 ), sin( shangle2 ), shangle2 };
        shr.set_direction( pseudo_dir );
        shr2.set_direction( pseudo_dir2 );

        _shower_ll_v.emplace_back( std::move(shr) );
        _secondshower_ll_v.emplace_back( std::move(shr2) );

        pixcluster.reserve( pixlist_v.size() );
        for ( auto& pix : pixlist_v ) {
          larlite::larflow3dhit lfpix;
          lfpix.resize(2,0);
          lfpix[0] = pix[0];
          lfpix[1] = pix[1];
          pixcluster.emplace_back( std::move(lfpix) );
        }
        _shower_pixcluster_v.emplace_back( std::move(pixcluster) );

        if (makeDisp&&p==2&& _truefid==1){

          TH2F* YPlaneADC_h = new TH2F("yplaneadc","yplaneadc",3456,0,3456.,1008,0,1008.);
          // TH2F* YPlaneSSNet_h = new TH2F("yplanessnet","yplanessnet",3456,0,3456.,1008,0,1008.);
          // TH2F* YPlaneSeg_h = new TH2F("yplaneseg","yplaneseg",3456,0,3456.,1008,0,1008.);
          // TH2F* YPlaneMasked_h = new TH2F("yplanemask","yplanemask",3456,0,3456.,1008,0,1008.);
          TH2F* vtx_h = new TH2F("vtx_h","vtx_h",3456,0,3456.,1008,0,1008.);
          TH2F* triangle1_h = new TH2F("triangle1","triangle1",3456,0,3456.,1008,0,1008.);
          TH2F* triangle2_h = new TH2F("triangle2","triangle2",3456,0,3456.,1008,0,1008.);
          TH2F* shr_h = new TH2F("shr_h","shr_h",3456,0,3456.,1008,0,1008.);
          //
          // for (int ii =0;ii<(int)masked_adc_vvv[p].size();ii++){
          //   YPlaneSeg_h->Fill(masked_adc_vvv[p][ii][1],masked_adc_vvv[p][ii][0],ev_segment->Image2DArray()[2].pixel(masked_adc_vvv[p][ii][0],masked_adc_vvv[p][ii][1]));
          // }
          std::cout<<"made hists"<<std::endl;
          for (int ii =0;ii<(int)adc_sparse_vvv[p].size();ii++){
            YPlaneADC_h->Fill(adc_sparse_vvv[p][ii][1],adc_sparse_vvv[p][ii][0],adc_sparse_vvv[p][ii][2]);
          }
          std::cout<<"FILLED ADC"<<std::endl;
          // for (int ii =0;ii<(int)masked_adc_vvv[p].size();ii++){
          //   YPlaneMasked_h->Fill(masked_adc_vvv[p][ii][1],masked_adc_vvv[p][ii][0],masked_adc_vvv[p][ii][2]);
          // }
          //
          // for (int ii =0;ii<(int)ssnet_sparse_vvv[p].size();ii++){
          //   YPlaneSSNet_h->Fill(ssnet_sparse_vvv[p][ii][1],ssnet_sparse_vvv[p][ii][0],ssnet_sparse_vvv[p][ii][2]);
          // }


          gStyle->SetOptStat(0);

          // vtx_h->Fill(vtx2d_true[3],vtx2d_true[0]);
          std::cout<<"BEFORE FILL"<<std::endl;
          std::vector<double> tmptruestart1 (3,0);
          tmptruestart1[0]= (double)_true_shower_start_vv[0][0];
          tmptruestart1[1]= (double)_true_shower_start_vv[0][1];
          tmptruestart1[2]= (double)_true_shower_start_vv[0][2];
          std::vector<double> tmptruestart2 (3,0);
          tmptruestart2[0]= (double)_true_shower_start_vv[1][0];
          tmptruestart2[1]= (double)_true_shower_start_vv[1][1];
          tmptruestart2[2]= (double)_true_shower_start_vv[1][2];
          std::vector<int> true_shower1_2d = Utils.getProjectedPixel(tmptruestart1, ev_adc->Image2DArray()[2].meta(), 3);
          std::vector<int> true_shower2_2d = Utils.getProjectedPixel(tmptruestart2, ev_adc->Image2DArray()[2].meta(), 3);
          shr_h->Fill(true_shower1_2d[3],true_shower1_2d[0]);
          shr_h->Fill(true_shower2_2d[3],true_shower2_2d[0]);
          std::cout<<"FILLED"<<std::endl;
          triangle1_h->Fill(vtx_pix[0],vtx_pix[1]);
          triangle1_h->Fill(vtx_pix[0] + shlength*cos(shangle+shopen),vtx_pix[1] + shlength*sin(shangle+shopen));
          triangle1_h->Fill(vtx_pix[0] + shlength*cos(shangle-shopen),vtx_pix[1] + shlength*sin(shangle-shopen));
          triangle2_h->Fill(vtx_pix2[0],vtx_pix2[1]);
          triangle2_h->Fill(vtx_pix2[0] + shlength2*cos(shangle2+shopen2),vtx_pix2[1] + shlength2*sin(shangle2+shopen2));
          triangle2_h->Fill(vtx_pix2[0] + shlength2*cos(shangle2-shopen2),vtx_pix2[1] + shlength2*sin(shangle2-shopen2));

          triangle1_h->SetMarkerStyle(kFullCircle);
          triangle1_h->SetMarkerColor(kRed);
          triangle2_h->SetMarkerStyle(kFullCircle);
          triangle2_h->SetMarkerColor(kRed);
          vtx_h->SetMarkerStyle(kStar);
          vtx_h->SetMarkerColor(kGreen);
          shr_h->SetMarkerStyle(kStar);
          shr_h->SetMarkerColor(kMagenta);
          TLine* line1_1 = new TLine(vtx_pix[0],vtx_pix[1],vtx_pix[0] + shlength*cos(shangle+shopen),vtx_pix[1] + shlength*sin(shangle+shopen));
          TLine* line1_2 = new TLine(vtx_pix[0],vtx_pix[1],vtx_pix[0] + shlength*cos(shangle-shopen),vtx_pix[1] + shlength*sin(shangle-shopen));
          TLine* line1_3 = new TLine(vtx_pix[0] + shlength*cos(shangle+shopen),vtx_pix[1] + shlength*sin(shangle+shopen),vtx_pix[0] + shlength*cos(shangle-shopen),vtx_pix[1] + shlength*sin(shangle-shopen));
          TLine* line2_1 = new TLine(vtx_pix2[0],vtx_pix2[1],vtx_pix2[0] + shlength2*cos(shangle2+shopen2),vtx_pix2[1] + shlength2*sin(shangle2+shopen2));
          TLine* line2_2 = new TLine(vtx_pix2[0],vtx_pix2[1],vtx_pix2[0] + shlength2*cos(shangle2-shopen2),vtx_pix2[1] + shlength2*sin(shangle2-shopen2));
          TLine* line2_3 = new TLine(vtx_pix2[0] + shlength2*cos(shangle2+shopen2),vtx_pix2[1] + shlength2*sin(shangle2+shopen2),vtx_pix2[0] + shlength2*cos(shangle2-shopen2),vtx_pix2[1] + shlength2*sin(shangle2-shopen2));
          line1_1->SetLineColor(kRed);
          line1_2->SetLineColor(kRed);
          line1_3->SetLineColor(kRed);
          line2_1->SetLineColor(kMagenta - 3);
          line2_2->SetLineColor(kMagenta - 3);
          line2_3->SetLineColor(kMagenta - 3);

          TCanvas can1("can", "histograms ", 3456, 1008);
          YPlaneADC_h->Draw("colz");
          vtx_h->Draw("SAME");
          shr_h->Draw("SAME");
          triangle1_h->Draw("SAME");
          line1_1->Draw("SAME");
          line1_2->Draw("SAME");
          line1_3->Draw("SAME");
          triangle2_h->Draw("SAME");
          line2_1->Draw("SAME");
          line2_2->Draw("SAME");
          line2_3->Draw("SAME");
          can1.SaveAs(Form("ncpi0_adc_%d_%d.png",(int)entry,(int)ivtx));

          // TCanvas can2("can", "histograms ", 3456, 1008);
          // YPlaneMasked_h->Draw("colz");
          // vtx_h->Draw("SAME");
          // triangle1_h->Draw("SAME");
          // line1_1->Draw("SAME");
          // line1_2->Draw("SAME");
          // line1_3->Draw("SAME");
          // triangle2_h->Draw("SAME");
          // line2_1->Draw("SAME");
          // line2_2->Draw("SAME");
          // line2_3->Draw("SAME");
          // can2.SaveAs(Form("ncpi0_masked_%d.png",(int)entry));

          // TCanvas can3("can", "histograms ", 3456, 1008);
          // YPlaneSSNet_h->Draw("colz");
          // vtx_h->Draw("SAME");
          // triangle1_h->Draw("SAME");
          // line1_1->Draw("SAME");
          // line1_2->Draw("SAME");
          // line1_3->Draw("SAME");
          // triangle2_h->Draw("SAME");
          // line2_1->Draw("SAME");
          // line2_2->Draw("SAME");
          // line2_3->Draw("SAME");
          // can3.SaveAs(Form("ncpi0_ssnet_%d.png",(int)entry));
          //
          // TCanvas can4("can", "histograms ", 3456, 1008);
          // YPlaneSeg_h->Draw("colz");
          // vtx_h->Draw("SAME");
          // triangle1_h->Draw("SAME");
          // line1_1->Draw("SAME");
          // line1_2->Draw("SAME");
          // line1_3->Draw("SAME");
          // triangle2_h->Draw("SAME");
          // line2_1->Draw("SAME");
          // line2_2->Draw("SAME");
          // line2_3->Draw("SAME");
          // can4.SaveAs(Form("ncpi0_seg_%d.png",(int)entry));

          delete vtx_h;
          delete shr_h;
          delete triangle1_h;
          delete triangle2_h;
          delete YPlaneADC_h;
          // delete YPlaneSeg_h;
          // delete YPlaneSSNet_h;
          // delete YPlaneMasked_h;
          delete line1_1;
          delete line1_2;
          delete line1_3;
          delete line2_1;
          delete line2_2;
          delete line2_3;

        }

        //fill triangle points vector

        if (shlength>0){
          std::vector<float> tmp_triangle_v (6,-1);
          tmp_triangle_v[0]=vtx_pix[0];
          tmp_triangle_v[1]=vtx_pix[1];
          tmp_triangle_v[2]=vtx_pix[0] + shlength*cos(shangle+shopen);
          tmp_triangle_v[3]=vtx_pix[1] + shlength*sin(shangle+shopen);
          tmp_triangle_v[4]=vtx_pix[0] + shlength*cos(shangle-shopen);
          tmp_triangle_v[5]=vtx_pix[1] + shlength*sin(shangle-shopen);
          shower_points_vv[p*2] = tmp_triangle_v;
        }
        if (shlength2>0){
          std::vector<float> tmp_triangle_v (6,-1);
          tmp_triangle_v[0]=vtx_pix2[0];
          tmp_triangle_v[1]=vtx_pix2[1];
          tmp_triangle_v[2]=vtx_pix2[0] + shlength2*cos(shangle2+shopen2);
          tmp_triangle_v[3]=vtx_pix2[1] + shlength2*sin(shangle2+shopen2);
          tmp_triangle_v[4]=vtx_pix2[0] + shlength2*cos(shangle2-shopen2);
          tmp_triangle_v[5]=vtx_pix2[1] + shlength2*sin(shangle2-shopen2);
          shower_points_vv[(p*2)+1] = tmp_triangle_v;
        }


      }//end of plane loop

      // store floats
      _shower_start_2d_vvv.push_back( shower_start_2d_vv );
      _shower_gap_vv.push_back( shower_gap_v );
      _shower_energy_vv.push_back( shower_energy_v );
      _shower_sumQ_vv.push_back( shower_sumQ_v );
      _shower_shlength_vv.push_back( shower_shlength_v );
      _shower_shangle_vv.push_back( shower_shangle_v );
      _shower_shopen_vv.push_back( shower_shopen_v );
      _secondshower_start_2d_vvv.push_back( secondshower_start_2d_vv );
      _secondshower_gap_vv.push_back( shower_gap_v );
      _secondshower_energy_vv.push_back( secondshower_energy_v );
      _secondshower_sumQ_vv.push_back( secondshower_sumQ_v );
      _secondshower_shlength_vv.push_back( secondshower_shlength_v );
      _secondshower_shangle_vv.push_back( secondshower_shangle_v );
      _secondshower_shopen_vv.push_back( secondshower_shopen_v );
      _smallQ.push_back( smallq_v );
      _smallQ2.push_back( smallq2_v );

      //do 3d shower match
      if (_second_shower){
        // std::cout<<"HERE"<<std::endl;
        std::vector<std::vector<float>> matchy1;
        std::vector<std::vector<float>> matchy2;
        std::vector<std::vector<int>> bestmatchy1;
        std::vector<std::vector<int>> bestmatchy2;
        matchy1 = SecondShower.Match_3D(shower_points_vv, masked_adc_vvv,1);
        matchy2 = SecondShower.Match_3D(shower_points_vv, masked_adc_vvv,2);
        bestmatchy1 = SecondShower.ChooseBestMatch(matchy1);
        bestmatchy2 = SecondShower.ChooseBestMatch(matchy2);
        _match_y1_vv.push_back(matchy1);
        _match_y2_vv.push_back(matchy2);
        _bestmatch_y1_vv.push_back(bestmatchy1);
        _bestmatch_y2_vv.push_back(bestmatchy2);


        bool RunMassCalc = SecondShower.RunMassCalc(matchy1,matchy2,
              _shower_energy_vv[ivtx],_secondshower_energy_vv[ivtx]);
        if (RunMassCalc) std::cout<<"USE THIS EVENT!! "<<RunMassCalc<<std::endl;
        int useformass;
        if (RunMassCalc) useformass =1;
        else useformass = 0;
        _useformass.push_back(useformass);

        if (RunMassCalc){
          //get 3D overlap points
          //which shower do I use? 0 = u1, 1 = u2, 2 = v1, 3 = v2
          int showertouse1 = -1;
          if (bestmatchy1[0][0] == 2 )showertouse1 = 0;
          if (bestmatchy1[0][1] == 2 )showertouse1 = 1;
          if (bestmatchy1[1][0] == 2 )showertouse1 = 2;
          if (bestmatchy1[1][1] == 2 )showertouse1 = 3;
          int showertouse2 = -1;
          if (bestmatchy2[0][0] == 2 )showertouse2 = 0;
          if (bestmatchy2[0][1] == 2 )showertouse2 = 1;
          if (bestmatchy2[1][0] == 2 )showertouse2 = 2;
          if (bestmatchy2[1][1] == 2 )showertouse2 = 3;
          // std::cout<<"Using showers: y1+"<<showertouse1<<" and y2+"<<showertouse2<<std::endl;
          if (showertouse1 == showertouse2){
            if (bestmatchy2[0][0] == 1 && showertouse1 !=0 )showertouse2 = 0;
            if (bestmatchy2[0][1] == 1 && showertouse1 !=1 )showertouse2 = 1;
            if (bestmatchy2[1][0] == 1 && showertouse1 !=2 )showertouse2 = 2;
            if (bestmatchy2[1][1] == 1 && showertouse1 !=3 )showertouse2 = 3;
            std::cout<<"SWITCHING TO USE y2: "<<showertouse2<<std::endl;
          }
          if (showertouse1 == showertouse2){
            if (bestmatchy1[0][0] == 1 && showertouse2 !=0 )showertouse1 = 0;
            if (bestmatchy1[0][1] == 1 && showertouse2 !=1 )showertouse1 = 1;
            if (bestmatchy1[1][0] == 1 && showertouse2 !=2 )showertouse1 = 2;
            if (bestmatchy1[1][1] == 1 && showertouse2 !=3 )showertouse1 = 3;
            // std::cout<<"SWITCHING TO USE y1: "<<showertouse1<<std::endl;
          }

          //only continue if matched to differentshowers
          if (showertouse1!=showertouse2){
            std::vector<std::vector<float>> FirstShowerPts;
            std::vector<std::vector<float>> SecondShowerPts;
            int planetomatch1 = 0;
            if (showertouse1 > 1)planetomatch1 =1;
            int planetomatch2 = 0;
            if (showertouse2 > 1)planetomatch2 =1;

            FirstShowerPts = SecondShower.Get3DPoints(shower_points_vv[showertouse1],
                              shower_points_vv[4],masked_adc_vvv, planetomatch1,
                              wire_meta[0]);
            SecondShowerPts = SecondShower.Get3DPoints(shower_points_vv[showertouse2],
                              shower_points_vv[5],masked_adc_vvv, planetomatch2,
                              wire_meta[0]);

            std::vector<float> parameters = SecondShower.GetOpeningAngle(FirstShowerPts,SecondShowerPts,_vtx_pos_vv[ivtx]);
            alpha = parameters[0];
            disttoint = (parameters[1]);
            impact1 = (parameters[2]);
            impact2 = (parameters[3]);

            std::vector<std::vector<float>> first = SecondShower.GetPCA(FirstShowerPts);
            std::vector<std::vector<float>> second = SecondShower.GetPCA(SecondShowerPts);
            std::vector<float> firstcenter = first[1];
            std::vector<float> secondcenter = second[1];

            //change direction to diff between pca center and vertex
            firstdirection[0] = firstcenter[0]-_vtx_pos_vv[ivtx][0];
            firstdirection[1] = firstcenter[1]-_vtx_pos_vv[ivtx][1];
            firstdirection[2] = firstcenter[2]-_vtx_pos_vv[ivtx][2];
            seconddirection[0] = secondcenter[0]-_vtx_pos_vv[ivtx][0];
            seconddirection[1] = secondcenter[1]-_vtx_pos_vv[ivtx][1];
            seconddirection[2] = secondcenter[2]-_vtx_pos_vv[ivtx][2];


            if (alpha >.349&&(alpha<3.14 || alpha >3.15)&&alpha<5.93){
               pi0mass = std::sqrt( 4 * _shower_energy_vv[ivtx][2]*_secondshower_energy_vv[ivtx][2] *
                        (sin(alpha/2.0) *sin(alpha/2.0)));
              std::cout<<"Pi0 Mass!!! "<<pi0mass<<std::endl;
            }

          }//end of if showers are different
        }//end of running mass calculation
      }
      _firstdirection.push_back(firstdirection);
      _seconddirection.push_back(seconddirection);
      _pi0mass.push_back(pi0mass);
      _disttoint.push_back(disttoint);
      _impact1.push_back(impact1);
      _impact2.push_back(impact2);
      _alpha.push_back(alpha);

      if (_second_shower && _use_ncpi0){
        //{u1,u2,v1,v2,y1,y2}
        std::vector<float> ShowerTruthMatch_pur_v={-1,-1,-1,-1,-1,-1};
        std::vector<float> ShowerTruthMatch_eff_v={-1,-1,-1,-1,-1,-1};

        std::vector<std::vector<float>> truthmatches(6,std::vector<float> (2,0));
        truthmatches[0] = (SecondShower.TruthMatchNCPi0(shower_points_vv[0], masked_adc_vvv[0],
              ev_segment, ev_instance, 0));
        truthmatches[1] = (SecondShower.TruthMatchNCPi0(shower_points_vv[1], masked_adc_vvv[0],
              ev_segment, ev_instance, 0));
        truthmatches[2] = (SecondShower.TruthMatchNCPi0(shower_points_vv[2], masked_adc_vvv[1],
              ev_segment, ev_instance, 1));
        truthmatches[3] = (SecondShower.TruthMatchNCPi0(shower_points_vv[3], masked_adc_vvv[1],
              ev_segment, ev_instance, 1));
        truthmatches[4] = (SecondShower.TruthMatchNCPi0(shower_points_vv[4], masked_adc_vvv[2],
              ev_segment, ev_instance, 2));
        truthmatches[5] = (SecondShower.TruthMatchNCPi0(shower_points_vv[5], masked_adc_vvv[2],
              ev_segment, ev_instance, 2));
        for(int idx = 0;idx <6;idx++){
          ShowerTruthMatch_pur_v[idx] = truthmatches[idx][0];
          ShowerTruthMatch_eff_v[idx] = truthmatches[idx][1];

        }
        _ShowerTruthMatch_eff_vv.push_back(ShowerTruthMatch_eff_v);
        _ShowerTruthMatch_pur_vv.push_back(ShowerTruthMatch_pur_v);
        //
        // std::cout<<"U PLANE 1 Purity: "<<ShowerTruthMatch_pur_v[0]<<std::endl;
        // std::cout<<"U PLANE 2 Purity: "<<ShowerTruthMatch_pur_v[1]<<std::endl;
        // std::cout<<"U Plane efficiency: "<<ShowerTruthMatch_eff_v[0]+ShowerTruthMatch_eff_v[1]<<std::endl;
        // std::cout<<"V PLANE 1 Purity: "<<ShowerTruthMatch_pur_v[2]<<std::endl;
        // std::cout<<"V PLANE 2 Purity: "<<ShowerTruthMatch_pur_v[3]<<std::endl;
        // std::cout<<"V Plane efficiency: "<<ShowerTruthMatch_eff_v[2]+ShowerTruthMatch_eff_v[3]<<std::endl;
        // std::cout<<"Y PLANE 1 Purity: "<<ShowerTruthMatch_pur_v[4]<<std::endl;
        // std::cout<<"Y PLANE 2 Purity: "<<ShowerTruthMatch_pur_v[5]<<std::endl;
        // std::cout<<"Y Plane efficiency: "<<ShowerTruthMatch_eff_v[4]+ShowerTruthMatch_eff_v[5]<<std::endl;
      }//end of if sec and ncpi0

      if ( _use_nueint){
        //{u1,v1,y1}
        std::vector<float> ShowerTruthMatch_pur_v={-1,-1,-1};
        std::vector<float> ShowerTruthMatch_eff_v={-1,-1,-1};
        std::vector<std::vector<float>> truthmatches(3,std::vector<float> (2,0));
        truthmatches[0] = (SecondShower.TruthMatchNueint(shower_points_vv[0], masked_adc_vvv[0],
              ev_segment, ev_instance,0));
        truthmatches[1] = (SecondShower.TruthMatchNueint(shower_points_vv[2], masked_adc_vvv[1],
              ev_segment, ev_instance,1));
        truthmatches[2] = (SecondShower.TruthMatchNueint(shower_points_vv[4], masked_adc_vvv[2],
              ev_segment, ev_instance,2));

        for(int idx = 0;idx <3;idx++){
          ShowerTruthMatch_pur_v[idx] = truthmatches[idx][0];
          ShowerTruthMatch_eff_v[idx] = truthmatches[idx][1];
        }

        _ShowerTruthMatch_eff_vv.push_back(ShowerTruthMatch_eff_v);
        _ShowerTruthMatch_pur_vv.push_back(ShowerTruthMatch_pur_v);
        // std::cout<<"U PLANE Purity: "<<_ShowerTruthMatch_pur_v[0]<<std::endl;
        // std::cout<<"U Plane efficiency: "<<_ShowerTruthMatch_eff_v[0]<<std::endl;
        // std::cout<<"V PLANE Purity: "<<_ShowerTruthMatch_pur_v[1]<<std::endl;
        // std::cout<<"V Plane efficiency: "<<_ShowerTruthMatch_eff_v[1]<<std::endl;
        // std::cout<<"Y PLANE Purity: "<<_ShowerTruthMatch_pur_v[2]<<std::endl;
        // std::cout<<"Y Plane efficiency: "<<_ShowerTruthMatch_eff_v[2]<<std::endl;
      }//end of if nueint


    }//end of loop through vertices

    //truth matching stuff-----------------------------------------------------
    //only do this if true is in InsideFiducial
    if (_truefid == 1 && ev_mcshower->size() == 2&&_use_mc){
      //loop through vertices again
      std::vector<double> tmptruestart1 (3,0);
      tmptruestart1[0]= (double)_true_shower_start_vv[0][0];
      tmptruestart1[1]= (double)_true_shower_start_vv[0][1];
      tmptruestart1[2]= (double)_true_shower_start_vv[0][2];
      std::vector<double> tmptruestart2 (3,0);
      tmptruestart2[0]= (double)_true_shower_start_vv[1][0];
      tmptruestart2[1]= (double)_true_shower_start_vv[1][1];
      tmptruestart2[2]= (double)_true_shower_start_vv[1][2];
      std::vector<int> true_shower1_tmp_2dstart = Utils.getProjectedPixel(tmptruestart1, ev_adc->Image2DArray()[2].meta(), 3);
      std::vector<int> true_shower2_tmp_2dstart = Utils.getProjectedPixel(tmptruestart2, ev_adc->Image2DArray()[2].meta(), 3);

      for (int vtx =0;vtx<(int)_vtx_pos_vv.size();vtx++){
        std::vector<float> tmp_distances_first = SecondShower.RecoTrueDistances(
              true_shower1_tmp_2dstart, true_shower2_tmp_2dstart, _shower_start_2d_vvv[vtx][2]);
        std::vector<float> tmp_distances_sec = SecondShower.RecoTrueDistances(
              true_shower1_tmp_2dstart, true_shower2_tmp_2dstart, _secondshower_start_2d_vvv[vtx][2]);

        //figure out smallest out of 4 options.
        float tmp_small = 100000;
        int tmp_idx = -1;
        if (tmp_distances_first[0]< tmp_small){
          tmp_small = tmp_distances_first[0];
          tmp_idx = 0;
        }
        if (tmp_distances_first[1]< tmp_small){
          tmp_small = tmp_distances_first[1];
          tmp_idx = 1;
        }
        if (tmp_distances_sec[0]< tmp_small){
          tmp_small = tmp_distances_sec[0];
          tmp_idx = 2;
        }
        if (tmp_distances_sec[1]< tmp_small){
          tmp_small = tmp_distances_sec[1];
          tmp_idx = 3;
        }

        //first do condition of first reco = first true
        if (tmp_idx==0||tmp_idx==1){
          _firstdirection_true.push_back(_true_shower_dir_vv[0]);
          _seconddirection_true.push_back(_true_shower_dir_vv[1]);
          _shower_start_2d_true_vvv.push_back(true_shower1_tmp_2dstart);
          _secondshower_start_2d_true_vvv.push_back(true_shower2_tmp_2dstart);
          _shower_energy_true_vv.push_back(_true_energy_vv[0]);
          _secondshower_energy_true_vv.push_back(_true_energy_vv[1]);
          _shower_recotrue_dist_v.push_back(tmp_distances_first[0]);
          _secondshower_recotrue_dist_v.push_back(tmp_distances_sec[1]);
        }
        //now the case where first reco = second true
        else {
          _firstdirection_true.push_back(_true_shower_dir_vv[1]);
          _seconddirection_true.push_back(_true_shower_dir_vv[0]);
          _shower_start_2d_true_vvv.push_back(true_shower2_tmp_2dstart);
          _secondshower_start_2d_true_vvv.push_back(true_shower1_tmp_2dstart);
          _shower_energy_true_vv.push_back(_true_energy_vv[1]);
          _secondshower_energy_true_vv.push_back(_true_energy_vv[0]);
          _shower_recotrue_dist_v.push_back(tmp_distances_first[1]);
          _secondshower_recotrue_dist_v.push_back(tmp_distances_sec[0]);
        }

      }//end of loop through vtx
    }

    else if (_truefid == 1 && ev_mcshower->size() == 1 && _use_mc){
      //loop through vertices again
      std::vector<double> tmptruestart1 (3,0);
      tmptruestart1[0]= (double)_true_shower_start_vv[0][0];
      tmptruestart1[1]= (double)_true_shower_start_vv[0][1];
      tmptruestart1[2]= (double)_true_shower_start_vv[0][2];
      std::vector<int> true_shower1_tmp_2dstart = Utils.getProjectedPixel(tmptruestart1, ev_adc->Image2DArray()[2].meta(), 3);

      for (int vtx =0;vtx<(int)_vtx_pos_vv.size();vtx++){
        std::vector<float> tmp_distances_first = SecondShower.RecoTrueDistances(
              true_shower1_tmp_2dstart,{-9999,-9999,-9999}, _shower_start_2d_vvv[vtx][2]);
        _firstdirection_true.push_back(_true_shower_dir_vv[0]);
        _shower_start_2d_true_vvv.push_back(true_shower1_tmp_2dstart);
        _shower_energy_true_vv.push_back(_true_energy_vv[0]);
        _shower_recotrue_dist_v.push_back(tmp_distances_first[0]);
      }//end of loop through vtx
    }
    //-------------------------------------------------------------------------
    //OutFile->cd();
    if (_use_mc)
      _numshowers = ev_mcshower->size();
    _ana_tree->Fill();
    // _ana_tree->Write();
    return true;
  }//end of process function

  /*
   * store larlite objects
   *
   * store shower info in _shower_ll_v and _shower_pixcluster_v data members.
   *
   * @param[in]  ioll   larlite for output
   *
   */
  void SSNetShowerReco::store_in_larlite( larlite::storage_manager& ioll ) {

    // entrydata["disttoint"] = []
    //     entrydata["impact1"] = []
    //     entrydata["impact2"] = []
    //     entrydata["alpha"] = []
    //     entrydata["firstdirection"] = []
    //     entrydata["seconddirection"] = []

    larlite::event_shower* evout_shower
      = (larlite::event_shower*)ioll.get_data( larlite::data::kShower, _ssnet_shower_tree_name );
    larlite::event_shower* evout_shower2
      = (larlite::event_shower*)ioll.get_data( larlite::data::kShower, _ssnet_shower_tree_name+"_sec" );
    larlite::event_larflowcluster* evout_lfcluster
      = (larlite::event_larflowcluster*)ioll.get_data( larlite::data::kLArFlowCluster, _ssnet_shower_tree_name );

    if ( _shower_ll_v.size()!=_shower_pixcluster_v.size() ) {
      throw std::runtime_error("[SSNetShwerReco::store_in_larlite] number of larlite shower objects and pixcluster not the same");
    }

    for ( size_t i=0; i<_shower_ll_v.size(); i++ ) {
      evout_shower->push_back(    _shower_ll_v[i] );
      evout_shower2->push_back(    _secondshower_ll_v[i] );
      evout_lfcluster->push_back( _shower_pixcluster_v[i] );
    }


  }

}
}
