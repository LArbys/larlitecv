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
  SSNetShowerReco::SSNetShowerReco() {

    _adc_tree_name = "wire";
    _calib_adc_tree_name = "calibrated";
    _ssnet_shower_image_stem = "uburn"; // sometimes ubspurn (when files made at FNAL)
    _vertex_tree_name = "test";
    _track_tree_name  = "trackReco";
    _Qcut = 10;
    _SSNET_SHOWER_THRESHOLD = 0.05;
    _use_calibrated_pixelsum2mev = false;
    clear();
  }

  /**
   * clear result containers
   *
   */
  void SSNetShowerReco::clear() {
    _shower_energy_vv.clear();
    _shower_sumQ_vv.clear();
    _shower_shlength_vv.clear();
    _vtx_pos_vv.clear();
    _shower_ll_v.clear();
    _shower_pixcluster_v.clear();
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
                                                                     int vtx_col,
                                                                     int vtx_row,
                                                                     float shLen,
                                                                     float shOpen ) {

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
      for (int ii =0;ii<chargeMap.size();ii++){
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
    for (int ii =0;ii<chargeMap.size();ii++){
      if ( _isInside2(t1X,t1Y,t2X,t2Y,t3X,t3Y,chargeMap[ii][1],chargeMap[ii][0])){
        sumIn+=chargeMap[ii][2];
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
      _enclosedCharge( chargeMap, ang, sumQ, triangle,
                       vtx_col,vtx_row,scanLen,scanOpen);
      //std::cout << " find-dir ang=" << ang*180.0/TMath::Pi() << " deg; =" << ang << " rad;  sumQ=" << sumQ << std::endl;
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
      _enclosedCharge( chargeMap, ang, sumQ, triangle,
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
    std::vector<float> lengths( (int)350/step, 0.0);
    for (int i=0; i<(int)350/step; i++ ) {
      lengths[i]  = step*i;
    }
    std::vector<float> charge;
    charge.reserve( lengths.size() );

    for ( auto const& length : lengths ) {
      float sumQ;
      triangle_t triangle;
      _enclosedCharge(chargeMap,theta,
                      sumQ, triangle,
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
      _enclosedCharge(chargeMap,theta,sumQ,triangle,
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
            float threshold ){
     //function to make sparse objects from the larcv Image2d
     //output: vector of (row,col, value)
     std::vector<std::vector<float>> output_vv;
     for (int r = 0; r<input_img.meta().rows();r++){
       for (int c= 0; c<input_img.meta().cols();c++){
         float pix_value = input_img.pixel(r,c);
         if (pix_value > threshold){
           std::vector<float> tmp_v = {(float) r, (float) c, pix_value};
           output_vv.push_back(tmp_v);
         }
       }//end of col loop
     }//end of row loop
     return output_vv;
   }// end of function

  bool SSNetShowerReco::process( larcv::IOManager& iolcv, larlite::storage_manager& ioll, int entry ) {//, larcv::IOManager& ioimgs ) {

    // get adc image (larcv)
    // get ssnet image (larcv)
    // get vertex (larcv)
    // get tracks
    // store in root ana tree, store in json file

    // parameters for later
    // std::string adc_tree_name = "wire";
    // std::string ssnet_shower_image_stem = "uburn"; // sometimes ubspurn (when files made at FNAL)
    // std::string vertex_tree_name = "test";
    // std::string track_tree_name  = "trackReco";
    // float Qcut = 10;
    // float Scut = 0.05;
    // float SSNET_SHOWER_THRESHOLD = 0.5;

    // clear result container
    clear();

    larcv::EventImage2D* ev_adc
      = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _adc_tree_name );
    const std::vector<larcv::Image2D>& adc_v = ev_adc->Image2DArray();

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
    std::vector<std::vector<float>> ssnet_u_shower_sparse = MakeImage2dSparse(ev_shower_score[0]->Image2DArray()[0],_SSNET_SHOWER_THRESHOLD);
    std::vector<std::vector<float>> ssnet_v_shower_sparse = MakeImage2dSparse(ev_shower_score[1]->Image2DArray()[0],_SSNET_SHOWER_THRESHOLD);
    std::vector<std::vector<float>> ssnet_y_shower_sparse = MakeImage2dSparse(ev_shower_score[2]->Image2DArray()[0],_SSNET_SHOWER_THRESHOLD);
    std::vector<std::vector<std::vector<float>>> ssnet_sparse_vvv = {ssnet_u_shower_sparse,ssnet_v_shower_sparse,ssnet_y_shower_sparse};
    larcv::EventPGraph* ev_vtx
      = (larcv::EventPGraph*)iolcv.get_data( larcv::kProductPGraph, _vertex_tree_name );

    std::vector<std::vector<std::vector<float>>> masked_adc_vvv;

    std::cout<<"Size of input sparse vectors: \n";
    std::cout<<"ADC: "<<adc_sparse_vvv[0].size()<<","<<adc_sparse_vvv[1].size()<<","<<adc_sparse_vvv[2].size()<<"\n";
    std::cout<<"SSNET: "<<ssnet_sparse_vvv[0].size()<<","<<ssnet_sparse_vvv[1].size()<<","<<ssnet_sparse_vvv[2].size()<<"\n";
    for ( size_t p=0; p<3; p++ ) {
      // mask the ADC image using ssnet
      int npixels_adc = adc_sparse_vvv[p].size();
      int npixels_ssnet = ssnet_sparse_vvv[p].size();
      std::vector<std::vector<float>> masked_plane_v;
      for ( int adc_i =0;adc_i<npixels_adc;adc_i++ ) {
        for ( int ssnet_i =0;ssnet_i<npixels_ssnet;ssnet_i++){
          if (ssnet_sparse_vvv[p][ssnet_i][0] == adc_sparse_vvv[p][adc_i][0] && ssnet_sparse_vvv[p][ssnet_i][1] == adc_sparse_vvv[p][adc_i][1]){
              masked_plane_v.push_back({(float) adc_sparse_vvv[p][adc_i][0],adc_sparse_vvv[p][adc_i][1],adc_sparse_vvv[p][adc_i][2]});
          }//end of if statements
        }//end of loop through ssnet
      }//end of loop throug adc
      masked_adc_vvv.push_back(masked_plane_v);
    }//end of plane loop
    std::cout<<"Size of masked adc sparse vectors: \n";
    std::cout<<"Masked ADC: "<<masked_adc_vvv[0].size()<<","<<masked_adc_vvv[1].size()<<","<<masked_adc_vvv[2].size()<<"\n";
    // for(int ii =0;ii<masked_adc_vvv[0].size();ii++){
    //   std::cout<<masked_adc_vvv[0][ii][0]<<","<<masked_adc_vvv[0][ii][1]<<std::endl;
    // }
    // get candidate vertices, if in fiducial volume, keep
    _vtx_pos_vv.clear();
    for ( auto const& pgraph : ev_vtx->PGraphArray() ) {
      if ( pgraph.ParticleArray().size()==0 ) continue; // dont expect this
      auto const& roi = pgraph.ParticleArray().front();
      std::vector<double> vtx3d = { roi.X(), roi.Y(), roi.Z() };
      std::cout << "Vertex Pos (" << vtx3d[0] << "," << vtx3d[1] << "," << vtx3d[2] << ")" << std::endl;
      if ((vtx3d[0]>    0.001) && (vtx3d[0] <  255.999) && (vtx3d[1]> -116.499) && (vtx3d[1] < 116.499)
          && (vtx3d[2]>    0.001) && (vtx3d[2] < 1036.999)){

        _vtx_pos_vv.push_back( vtx3d );


      }//end of if statement
      else std::cout<<"Outside fiducial volume!"<<std::endl;
    }//end of loop over vertices/pgraph

    // now get shower energy per vertex, per plane

    for ( size_t ivtx=0; ivtx<_vtx_pos_vv.size(); ivtx++ ) {


      std::cout << "[SSNetShowerReco] Reconstruct vertex[" << ivtx << "] "
                << "pos=(" << _vtx_pos_vv[ivtx][0] << "," << _vtx_pos_vv[ivtx][1] << "," << _vtx_pos_vv[ivtx][2] << ")"
                << std::endl;

      //get 2d vertex position
      std::vector<int> imgcoord_v(4);  // (U,V,Y,tick)
      imgcoord_v[3] = 3200 + _vtx_pos_vv[ivtx][0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
      for ( size_t p=0; p<3; p++ ) {
        imgcoord_v[p] = larutil::Geometry::GetME()->NearestWire( _vtx_pos_vv[ivtx], (int)p );
      }

      std::vector<float> shower_energy_v(3,0);
      std::vector<float> shower_sumQ_v(3,0);
      std::vector<float> shower_shlength_v(3,0);

      for ( size_t p=0; p<3; p++ ) {
        //std::cout << "[SSNetShowerReco]   Plane [" << p << "]" << std::endl;
        int vtx_pix[2] = { (int)wire_meta[p].col(imgcoord_v[p]) , (int)wire_meta[p].row(imgcoord_v[3])};
        //std::cout << "[SSNetShowerReco]     vertex pixel: (" << vtx_pix[0] << "," << vtx_pix[1] << ")" << std::endl;
        float shangle  = _findDir( masked_adc_vvv[p], vtx_pix[0], vtx_pix[1] );
        //std::cout << "[SSNetShowerReco]     shower angle: " << shangle << std::endl;
        float shlength = _findLen( masked_adc_vvv[p],  shangle, vtx_pix[0], vtx_pix[1] );
        //std::cout << "[SSNetShowerReco]     shower length: " << shlength << std::endl;
        float shopen   = _findOpen( masked_adc_vvv[p], shangle, shlength, vtx_pix[0], vtx_pix[1] );
        //std::cout << "[SSNetShowerReco]     open angle: " << shopen << std::endl;
        float sumQ;
        triangle_t tri;
        std::vector< std::vector<int> > pixlist_v =
          _enclosedCharge( masked_adc_vvv[p], shangle, sumQ, tri, vtx_pix[0], vtx_pix[1], shlength, shopen );

        float reco_energy = 0;
        if ( _use_calibrated_pixelsum2mev )
          reco_energy = sumQ*0.01324 + 37.83337; // calibrated images
        else
          reco_energy = sumQ*0.013456 + 2.06955; // uncalibrated images

        std::cout << "[SSNetShowerReco] plane[" << p << "] final sumQ=" << sumQ << " reco=" << reco_energy << std::endl;

        shower_energy_v[p] = reco_energy;
        shower_sumQ_v[p] = sumQ;
        shower_shlength_v[p] = shlength;

        // make larlite
        // -------------
        larlite::shower shr;  // stores shower parameters
        larlite::larflowcluster pixcluster; // stores pixels

        // store vertex index
        shr.set_id( ivtx );

        // store vertex
        TVector3 vtx3 = { _vtx_pos_vv[ivtx][0], _vtx_pos_vv[ivtx][1], _vtx_pos_vv[ivtx][2] };
        shr.set_start_point( vtx3 );

        // store plane
        shr.set_total_best_plane( p );

        // set energy
        shr.set_total_energy( reco_energy );

        // set sumq
        shr.set_total_energy_err( sumQ );

        // set length
        shr.set_length( shlength );

        // opening angle
        shr.set_opening_angle( shopen );

        // for direction store 3 vector
        // (x,y) are the wire and tick angle. (z) is the angle itself
        TVector3 pseudo_dir = { cos( shangle ), sin( shangle ), shangle };
        shr.set_direction( pseudo_dir );

        _shower_ll_v.emplace_back( std::move(shr) );

        pixcluster.reserve( pixlist_v.size() );
        for ( auto& pix : pixlist_v ) {
          larlite::larflow3dhit lfpix;
          lfpix.resize(2,0);
          lfpix[0] = pix[0];
          lfpix[1] = pix[1];
          pixcluster.emplace_back( std::move(lfpix) );
        }
        _shower_pixcluster_v.emplace_back( std::move(pixcluster) );
      }//end of plane loop

      // store floats
      _shower_energy_vv.push_back( shower_energy_v );
      _shower_sumQ_vv.push_back( shower_sumQ_v );
      _shower_shlength_vv.push_back( shower_shlength_v );

    }//end of loop through vertices

    TH2F* YPlaneADC_h = new TH2F("yplaneadc","yplaneadc",3456,0,3456.,1008,0,1008.);
    TH2F* YPlaneSSNet_h = new TH2F("yplanessnet","yplanessnet",3456,0,3456.,1008,0,1008.);
    TH2F* YPlaneMasked_h = new TH2F("yplanemasked","yplanemasked",3456,0,3456.,1008,0,1008.);


    // for (int ii =0;ii<adc_sparse_vvv[2].size();ii++){
    //   YPlaneADC_h->Fill(adc_sparse_vvv[2][ii][1],adc_sparse_vvv[2][ii][0],adc_sparse_vvv[2][ii][2]);
    // }
    // for (int ii =0;ii<ssnet_sparse_vvv[2].size();ii++){
    //   YPlaneSSNet_h->Fill(ssnet_sparse_vvv[2][ii][1],ssnet_sparse_vvv[2][ii][0],ssnet_sparse_vvv[2][ii][2]);
    // }
    // for (int ii =0;ii<masked_adc_vvv[2].size();ii++){
    //   YPlaneMasked_h->Fill(masked_adc_vvv[2][ii][1],masked_adc_vvv[2][ii][0],masked_adc_vvv[2][ii][2]);
    // }

    for (int r =0;r<adc_v[2].meta().rows();r++){
      for (int c =0;c<adc_v[2].meta().cols();c++){
        YPlaneADC_h->Fill((float)c,(float)r, adc_v[2].pixel(r,c));
      }
    }

    gStyle->SetOptStat(0);

    TCanvas can1("can", "histograms ", 3456, 1008);
    YPlaneADC_h->Draw("colz");
    can1.SaveAs(Form("sparsetest_comb_adc_%d.png",(int)entry));
    //
    // TCanvas can2("can", "histograms ", 3456, 1008);
    // YPlaneSSNet_h->Draw("colz");
    // can2.SaveAs(Form("sparsetest_comb_ssnet_%d.png",(int)entry));
    //
    // TCanvas can3("can", "histograms ", 3456, 1008);
    // YPlaneMasked_h->Draw("colz");
    // can3.SaveAs(Form("sparsetest_comb_masked_%d.png",(int)entry));

    return true;
  }

  /*
   * store larlite objects
   *
   * store shower info in _shower_ll_v and _shower_pixcluster_v data members.
   *
   * @param[in]  ioll   larlite for output
   *
   */
  void SSNetShowerReco::store_in_larlite( larlite::storage_manager& ioll ) {

    larlite::event_shower* evout_shower
      = (larlite::event_shower*)ioll.get_data( larlite::data::kShower, "ssnetshowerreco" );
    larlite::event_larflowcluster* evout_lfcluster
      = (larlite::event_larflowcluster*)ioll.get_data( larlite::data::kLArFlowCluster, "ssnetshowerreco" );

    if ( _shower_ll_v.size()!=_shower_pixcluster_v.size() ) {
      throw std::runtime_error("[SSNetShwerReco::store_in_larlite] number of larlite shower objects and pixcluster not the same");
    }

    for ( size_t i=0; i<_shower_ll_v.size(); i++ ) {
      evout_shower->push_back(    _shower_ll_v[i] );
      evout_lfcluster->push_back( _shower_pixcluster_v[i] );
    }

  }

}
}
