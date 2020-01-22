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
   * use triple product to get area of triangle
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
   * get the enclosed charge inside the image.
   */
  void SSNetShowerReco::_enclosedCharge( const larcv::Image2D& chargeMap,
                                         float theta,
                                         float& sumIn,
                                         std::vector< std::vector<float> >& triangle,
                                         int vtx_col,
                                         int vtx_row,
                                         float shLen,
                                         float shOpen ) {
                                         
    std::vector<int> vtx = { vtx_col, vtx_row };
    float t1X = vtx[0];
    float t1Y = vtx[1];
    
    float t2X = vtx[0] + shLen*cos(theta+shOpen);
    float t2Y = vtx[1] + shLen*sin(theta+shOpen);
    
    float t3X = vtx[0] + shLen*cos(theta-shOpen);
    float t3Y = vtx[1] + shLen*sin(theta-shOpen);
    
    sumIn = 0;
    for ( size_t r=0; r<chargeMap.meta().rows(); r++ ) {
      for ( size_t c=0; c<chargeMap.meta().cols(); c++ ) {
        if ( _isInside(t1X,t1Y,t2X,t2Y,t3X,t3Y,(float)c,(float)r) )
          sumIn += chargeMap.pixel(r,c);
      }
    }

    triangle.resize(3);
    for (int i=0; i<3; i++ ) {
      triangle[i].resize(2);
    }
    triangle[0] = { t1X, t1Y };
    triangle[1] = { t2X, t2Y };
    triangle[2] = { t3X, t3Y };
  
    return;
  }
    

  float SSNetShowerReco::_findDir( const larcv::Image2D& chargeMap,
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

      coarseAngs[i]  = 2*TMath::Pi() / (float)coarseSteps * i;
      float ang = coarseAngs[i];
      float sumQ;
      std::vector< std::vector<float> > triangle;
      _enclosedCharge( chargeMap, ang, sumQ, triangle,
                       vtx_col,vtx_row,scanLen,scanOpen);
      if ( sumQ > maxCharge ) {
        maxCharge = sumQ;
        bestDir   = ang;
        bestTri   = triangle;
      }
    }
            
    // Fine Tune
    float fineBestDir = -9999;

    for (int i=-5; i<=5; i++ ) {
      //angs = [bestDir + pi/180*i for i in range(-5,5)]
      float ang = bestDir + TMath::Pi()*i/180.0;
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


  float SSNetShowerReco::_findLen( const larcv::Image2D& chargeMap,
                                   float theta,
                                   int vtx_col, int vtx_row,
                                   float scanOpen ) {
    
    const int step     = 10;
    std::vector<float> lengths( (int)350/step);
    for (int i=0; i<(int)350/step; i++ ) {
      lengths[i]  = step*i;
    }
    std::vector<float> charge;
    charge.reserve( lengths.size() );
    
    for ( auto const& length : lengths ) {
      float sumQ;
      std::vector< std::vector<float> > triangle;
      _enclosedCharge(chargeMap,theta,
                      sumQ, triangle,
                      vtx_col, vtx_row,
                      length,scanOpen);
      charge.push_back( sumQ );
    }

    float maxCharge = *std::max_element( charge.begin(), charge.end() );
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

  float SSNetShowerReco::_findOpen( const larcv::Image2D& chargeMap,
                                    float theta,
                                    float length,
                                    int vtx_col, int vtx_row ) {
    
    float step    = 0.02;
    std::vector<float> opens( 14, 0 );
    for (int i=1; i<15; i++ ) {
      opens[i] = step*i;
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


  bool SSNetShowerReco::process( larcv::IOManager& iolcv, larlite::storage_manager& ioll ) {

    // get adc image (larcv)
    // get ssnet image (larcv)
    // get vertex (larcv)
    // get tracks
    // store in root ana tree, store in json file

    // parameters for later
    std::string adc_tree_name = "wire";
    std::string ssnet_shower_image_stem = "ubspurn"; // sometimes uburn (when files made at FNAL)
    std::string vertex_tree_name = "test";
    std::string track_tree_name  = "trackReco";
    float Qcut = 10;
    float Scut = 0.05;
    float SSNET_SHOWER_THRESHOLD = 0.5;

    larcv::EventImage2D* ev_adc
      = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, adc_tree_name );
    const std::vector<larcv::Image2D>& adc_v = ev_adc->Image2DArray();

    larcv::EventImage2D* ev_shower_score[3] = { nullptr };
    for ( size_t p=0; p<3; p++ ) {
      char treename[50];
      sprintf( treename, "%s_plane%d", ssnet_shower_image_stem.c_str(), (int)p );
      ev_shower_score[p] = 
        (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, treename );
    }

    larcv::EventPGraph* ev_vtx
      = (larcv::EventPGraph*)iolcv.get_data( larcv::kProductPGraph, vertex_tree_name );

    
    // get candidate vertices, make crops around said vertex
    std::vector< std::vector<larcv::Image2D> > adc_vtxcrop_vv;
    std::vector< std::vector<larcv::Image2D> > shower_vtxcrop_vv;
    std::vector< std::vector<int> >            vtx_incrop_vv;
    for ( auto const& pgraph : ev_vtx->PGraphArray() ) {
      if ( pgraph.ParticleArray().size()==0 ) continue; // dont expect this
      auto const& roi = pgraph.ParticleArray().front();
      std::vector<double> vtx3d = { roi.X(), roi.Y(), roi.Z() };
      
      std::vector<int> imgcoord_v(4);  // (U,V,Y,tick)
      imgcoord_v[3] = 3200 + vtx3d[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;      
      for ( size_t p=0; p<3; p++ ) {
        imgcoord_v[p] = larutil::Geometry::GetME()->NearestWire( vtx3d, (int)p );
      }

      // define crop
      std::vector<larcv::Image2D> crop_v;
      std::vector<larcv::Image2D> sscrop_v;
      for ( size_t p=0; p<3; p++ ) {
        auto const& meta = adc_v[p].meta();
        int tbounds[2] = { (int)(imgcoord_v[3]-256*meta.pixel_height()),
                           (int)(imgcoord_v[3]+256*meta.pixel_height()) };
        int wbounds[2] = { (int)(imgcoord_v[p]-256*meta.pixel_width()),
                           (int)(imgcoord_v[p]+256*meta.pixel_width()) };
        
        // correct bounds
        if ( tbounds[0]<=meta.min_y() ) {
          tbounds[0] = meta.min_y();
          tbounds[1] = tbounds[0] + 512*meta.pixel_height();
        }
        if ( tbounds[1]>=meta.max_y() ) {
          tbounds[1] = meta.max_y();
          tbounds[0] = tbounds[1] - 512*meta.pixel_height();
        }
        if ( wbounds[0]<=meta.min_x() ) {
          wbounds[0] = meta.min_x();
          wbounds[1] = wbounds[0] + 512*meta.pixel_width();
        }
        if ( wbounds[1]>=meta.max_x() ) {
          wbounds[1] = meta.max_x();
          wbounds[0] = wbounds[1] - 512*meta.pixel_width();
        }


        // define meta for crop
        larcv::ImageMeta cropmeta( 512*meta.pixel_width(), 512*meta.pixel_height(),
                                   512, 512, wbounds[0], tbounds[1],
                                   meta.plane() );
        larcv::Image2D crop = adc_v[p].crop(cropmeta);
        larcv::Image2D sscrop = ev_shower_score[p]->Image2DArray()[p].crop(cropmeta);

        // mask the ADC image using ssnet
        for ( size_t r=0; r<crop.meta().rows(); r++ ) {
          for ( size_t c=0; c<crop.meta().cols(); c++ ) {
            if ( sscrop.pixel(r,c)<SSNET_SHOWER_THRESHOLD ) {
              crop.set_pixel(r,c,0.0);
            }
          }
        }
        
        crop_v.emplace_back( std::move(crop) );
        sscrop_v.emplace_back( std::move(sscrop) );
        
      }//end of plane loop

      
      adc_vtxcrop_vv.emplace_back( std::move(crop_v) );
      shower_vtxcrop_vv.emplace_back( std::move(sscrop_v) );
      vtx_incrop_vv.push_back( imgcoord_v );
    }//end of loop over vertices/pgraph

    //
    
    // now get shower energy per vertex, per plane

    std::vector< std::vector<float> > shower_energy_vv;
    for ( size_t ivtx=0; ivtx<adc_vtxcrop_vv.size(); ivtx++ ) {
      auto& crop_v     = adc_vtxcrop_vv[ivtx];
      auto& imgcoord_v = vtx_incrop_vv[ivtx];

      std::vector<float> shower_energy_v(3,0);
      for ( size_t p=0; p<3; p++ ) {
        int vtx_pix[2] = { crop_v[p].meta().col( imgcoord_v[p] ),
                           crop_v[p].meta().row( imgcoord_v[3] ) };
        float shangle  = _findDir( crop_v[p], vtx_pix[0], vtx_pix[1] );
        float shlength = _findLen( crop_v[p],  shangle, vtx_pix[0], vtx_pix[1] );
        float shopen   = _findOpen( crop_v[p], shangle, shlength, vtx_pix[0], vtx_pix[1] );
        float sumQ;
        triangle_t tri;
        _enclosedCharge( crop_v[p], shangle, sumQ, tri, vtx_pix[0], vtx_pix[1], shopen, shlength );

        float reco_energy = sumQ*0.0115 + 50.0;
        shower_energy_v[p] = reco_energy;
      }
      
      shower_energy_vv.push_back( shower_energy_v );
      
    }


    // ok now store

    return true;
  }

}
}
