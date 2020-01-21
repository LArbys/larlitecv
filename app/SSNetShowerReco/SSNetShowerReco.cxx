#include "SSNetShowerReco.h"

#include "TMath.h"

#include "DataFormat/EventImage2D.h"
#include "DataFormat/track.h"

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


  void SSNetShowerReco::process( larcv::IOManager& iolcv, larlite::storage_manager& ioll ) {

    // get adc image
    // get ssnet image
    // get vertex
    // get tracks
  
  /*
Qcut = 10
Scut = 0.05

# IN WHATEVER FASHION IS EASIEST, PRODUCE THE VARIABLE "showerRawVals" WHICH MUST BE A LIST OF LISTS, WITH FORMAT
# [ [X_PIXEL_POSITION, Y_PIXEL_POSITION , ADC_CHARGE] ... ] WHERE ONLY PIXELS WITH SSNET SHOWER SCORE > "Scut" AND ACD > "Qcut"
# ARE NON_ZERO. X AND Y ARE IN 2D IMAGE COORDINATES. I'VE BEEN USING 512X512 IMAGES WITH THE VERTEX CENTERED.
# SIMILARLY DEFINE THE VARIABLES "vtx_x" and "vtx_y" WHICH IS THE VERTEX POSITION IN 2D IMAGE COORDINATES

shangle       = FindDir(showerRawVals,vtx=[vtx_x,vtx_y])  
shlength      = FindLen(showerRawVals,shangle,vtx=[vtx_x,vtx_y])
shopen        = FindOpen(showerRawVals,shangle,length=shlength,vtx=[vtx_x,vtx_y])
sumQ,triangle = EnclosedCharge(showerRawVals,shangle,vtx=[vtx_x,vtx_y],shOpen=shopen,shLen=length)

reco_shower_energy = sumQ*0.0115 + 50
*/

    return;
  }

}
}
