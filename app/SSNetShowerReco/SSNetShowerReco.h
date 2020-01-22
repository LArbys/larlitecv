#ifndef __LLCV_SSNET_SHOWER_RECO_H__
#define __LLCV_SSNET_SHOWER_RECO_H__

#include <vector>
#include "DataFormat/Image2D.h"
#include "DataFormat/IOManager.h"
#include "DataFormat/storage_manager.h"

namespace larlitecv {
namespace ssnetshowerreco {

  class SSNetShowerReco {

  public:

    typedef std::vector< std::vector<float> > triangle_t;
    
    SSNetShowerReco() {};
    virtual ~SSNetShowerReco() {};

    bool process( larcv::IOManager& iocv, larlite::storage_manager& ioll );

    ::std::vector< std::vector<float> > getVertexShowerEnergies() const { return _shower_energy_vv; };
    int   numVertices() const { return _shower_energy_vv.size(); };
    float getVertexShowerEnergy( int vtxid, int plane ) const { return _shower_energy_vv[vtxid][plane]; };

  protected:
    
    float _area( float x1, float y1,
                 float x2, float y2,
                 float x3, float y3 );

    float _sign( float x1, float y1,
                 float x2, float y2,
                 float x3, float y3 );
    
    bool _isInside(float x1, float y1,
                   float x2, float y2,
                   float x3, float y3,
                   float x, float y );

    bool _isInside2(float x1, float y1,
                    float x2, float y2,
                    float x3, float y3,
                    float x, float y );
    
    void _enclosedCharge( const larcv::Image2D& chargeMap,
                          float theta,
                          float& sumIn,
                          std::vector< std::vector<float> >& triangle,
                          int vtx_col=255,
                          int vtx_row=255,
                          float shLen = 100.0,
                          float shOpen = 0.2);

    float _findDir( const larcv::Image2D& chargeMap,
                    int vtx_col=255,
                    int vtx_row=255,
                    float scanLen = 50,
                    float scanOpen=0.05 );

    float _findLen( const larcv::Image2D& chargeMap,
                    float theta,
                    int vtx_col=255, int vtx_row=255,
                    float scanOpen=0.2);

    float _findOpen( const larcv::Image2D& chargeMap,
                     float theta,
                     float length,
                     int vtx_col=255, int vtx_row=255 );

    
    
    
    std::vector< std::vector<float> > _shower_energy_vv;
    
    

  };
  
  
}
}

#endif
