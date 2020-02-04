#ifndef __LLCV_SSNET_SHOWER_RECO_H__
#define __LLCV_SSNET_SHOWER_RECO_H__

#include <vector>
#include <string>

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/IOManager.h"

// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/shower.h"
#include "DataFormat/larflowcluster.h"

namespace larlitecv {
namespace ssnetshowerreco {

  class SSNetShowerReco {

  public:

    typedef std::vector< std::vector<float> > triangle_t;
    
    SSNetShowerReco();
    virtual ~SSNetShowerReco() {};

    bool process( larcv::IOManager& iocv, larlite::storage_manager& ioll );//, larcv::IOManager& ioimgs );

    // get functions
    // -------------
    
    ::std::vector< std::vector<float> > getVertexShowerEnergies() const { return _shower_energy_vv; };
    int   numVertices() const { return _shower_energy_vv.size(); };
    float getVertexShowerEnergy( int vtxid, int plane ) const { return _shower_energy_vv[vtxid][plane]; };
    float getVertexShowerSumQ( int vtxid, int plane ) const { return _shower_sumQ_vv[vtxid][plane]; };
    float getVertexShowerShlength( int vtxid, int plane ) const { return _shower_shlength_vv[vtxid][plane]; };
    ::std::vector<double> getVertexPos( int vtxid ) const { return _vtx_pos_vv[vtxid]; };

    larlite::shower& getShowerObject( int vtxid, int plane ) { return _shower_ll_v.at( 3*vtxid+plane ); };
    larlite::larflowcluster& getShowerPixelList( int vtxid, int plane ) { return _shower_pixcluster_v.at( 3*vtxid+plane ); };

    void store_in_larlite( larlite::storage_manager& ioll );
    void use_calibrated_pixsum2mev( bool use=true ) { _use_calibrated_pixelsum2mev = use; };

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
    
    std::vector< std::vector<int> > _enclosedCharge( const larcv::Image2D& chargeMap,
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

    // parameters
    // -----------
    std::string _adc_tree_name;
    std::string _calib_adc_tree_name;
    std::string _ssnet_shower_image_stem; // sometimes ubspurn (when files made at FNAL)
    std::string _vertex_tree_name;
    std::string _track_tree_name;
    float _Qcut;
    float _SSNET_SHOWER_THRESHOLD;
    bool  _use_calibrated_pixelsum2mev;
    
  public:
    // set methods
    // ------------
    void set_adc_treename( std::string name )          { _adc_tree_name = name; };
    void set_ssnet_shower_stemname( std::string name ) { _ssnet_shower_image_stem = name; };
    void set_vertex_treename( std::string name )       { _vertex_tree_name = name; };
    void set_track_treename( std::string name )        { _track_tree_name = name; };
    void set_Qcut( float cut )                         { _Qcut = cut; };
    void set_SSNet_threshold( float threshold )        { _SSNET_SHOWER_THRESHOLD = threshold; };


  protected:
    // variables filled
    // -----------------
    std::vector< std::vector<float> >  _shower_energy_vv;
    std::vector< std::vector<float> >  _shower_sumQ_vv;
    std::vector< std::vector<float> >  _shower_shlength_vv;
    std::vector< std::vector<double> > _vtx_pos_vv;
    std::vector< larlite::shower >         _shower_ll_v;
    std::vector< larlite::larflowcluster > _shower_pixcluster_v;
  public:
    void clear();
    

  };
  
  
}
}

#endif
