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

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TLine.h"

#include "SecondShower.h"
#include "Utils.h"

namespace larlitecv {
namespace ssnetshowerreco {

  class SSNetShowerReco {

  public:

    typedef std::vector< std::vector<float> > triangle_t;

    SSNetShowerReco();
    virtual ~SSNetShowerReco() {};

    bool process( larcv::IOManager& iocv, larlite::storage_manager& ioll, int entry);//, larcv::IOManager& ioimgs );

    // get functions
    // -------------

    ::std::vector< std::vector<float> > getVertexShowerEnergies() const { return _shower_energy_vv; };
    int   numVertices() const { return _shower_energy_vv.size(); };
    int   numShowers() const { return _true_energy_vv.size(); };
    float getVertexShowerEnergy( int vtxid, int plane ) const { return _shower_energy_vv[vtxid][plane]; };
    float getVertexShowerSumQ( int vtxid, int plane ) const { return _shower_sumQ_vv[vtxid][plane]; };
    float getVertexShowerShlength( int vtxid, int plane ) const { return _shower_shlength_vv[vtxid][plane]; };
    float getVertexSecondShowerEnergy( int vtxid, int plane ) const { return _secondshower_energy_vv[vtxid][plane]; };
    float getVertexSecondShowerSumQ( int vtxid, int plane ) const { return _secondshower_sumQ_vv[vtxid][plane]; };
    float getVertexSecondShowerShlength( int vtxid, int plane ) const { return _secondshower_shlength_vv[vtxid][plane]; };
    ::std::vector<double> getVertexPos( int vtxid ) const { return _vtx_pos_vv[vtxid]; };
    double getTrueShowerEnergy(int shower) const { return _true_energy_vv[shower]; };
    ::std::vector<double> getTrueShowerStarts(int shower) const { return _true_shower_start_vv[shower]; };
    double getRemainingADC() const { return _remaining_adc; };
    float getOverlapFraction1(int plane, int shower ) const { return _match_y1_vv[plane][shower]; };
    float getOverlapFraction2(int plane, int shower ) const { return _match_y2_vv[plane][shower]; };
    int numpointsU() const {return _uplane_profile_vv.size(); };
    int numpointsV() const {return _vplane_profile_vv.size(); };
    int numpointsY() const {return _yplane_profile_vv.size(); };
    float getUPlaneShowerProfile( int pt, int ii ) const { return _uplane_profile_vv[pt][ii]; };
    float getVPlaneShowerProfile( int pt, int ii ) const { return _vplane_profile_vv[pt][ii]; };
    float getYPlaneShowerProfile( int pt, int ii ) const { return _yplane_profile_vv[pt][ii]; };
    float getShowerTruthMatch (int shower) const {return _ShowerTruthMatch_v[shower]; };

    larlite::shower& getShowerObject( int vtxid, int plane ) { return _shower_ll_v.at( 3*vtxid+plane ); };
    larlite::larflowcluster& getShowerPixelList( int vtxid, int plane ) { return _shower_pixcluster_v.at( 3*vtxid+plane ); };

    void store_in_larlite( larlite::storage_manager& ioll );
    void use_calibrated_pixsum2mev( bool use=true ) { _use_calibrated_pixelsum2mev = use; };
    void use_second_shower( bool use=true ) { _second_shower = use; };
    void use_ncpi0( bool use=true ) { _use_ncpi0 = use; };
    void use_nueint( bool use=true ) { _use_nueint = use; };

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

    std::vector< std::vector<int> > _enclosedCharge( std::vector<std::vector<float>> chargeMap,
                                                     float theta,
                                                     float& sumIn,
                                                     std::vector< std::vector<float> >& triangle,
                                                     int vtx_col=255,
                                                     int vtx_row=255,
                                                     float shLen = 100.0,
                                                     float shOpen = 0.3);

    float _findDir( std::vector<std::vector<float>> chargeMap,
                    int vtx_col=255,
                    int vtx_row=255,
                    float scanLen = 50,
                    float scanOpen=0.05 );

    float _findLen( std::vector<std::vector<float>> chargeMap,
                    float theta,
                    int vtx_col=255, int vtx_row=255,
                    float scanOpen=0.2);

    float _findOpen( std::vector<std::vector<float>> chargeMap,
                     float theta,
                     float length,
                     int vtx_col=255, int vtx_row=255 );

     std::vector<std::vector<float>> MakeImage2dSparse(const larcv::Image2D& input_img,
              float threshold = 0 );


    // parameters
    // -----------
    std::string _adc_tree_name;
    std::string _calib_adc_tree_name;
    std::string _ssnet_shower_image_stem; // sometimes ubspurn (when files made at FNAL)
    std::string _vertex_tree_name;
    std::string _partroi_tree_name;
    std::string _track_tree_name;
    std::string _mcshower_tree_name;
    std::string _instance_tree_name;
    std::string _segment_tree_name;

    float _Qcut;
    float _SSNET_SHOWER_THRESHOLD;
    bool  _use_calibrated_pixelsum2mev;
    bool _use_ncpi0;
    bool _use_nueint;
    bool _second_shower;
    float _second_shower_adc_threshold;

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
    std::vector< std::vector<float> >  _secondshower_energy_vv;
    std::vector< std::vector<float> >  _secondshower_sumQ_vv;
    std::vector< std::vector<float> >  _secondshower_shlength_vv;
    std::vector< std::vector<double> > _vtx_pos_vv;
    std::vector< larlite::shower >         _shower_ll_v;
    std::vector< larlite::larflowcluster > _shower_pixcluster_v;
    std::vector<double>  _true_energy_vv;
    std::vector<std::vector<double>> _true_shower_start_vv;
    double _remaining_adc;
    std::vector<std::vector<float>> _match_y1_vv;
    std::vector<std::vector<float>> _match_y2_vv;
    std::vector<std::vector<int>> _bestmatch_y1_vv;
    std::vector<std::vector<int>> _bestmatch_y2_vv;
    std::vector<std::vector<float>> _uplane_profile_vv;
    std::vector<std::vector<float>> _vplane_profile_vv;
    std::vector<std::vector<float>> _yplane_profile_vv;
    std::vector<float> _ShowerTruthMatch_v;
  public:
    void clear();


  };


}
}

#endif
