#ifndef __LLCV_SSNET_SHOWER_RECO_H__
#define __LLCV_SSNET_SHOWER_RECO_H__

#include <vector>
#include <string>
#include <cstring>

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

    SSNetShowerReco(bool make_larlite,std::string outputname);
    virtual ~SSNetShowerReco() {};

    void initialize();
    bool process( larcv::IOManager& iocv, larlite::storage_manager& ioll, int entry);//, larcv::IOManager& ioimgs );
    void setupAnaTree();
    void finalize();
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
    float getOverlapFraction1(int vtx, int plane, int shower ) const { return _match_y1_vv[vtx][plane][shower]; };
    float getOverlapFraction2(int vtx, int plane, int shower ) const { return _match_y2_vv[vtx][plane][shower]; };
    int numpointsU() const {return _uplane_profile_vv.size(); };
    int numpointsV() const {return _vplane_profile_vv.size(); };
    int numpointsY() const {return _yplane_profile_vv.size(); };
    float getUPlaneShowerProfile( int pt, int ii ) const { return _uplane_profile_vv[pt][ii]; };
    float getVPlaneShowerProfile( int pt, int ii ) const { return _vplane_profile_vv[pt][ii]; };
    float getYPlaneShowerProfile( int pt, int ii ) const { return _yplane_profile_vv[pt][ii]; };
    float getShowerTruthMatchPur (int vtx, int shower) const {return _ShowerTruthMatch_pur_vv[vtx][shower]; };
    float getShowerTruthMatchEff (int vtx, int shower) const {return _ShowerTruthMatch_eff_vv[vtx][shower]; };
    int getUseForMass(int vtx) const {return _useformass[vtx];};
    float getPi0Mass(int vtx) const {return _pi0mass[vtx];};
    float getDistToInt(int vtx) const {return _disttoint[vtx];};
    float getImpact1(int vtx) const {return _impact1[vtx];};
    float getImpact2(int vtx) const {return _impact2[vtx];};
    int getHasPi0() const {return _haspi0;};
    int getCCNC() const {return _ccnc;};
    float getAlpha(int vtx) const {return _alpha[vtx];};
    float getFirstDirection(int vtx, int dir) const {return _firstdirection[vtx][dir];};
    float getSecondDirection(int vtx, int dir) const {return _seconddirection[vtx][dir];};

    larlite::shower& getShowerObject( int vtxid, int plane ) { return _shower_ll_v.at( 3*vtxid+plane ); };
    larlite::larflowcluster& getShowerPixelList( int vtxid, int plane ) { return _shower_pixcluster_v.at( 3*vtxid+plane ); };

    void store_in_larlite( larlite::storage_manager& ioll );
    void use_calibrated_pixsum2mev( bool use=true ) { _use_calibrated_pixelsum2mev = use; };
    void use_second_shower( bool use=true ) { _second_shower = use; };
    void use_ncpi0( bool use=true ) { _use_ncpi0 = use; };
    void use_nueint( bool use=true ) { _use_nueint = use; };
    void use_bnb( bool use=true ) { _use_bnb = use; };


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
                                                     bool calcE,
                                                     int vtx_col=255,
                                                     int vtx_row=255,
                                                     float shLen = 100.0,
                                                     float shOpen = 0.3);

    float _findDir( std::vector<std::vector<float>> chargeMap,
                    int vtx_col=255,
                    int vtx_row=255,
                    float scanLen = 200,
                    float scanOpen=0.1 );

    float _findLen( std::vector<std::vector<float>> chargeMap,
                    float theta,
                    int vtx_col=255, int vtx_row=255,
                    float scanOpen=0.3);

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
    std::string _ssnet_shower_tree_name;  // default: ssnetshowerreco
    std::string _vertex_tree_name;
    std::string _partroi_tree_name;
    std::string _track_tree_name;
    std::string _mcshower_tree_name;
    std::string _instance_tree_name;
    std::string _segment_tree_name;
    std::string _mctruth_name;
    std::string _thrumu_tree_name;

    bool  _make_larlite;
    float _Qcut;
    float _SSNET_SHOWER_THRESHOLD;
    bool  _use_calibrated_pixelsum2mev;
    bool _use_ncpi0;
    bool _use_nueint;
    bool _use_bnb;
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
    void set_output_treename( std::string name )       { _ssnet_shower_tree_name = name; };


  protected:
    // variables filled
    // -----------------
    int _run;
    int _subrun;
    int _event;
    std::vector<int> _vtxid;
    std::vector< std::vector<float> >  _shower_energy_vv;
    std::vector< std::vector<int> >  _shower_gap_vv;
    std::vector< std::vector< std::vector<int> > > _shower_start_2d_vvv;
    std::vector< std::vector<float> >  _shower_sumQ_vv;
    std::vector< std::vector<float> >  _shower_shlength_vv;
    std::vector< std::vector<float> >  _shower_shangle_vv;
    std::vector< std::vector<float> >  _shower_shopen_vv;
    std::vector< std::vector<float> >  _secondshower_energy_vv;
    std::vector< std::vector<int> >  _secondshower_gap_vv;
    std::vector< std::vector< std::vector<int> > > _secondshower_start_2d_vvv;
    std::vector< std::vector<float> >  _secondshower_sumQ_vv;
    std::vector< std::vector<float> >  _secondshower_shlength_vv;
    std::vector< std::vector<float> >  _secondshower_shangle_vv;
    std::vector< std::vector<float> >  _secondshower_shopen_vv;
    std::vector< std::vector<double> > _vtx_pos_vv;
    std::vector< larlite::shower >         _shower_ll_v;
    std::vector< larlite::shower >         _secondshower_ll_v;
    std::vector< larlite::larflowcluster > _shower_pixcluster_v;
    std::vector<double>  _true_energy_vv;
    std::vector<std::vector<double>> _true_shower_start_vv;
    double _remaining_adc;
    std::vector<std::vector<std::vector<float>>> _match_y1_vv;
    std::vector<std::vector<std::vector<float>>> _match_y2_vv;
    std::vector<std::vector<std::vector<int>>> _bestmatch_y1_vv;
    std::vector<std::vector<std::vector<int>>> _bestmatch_y2_vv;
    std::vector<std::vector<float>> _uplane_profile_vv;
    std::vector<std::vector<float>> _vplane_profile_vv;
    std::vector<std::vector<float>> _yplane_profile_vv;
    std::vector<std::vector<float>> _ShowerTruthMatch_pur_vv;
    std::vector<std::vector<float>> _ShowerTruthMatch_eff_vv;
    std::vector<int> _useformass;
    std::vector<float> _pi0mass;
    std::vector<float> _disttoint;
    std::vector<float> _impact1;
    std::vector<float> _impact2;
    std::vector<float> _alpha;
    std::vector<std::vector<float>> _seconddirection;
    std::vector<std::vector<float>> _firstdirection;
    int _haspi0;
    int _ccnc;
    TTree* _ana_tree;
    TFile* OutFile;
  public:
    void clear();
    TTree* getAnaTree() { return _ana_tree; };

  };


}
}

#endif
