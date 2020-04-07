#ifndef __LLCV_SECOND_SHOWER_H__
#define __LLCV_SECOND_SHOWER_H__

#include <vector>
#include <string>

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/IOManager.h"

// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/shower.h"
#include "DataFormat/larflowcluster.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcnu.h"
#include "DataFormat/track.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/potsummary.h"

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
#include "TPrincipal.h"
#include "TMatrixD.h"

#include "Utils.h"

namespace larlitecv {
namespace ssnetshowerreco {

  class SecondShower {
    //set of functions to run a second shower search with the shower reco

  public:

    typedef std::vector< std::vector<float> > triangle_t;

    //functions
    std::vector< std::vector<int> > _enclosedCharge( std::vector<std::vector<float>> chargeMap,
                                                     float theta,
                                                     float& sumIn,
                                                     std::vector< std::vector<float> >& triangle,
                                                     bool calcE,
                                                     int vtx_col=255,
                                                     int vtx_row=255,
                                                     float shLen = 100.0,
                                                     float shOpen = 0.3);
     bool _isInside2( float x1, float y1,float x2, float y2, float x3, float y3,
              float x, float y );
    float _sign( float x1, float y1,float x2, float y2,float x3, float y3 );
    int ChooseGapSize(int vtx_col,int vtx_row,float shangle,
      std::vector<std::vector<float>> crop_v, int showernum);

    //truth functions
    std::vector<double> SaveTrueEnergies(larlite::event_mcshower* ev_mcshower);
    std::vector<std::vector<double>> SaveTrueStarts(larlite::event_mcshower* ev_mcshower);
    std::vector<std::vector<double>> SaveTrueDirections(larlite::event_mcshower* ev_mcshower);
    std::vector<float> RecoTrueDistances(std::vector<int> true_shower1_tmp_2dstart,
          std::vector<int> true_shower2_tmp_2dstart, std::vector<int>shower_start_2d_v);

    std::vector<std::vector<float>> Match_3D(std::vector<std::vector<float>> triangle_vv,
        std::vector<std::vector<std::vector<float>>> img_vvv, int shower_num);
    std::vector<std::vector<int>> ChooseBestMatch(std::vector<std::vector<float>> match_vv);
    std::vector<std::vector<float>> SaveTrueProfile(int plane, larcv::EventImage2D* ev_segment,
          larcv::EventImage2D* ev_instance,larlite::event_mcshower* ev_mcshower,
          std::vector<std::vector<std::vector<float>>> img_vvv);
    std::vector<float> TruthMatchNCPi0(std::vector<float> triangle_v, std::vector<std::vector<float>> img_vv,
          larcv::EventImage2D* ev_segment, larcv::EventImage2D* ev_instance, int plane);
    std::vector<float> TruthMatchNueint(std::vector<float> triangle_v, std::vector<std::vector<float>> img_vv,
          larcv::EventImage2D* ev_segment, larcv::EventImage2D* ev_instance, int plane);
    bool RunMassCalc(std::vector<std::vector<float>> _match_y1_vv,
          std::vector<std::vector<float>>_match_y2_vv,
          std::vector<float> shower_energy_v,
          std::vector<float> secondshower_energy_v);
    std::vector<std::vector<float>> Get3DPoints(std::vector<float> shower_points1,
            std::vector<float> shower_points2, std::vector<std::vector<std::vector<float>>> masked_adc_vvv,
            int planeofmatch, larcv::ImageMeta wire_meta);
    std::vector<float> GetOpeningAngle(std::vector<std::vector<float>> firstshower, std::vector<std::vector<float>> secondshower,
            std::vector<double> vertex);
    std::vector<std::vector<float>> GetPCA(std::vector<std::vector<float>> shower);


  protected:


    // parameters
    // -----------



  };//end of second shower class


}//end of ssnetshowerreco namespace
}//end of larlitecv namespace

#endif
