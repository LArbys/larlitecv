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
                                                     int vtx_col=255,
                                                     int vtx_row=255,
                                                     float shLen = 100.0,
                                                     float shOpen = 0.3);
     bool _isInside2( float x1, float y1,float x2, float y2, float x3, float y3,
              float x, float y );
    float _sign( float x1, float y1,float x2, float y2,float x3, float y3 );
    int ChooseGapSize(int vtx_col,int vtx_row,float shangle,
      std::vector<std::vector<float>> crop_v, int showernum);
    std::vector<double> SaveTrueEnergies(larlite::event_mcshower* ev_mcshower,
                std::vector<double> _scex);
    std::vector<std::vector<double>> SaveTrueStarts(larlite::event_mcshower* ev_mcshower);

    std::vector<std::vector<float>> Match_3D(std::vector<std::vector<float>> triangle_vv,
        std::vector<std::vector<std::vector<float>>> img_vvv, int shower_num);
    std::vector<std::vector<int>> ChooseBestMatch(std::vector<std::vector<float>> match_vv);
    std::vector<std::vector<float>> SaveTrueProfile(int plane, larcv::EventImage2D* ev_segment,
          larcv::EventImage2D* ev_instance,larlite::event_mcshower* ev_mcshower,
          std::vector<std::vector<std::vector<float>>> img_vvv);
    float TruthMatchNCPi0(std::vector<float> triangle_v, std::vector<std::vector<float>> img_vv, 
          larcv::EventImage2D* ev_segment, larcv::EventImage2D* ev_instance, int plane);



  protected:


    // parameters
    // -----------



  };//end of second shower class


}//end of ssnetshowerreco namespace
}//end of larlitecv namespace

#endif
