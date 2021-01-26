#ifndef __LLCV_GET_TRUTH_VARIABLES_H__
#define __LLCV_GET_TRUTH_VARIABLES_H__

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
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"

//DQDXCalculator Package
#include "DQDXCalculator/UtilFunctions.h"
#include "DQDXCalculator/Visualize_Functions.h"
#include "DQDXCalculator/DQDXBuilder.h"

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

namespace larlitecv {
namespace ssnetshowerreco {

  class GetTruthVariables {

  public:

    GetTruthVariables();
    virtual ~GetTruthVariables() {};

    void initialize();
    bool process(larcv::IOManager& iolcv,larlite::storage_manager& ioll, int entry);
    void finalize();
    // get functions
    // -------------
    int NumParts() const { return _numparts; };
    double GetNeutrinoEnergy() const { return _enu_true; };
    float getSSNetShowerAverage() const { return _averageSSNetShowerScore; };
    int getParticlePDG(int i) const { return _particle_PDG_v[i]; };
    int getParticleStatus(int i) const { return _particle_status_v[i]; };
    double getPartMomentumX(int i) const { return _particle_momx_v[i]; };
    double getPartMomentumY(int i) const { return _particle_momy_v[i]; };
    double getPartMomentumZ(int i) const { return _particle_momz_v[i]; };
    double getPartE(int i) const { return _particle_E_v[i]; };

  protected:

    // parameters
    // -----------
    std::string _mctruth_name;
    std::string _mcshower_tree_name;
    std::string _ssnet_shower_image_stem; // sometimes ubspurn (when files made at FNAL)
    std::string _thrumu_tree_name;


  public:



  protected:
    // variables filled
    // -----------------
    int  _run;
    int  _subrun;
    int  _event;
    int _numparts;
    double _enu_true;
    float _averageSSNetShowerScore;
    std::vector<int> _particle_PDG_v;
    std::vector<int> _particle_status_v;
    std::vector<double> _particle_momx_v;
    std::vector<double> _particle_momy_v;
    std::vector<double> _particle_momz_v;
    std::vector<double> _particle_E_v;
    std::vector<double> _gamma_energy_v;
  public:
    void clear();
  };


}
}

#endif
