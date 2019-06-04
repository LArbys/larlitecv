#ifndef DEDXCALCULATOR_H
#define DEDXCALCULATOR_H

#include <unordered_map>

namespace llcv {

  class dEdxCalculator {
  public:
    dEdxCalculator();
    ~dEdxCalculator() {}
    float ProtonEstimate(float rr);
    float ProtonEstimateFit(float rr);
    float MuonEstimate(float rr);

  private:
    std::unordered_map<float, float> _proton_m;
    std::unordered_map<float, float> _proton_fit_m;
    std::unordered_map<float, float> _muon_m;
  };


}

#endif
