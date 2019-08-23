#ifndef __SELMICHELFLASH_CXX__
#define __SELMICHELFLASH_CXX__

#include "SelMichelFlash.h"
#include "LArUtil/Geometry.h"

#include <cmath>
#include <sstream>
#include <numeric>
#include <iomanip>

#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include "TFitResult.h"

namespace llcv {

  void SelMichelFlash::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    
    _ticks_from_flash  = pset.get<int>("TicksFromFlash");
    _ticks_in_fit      = pset.get<int>("TicksInFit");
    _baseline_estimate = pset.get<int>("TicksForBaseline");
    _nclose_pmt        = pset.get<int>("NClosePMT");

    auto geo = larutil::Geometry::GetME();    

    _opch_to_xyz.resize(32);
    std::vector<double> xyz(3,-1);
    for(int opch=0; opch<32; ++opch) {
      auto opid = geo->OpDetFromOpChannel(opch);
      geo->GetOpDetPosition(opid,xyz);
      _opch_to_xyz[opch].resize(3,-1);
      _opch_to_xyz[opch][0] = xyz[0];
      _opch_to_xyz[opch][1] = xyz[1];
      _opch_to_xyz[opch][2] = xyz[2];
    }

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelMichelFlash::Initialize() {

    _outtree = new TTree("michel","");
    AttachRSEV(_outtree);
    _outtree->Branch("opdigit_vv",&_opdigit_vv);
    _outtree->Branch("michel_x",&_michel_x,"michel_x/D");
    _outtree->Branch("michel_y",&_michel_y,"michel_y/D");
    _outtree->Branch("michel_z",&_michel_z,"michel_z/D");

    _outtree->Branch("close_pmt_v" , &_close_pmt_v);

    _outtree->Branch("start_tick"  , &_start_tick    , "start_tick/I");
    _outtree->Branch("end_tick"    , &_end_tick      , "end_tick/I");
    _outtree->Branch("baseline"    , &_baseline      , "baseline/F");
    _outtree->Branch("amplitude"   , &_amplitude     , "amplitude/F");
    _outtree->Branch("fit_constant", &_fit_constant  , "fit_constant/F");
    _outtree->Branch("fit_chi2"    , &_fit_chi2      , "fit_chi2/F");      
    _outtree->Branch("fit_Ndf"     , &_fit_Ndf       , "fit_Ndf/F");

  }
  
  double SelMichelFlash::Select() {
    Reset();

    LLCV_DEBUG() << "start" << std::endl;

    LLCV_DEBUG() << "Is michel?: " << Tree().Scalar<int>("ismichel") << std::endl;
    if (Tree().Scalar<int>("ismichel") == 0) { 
      LLCV_DEBUG() << "SKIP!" << std::endl;
      return 0.0;
    }

    _michel_x = Data().Vertex()->X();
    _michel_y = Data().Vertex()->Y();
    _michel_z = Data().Vertex()->Z();

    const auto& digits = Data().Digits();

    LLCV_DEBUG() << "GOT: " << digits.size() << std::endl;

    std::vector<float> ch_distance_v(32,kINVALID_FLOAT);

    for(const auto digit : digits) {
      if (digit->ChannelNumber() > 32) continue;
      size_t id = digit->ChannelNumber();
      _opdigit_vv[id] = (std::vector<short>)(*digit);
      ch_distance_v[id] = Distance2D(_michel_y,_opch_to_xyz[id][1],
				     _michel_z,_opch_to_xyz[id][2]);
    }
    
    //
    // find the closest PMTs, sum their waveforms
    //
    auto dist_id_v = argmin(ch_distance_v);
    std::vector<float> wf_sum_v(_opdigit_vv.front().size(),0.0);

    for(int id=0; id < _nclose_pmt; ++id) {
      auto dist_id = dist_id_v[id];
      LLCV_DEBUG() << "@id=" << id << " opch=" << dist_id << std::endl;

      _close_pmt_v[id] = dist_id;
      for(size_t did=0; did<wf_sum_v.size(); ++did)
	wf_sum_v[did] += _opdigit_vv[dist_id][did];
    }


    //
    // estimate the baseline
    //
    float baseline = 0.0;
    for(int id=0; id<_baseline_estimate; ++id)
      baseline += wf_sum_v[id];

    baseline /= (float)_baseline_estimate;
    std::cout << "baseline=" << baseline << std::endl;

    //
    // get the max value and location
    //
    auto wf_max_iter = std::max_element(wf_sum_v.begin(), wf_sum_v.end());
    size_t wf_max_id = std::distance(wf_sum_v.begin(), wf_max_iter);
    
    //
    // set the start and end for fit
    //
    size_t start_x = wf_max_id + _ticks_from_flash;
    size_t end_x   = start_x + _ticks_in_fit;
    float start_y = wf_sum_v.at(start_x);

    //
    // setup the fit function and fill data
    //
    std::stringstream ss;
    ss << std::setprecision(10);
    ss << baseline << "+(" << start_y << "-" << baseline << ")*";
    ss << "TMath::Exp([0]*(x-" << start_x << "))";
    
    std::vector<float> ydata_fit;
    std::vector<float> xdata_fit;
    for(size_t id=start_x; id<end_x; ++id) {
      xdata_fit.push_back((float)id);
      ydata_fit.push_back(wf_sum_v[id]);
    }

    TGraph tg((Int_t)xdata_fit.size(),xdata_fit.data(),ydata_fit.data());
    TF1 fit("expo_fit",ss.str().c_str(),xdata_fit.front(),xdata_fit.back());

    fit.SetParameter(0,-1.0);
    TFitResultPtr fit_res_ptr = tg.Fit(&fit,"SQN0");
    
    auto fit_constant = fit_res_ptr->Value(0);
    auto fit_chi2     = fit_res_ptr->Chi2();
    auto fit_Ndf      = fit_res_ptr->Ndf();

    _start_tick   = (int)   start_x;
    _end_tick     = (int)   end_x;
    _baseline     = (float) baseline;
    _amplitude    = (float) start_y;
    _fit_constant = (float) fit_constant;
    _fit_chi2     = (float) fit_chi2;
    _fit_Ndf      = (float) fit_Ndf;

    _outtree->Fill();
    return 0.0;
  }
  
  double SelMichelFlash::Distance2D(double x1, double x2,
				    double y1, double y2) {
    return std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  }


  std::vector<size_t> SelMichelFlash::argmin(const std::vector<float> &v) {
    
    std::vector<size_t> idx_v(v.size());
    std::iota(idx_v.begin(), idx_v.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(idx_v.begin(), idx_v.end(),
	      [&v](size_t i1, size_t i2) {
		return v[i1] < v[i2];
	      });
    
    return idx_v;
  }

  std::vector<size_t> SelMichelFlash::argmax(const std::vector<float> &v) {
    
    std::vector<size_t> idx_v(v.size());
    std::iota(idx_v.begin(), idx_v.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(idx_v.begin(), idx_v.end(),
	      [&v](size_t i1, size_t i2) {
		return v[i1] > v[i2];
	      });
    
    return idx_v;
  }
  
  void SelMichelFlash::Finalize() {
    _outtree->Write();
  }


  void SelMichelFlash::Reset() {
    _opdigit_vv.clear();
    _opdigit_vv.resize(32);

    _close_pmt_v.clear();
    _close_pmt_v.resize(_nclose_pmt,-1);

    _michel_x = kINVALID_DOUBLE;
    _michel_y = kINVALID_DOUBLE;
    _michel_z = kINVALID_DOUBLE;

    return;
  }

}


#endif
