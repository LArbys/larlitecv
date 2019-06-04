#ifndef __SELMICHELFLASH_H__
#define __SELMICHELFLASH_H__

#include "InterTool_Core/InterSelBase.h"

namespace llcv {
  
  class SelMichelFlash : public InterSelBase { 

  public:

  SelMichelFlash(std::string name="SelMichelFlash") : 
    InterSelBase(name),
      _outtree(nullptr) {}
    
    ~SelMichelFlash(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
    
  protected:

    TTree* _outtree;

    int _ticks_from_flash;
    int _ticks_in_fit;
    int _baseline_estimate;
    int _nclose_pmt;

    std::vector<std::vector<short> > _opdigit_vv;
    std::vector<std::vector<float> > _opch_to_xyz;

    double _michel_x;
    double _michel_y;
    double _michel_z;
    
    std::vector<int> _close_pmt_v;

    int _start_tick;
    int _end_tick; 
    float _amplitude;
    float _baseline;  
    float _fit_constant; 
    float _fit_chi2;
    float _fit_Ndf;

    double Distance2D(double x1, double x2,double y1, double y2);
    std::vector<size_t> argmin(const std::vector<float> &v);
    std::vector<size_t> argmax(const std::vector<float> &v);
    
    void Reset();

  };

}


#endif
