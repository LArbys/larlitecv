#ifndef __INTERPMT_H__
#define __INTERPMT_H__

//llcv
#include "LLCVBase/AnaBase.h"

#include "InterTool_Core/InterDriver.h"

namespace llcv {
  
  class InterPMT : public AnaBase {
    
  public:
    
  InterPMT(const std::string name="InterPMT") :  AnaBase(name) {}
    ~InterPMT() {}
    
    void configure(const larcv::PSet& cfg);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();

    InterDriver& Driver() { return _driver; }

  private:
    
    std::string _adc_img_prod;
    std::string _trk_img_prod;
    std::string _shr_img_prod;

    std::string _track_vertex_prod;
    std::string _opdigit_prod;

    float _epsilon;

    InterDriver _driver;

  };

}

#endif
