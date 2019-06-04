#ifndef __INTERMICHEL_H__
#define __INTERMICHEL_H__

//llcv
#include "LLCVBase/AnaBase.h"

#include "InterTool_Core/InterDriver.h"

namespace llcv {
  
  class InterMichel : public AnaBase {
    
  public:
    
  InterMichel(const std::string name="InterMichel") :  AnaBase(name) {}
    ~InterMichel() {}
    
    void configure(const larcv::PSet& cfg);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();

    InterDriver& Driver() { return _driver; }

  private:
    
    std::string _adc_img_prod;
    std::string _trk_img_prod;
    std::string _shr_img_prod;

    std::string _hit_prod;
    // std::string _opflash_prod;

    float _epsilon;

    InterDriver _driver;

  };

}

#endif
