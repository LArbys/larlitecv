#ifndef __SELMICHELSTUDY_H__
#define __SELMICHELSTUDY_H__

#include "InterTool_Core/InterSelBase.h"

namespace llcv {
  
  class SelMichelStudy : public InterSelBase { 

  public:

  SelMichelStudy(std::string name="SelMichelStudy") : InterSelBase(name) {}
    ~SelMichelStudy(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize() {}
    double Select();
    void Finalize();
    
    
  protected:
    
    
  };

}


#endif
