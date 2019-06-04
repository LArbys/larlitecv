#ifndef __INTERSELTOY_H__
#define __INTERSELTOY_H__

#include "InterTool_Core/InterSelBase.h"

namespace llcv {
  
  class InterSelToy : public InterSelBase { 

  public:

  InterSelToy(std::string name="InterSelToy") : InterSelBase(name) {}
    ~InterSelToy(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize() {}
    double Select();
    void Finalize();
    
    
  protected:
    
    
  };

}


#endif
