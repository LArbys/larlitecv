#ifndef __SELEXAMPLE_H__
#define __SELEXAMPLE_H__

#include "InterTool_Core/InterSelBase.h"

namespace llcv {
  
  class SelExample : public InterSelBase { 

  public:

  SelExample(std::string name="SelExample")
    :   InterSelBase(name)
      , _outtree(nullptr) {}
    
    ~SelExample(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
  private:
    TTree* _outtree;
    
  };

}


#endif
