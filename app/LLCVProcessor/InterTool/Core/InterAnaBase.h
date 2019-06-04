#ifndef __INTERANABASE_H__
#define __INTERANABASE_H__

#include "LLCVBase/llcv_base.h"

namespace llcv {

  class InterAnaBase : public llcv_base { 
  public:
  InterAnaBase(std::string name="InterAnaBase") : llcv_base(name) {}
    virtual ~InterAnaBase(){}
    
  };

}


#endif
