#ifndef __INTERTOOLTYPES_CXX__
#define __INTERTOOLTYPES_CXX__

#include "InterToolTypes.h"

namespace llcv {


  InterSpecType LeafToSpecType(const std::string& type) {
    if (type=="Int_t"   ) return kINT;
    if (type=="Float_t" ) return kFLOAT;
    if (type=="Double_t") return kDOUBLE;
    if (type=="vector<float>")          return kVFLOAT;
    if (type=="vector<vector<float> >") return kVVFLOAT;
    return kINTER_SPEC_TYPE_UNKNOWN;
  };


}



#endif
