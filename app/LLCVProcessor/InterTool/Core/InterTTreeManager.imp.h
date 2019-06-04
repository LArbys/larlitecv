#ifndef __INTERTTREEMANAGER_IMP_H__
#define __INTERTTREEMANAGER_IMP_H__

#include "InterTTreeManager.h"

namespace llcv {

  template <class T> T InterTTreeManager::Scalar(const std::string& str) const {

    static std::unordered_map<std::string, int>::const_iterator search_int;
    static std::unordered_map<std::string, float>::const_iterator search_float;
    static std::unordered_map<std::string, double>::const_iterator search_double;

    search_int = _spec._imap.find(str);
    if(search_int != _spec._imap.end())
      return (T)(search_int->second);

    search_float = _spec._fmap.find(str);
    if(search_float != _spec._fmap.end())
      return (T)(search_float->second);

    search_double = _spec._dmap.find(str);
    if(search_double != _spec._dmap.end())
      return (T)(search_double->second);
    
    LLCV_CRITICAL() << "Could not find scalar " << str << " in TTree" << std::endl;
    throw llcv_err("die");

    T res;      
    return res;
  }


}

#endif
