#ifndef __MATCHOBJECTALGOTIMEIOU_H__
#define __MATCHOBJECTALGOTIMEIOU_H__

#include "MatchObjectAlgoBase.h"

namespace llcv {

  class MatchObjectAlgoTimeIOU : public MatchObjectAlgoBase {

  public:
    
  MatchObjectAlgoTimeIOU() : MatchObjectAlgoBase("MathAlgoTimeIOU") {}
    ~MatchObjectAlgoTimeIOU() {}
    
    void Configure(const larcv::PSet& pset);
    void Initialize();
    
  protected:

    float Match(const Object2D& obj1, const Object2D& obj2);
    float Match(const Object2D& obj1, const Object2D& obj2, const Object2D& obj3);
    
  };
  
}

#endif
