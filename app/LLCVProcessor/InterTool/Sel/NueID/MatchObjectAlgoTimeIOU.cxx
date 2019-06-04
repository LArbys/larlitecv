#ifndef __MATCHOBJECTALGOTIMEIOU_CXX__
#define __MATCHOBJECTALGOTIMEIOU_CXX__

#include "MatchObjectAlgoTimeIOU.h"

namespace llcv {
  
  void MatchObjectAlgoTimeIOU::Configure(const larcv::PSet& pset) {
    this->set_verbosity((msg::Level_t)(pset.get<unsigned short>("Verbosity", (unsigned short)(this->logger().level()))));
    // this->set_verbosity((msg::Level_t)0);
    _threshold = pset.get<float>("Threshold",0.1);
    _require_plane2 = pset.get<bool>("RequirePlane2",false);
    _match_three_planes = pset.get<bool>  ("MatchThreePlanes",true);
    _three_planes_boost = pset.get<float> ("ThreePlanesBoost",1.1);
    _plane_two_boost    = pset.get<float> ("PlaneTwoBoost", 1.0);
  }

    
  void MatchObjectAlgoTimeIOU::Initialize() {}
  
  float MatchObjectAlgoTimeIOU::Match(const Object2D& obj1, const Object2D& obj2) {

    float score = -1.0 * larocv::kINVALID_FLOAT;

    const auto& ctor0 = obj1._line;
    const auto& ctor1 = obj2._line;

    //min @ max time tick
    int min_time_tick0 = larocv::kINVALID_INT;
    int max_time_tick0 = -1.0*larocv::kINVALID_INT;

    for(const auto& px : ctor0) {
      min_time_tick0 = std::min(min_time_tick0,px.x);
      max_time_tick0 = std::max(max_time_tick0,px.x);
    }

    int min_time_tick1 = larocv::kINVALID_INT;
    int max_time_tick1 = -1.0*larocv::kINVALID_INT;

    for(const auto& px : ctor1) {
      min_time_tick1 = std::min(min_time_tick1,px.x);
      max_time_tick1 = std::max(max_time_tick1,px.x);
    }

    auto min_time0 = min_time_tick0;
    auto min_time1 = min_time_tick1;

    auto max_time0 = max_time_tick0;
    auto max_time1 = max_time_tick1;

    // which one is higher in time?
    float max_time_top = larocv::kINVALID_FLOAT;
    float min_time_top = larocv::kINVALID_FLOAT;

    float max_time_bot = larocv::kINVALID_FLOAT;
    float min_time_bot = larocv::kINVALID_FLOAT;

    if (max_time0 > max_time1) {
      max_time_top = max_time0;
      min_time_top = min_time0;

      max_time_bot = max_time1;
      min_time_bot = min_time1;
    } else {
      max_time_top = max_time1;
      min_time_top = min_time1;

      max_time_bot = max_time0;
      min_time_bot = min_time0;
    }
    
    LLCV_DEBUG() << "top (min,max)=(" << min_time_top << "," << max_time_top << ")" << std::endl;
    LLCV_DEBUG() << "bot (min,max)=(" << min_time_bot << "," << max_time_bot << ")" << std::endl;

    if (min_time_top > max_time_bot) return score;

    auto max_time = max_time_top;
    auto min_time = std::min(min_time_top,min_time_bot);

    auto com_max_time = max_time_bot;
    auto com_min_time = std::max(min_time_top,min_time_bot);
    
    float timeiou = (com_max_time - com_min_time)  / (max_time - min_time);
    
    score = timeiou;
    
    return score;
  }
    
  float MatchObjectAlgoTimeIOU::Match(const Object2D& obj1, const Object2D& obj2, const Object2D& obj3) {
    
    float score = -1.0 * kINVALID_FLOAT;
    
    auto m01  = Match(obj1,obj2);
    auto m02  = Match(obj1,obj3);
    auto m12  = Match(obj2,obj3);
    
    LLCV_DEBUG() << "(m00,m01,m12)=("<<m01<<","<<m02<<","<<m12<<")" << std::endl;
    
    auto m_max = std::max({m01, m02, m12});
    auto m_min = std::min({m01, m02, m12});
    
    // if ((m_max - m_min) > 0.3)
    if ((m_max - m_min) > 0.1)
      return score;
  
    score = m_max;

    return score;
  }

}
#endif
