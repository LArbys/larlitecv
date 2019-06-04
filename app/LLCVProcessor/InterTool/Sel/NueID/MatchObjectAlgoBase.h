#ifndef __MATCHOBJECTALGOBASE_H__
#define __MATCHOBJECTALGOBASE_H__

#include "LLCVBase/llcv_base.h"
#include "Base/PSet.h"
#include "Object2D.h"
#include "LArOpenCV/ImageCluster/Base/BaseUtil.h"
#include "LArOpenCV/ImageCluster/Base/MatchBookKeeper.h"
#include <vector>
#include <array>

namespace llcv {

  class MatchObjectAlgoBase : public llcv_base {

  public:
    
  MatchObjectAlgoBase() : llcv_base("MatchObjectAlgoBase") {}
  MatchObjectAlgoBase(const std::string& name) : llcv_base(name) {}
    
    virtual ~MatchObjectAlgoBase() {}
    
    virtual void Configure(const larcv::PSet &pset) = 0;
    virtual void Initialize() = 0;
    
    void Register(const Object2D& obj, size_t plane);
    std::vector<std::vector<std::pair<size_t,size_t> > > MatchObjects(std::vector<float>& score_v);
    const Object2D* Object(size_t plane,size_t id);

    void ClearMatch();
    void ClearEvent();

  protected:
    virtual float Match(const Object2D& obj1, const Object2D& obj2) = 0;
    virtual float Match(const Object2D& obj1, const Object2D& obj2, const Object2D& obj3) = 0;
    
    bool  _match_three_planes;
    bool  _require_plane2;
    float _three_planes_boost;
    float _plane_two_boost;
    float _threshold;
    
  private:
    
    larocv::MatchBookKeeper _MatchBookKeeper;
    std::vector<size_t> _seed_v;

    std::array<std::vector<size_t>,3> _object_per_plane_vv;
    std::vector<const Object2D*> _object_ptr_v;
    std::vector<std::pair<size_t,size_t> > _object_id_to_plane_v;

    std::vector<unsigned int> _match_v;
    std::array<std::vector<const Object2D*>,3> _object_vv;
    std::array<bool,3> _valid_plane_v;

  };

}

#endif
