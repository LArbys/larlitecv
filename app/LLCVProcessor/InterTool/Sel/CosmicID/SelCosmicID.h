#ifndef __SELCOSMICID_H__
#define __SELCOSMICID_H__

#include "InterTool_Core/InterSelBase.h"

namespace llcv {
  
  class SelCosmicID : public InterSelBase { 

  public:

  SelCosmicID(std::string name="SelCosmicID") 
    :   InterSelBase(name)
      , _out_tree(nullptr) 
    {}

    ~SelCosmicID(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();

  private:
    void ResetTree();
    
    
  private:
    
    TTree* _out_tree;

    float _max_radius;
    float _min_radius;
    float _radius_step;

    std::vector<float> _radius_v;
    std::vector<std::vector<int> > _pt_xing_vv;
    std::vector<std::vector<int> > _connected_vv;
    std::vector<std::vector<float> > _distance_vv;

  };

}


#endif
