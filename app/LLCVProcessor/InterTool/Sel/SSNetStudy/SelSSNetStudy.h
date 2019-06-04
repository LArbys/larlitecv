#ifndef __SELSSNETSTUDY_H__
#define __SELSSNETSTUDY_H__

#include "InterTool_Core/InterSelBase.h"

namespace llcv {
  
  class SelSSNetStudy : public InterSelBase { 

  public:

  SelSSNetStudy(std::string name="SelSSNetStudy") : 
    InterSelBase(name), _out_tree(nullptr)  {}
    ~SelSSNetStudy(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
  private:
    void ResetTree();

  private:

    TTree* _out_tree;

    std::vector<std::vector<float> > _ssnet_track_frac_vv;
    std::vector<std::vector<float> > _ssnet_shower_frac_vv;
    std::vector<cv::Mat> _white_mat_v;

  };

}


#endif
