#ifndef TRACKTRUTHMATCH_H
#define TRACKTRUTHMATCH_H

#include "LLCVBase/AnaBase.h"

namespace llcv {
  
  class TrackTruthMatch : public AnaBase {
    
  public:
    
  TrackTruthMatch(const std::string name="TrackTruthMatch") : AnaBase(name) {}
    ~TrackTruthMatch() {}

    void configure(const larcv::PSet& cfg);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();

  private:

    std::string _trk_reco_prod;
    std::string _adc_img_prod; 
    std::string _seg_img_prod; 

    TTree* _tree;
    
    int _run;
    int _subrun;
    int _event;
    int _entry;
    int _vtxid;
    int _ntracks;

    std::vector<float> _unknownfrac_v; 
    std::vector<float> _electronfrac_v;
    std::vector<float> _gammafrac_v;   
    std::vector<float> _pizerofrac_v;  
    std::vector<float> _muonfrac_v;    
    std::vector<float> _kminusfrac_v;  
    std::vector<float> _piminusfrac_v; 
    std::vector<float> _protonfrac_v;  

    std::vector<int> _npx_v;
    std::vector<int> _npts_v;
    std::vector<int> _trk_type_v;
    std::vector<std::vector<int> > _trk_type_vv;

    void ClearVertex();
      
  };

}

#endif
