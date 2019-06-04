#ifndef TRACKPGRAPHMATCH_H
#define TRACKPGRAPHMATCH_H

#include "LLCVBase/AnaBase.h"

namespace llcv {
  
  class TrackPGraphMatch : public AnaBase {
    
  public:
    
  TrackPGraphMatch(const std::string name="TrackPGraphMatch") : AnaBase(name) {}
    ~TrackPGraphMatch() {}

    void configure(const larcv::PSet& cfg);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();

  private:

    std::string _trk_reco_prod;
    std::string _adc_img_prod; 
    std::string _pgraph_prod;
    std::string _pixel_prod;

    TTree* _tree;
    
    int _run;
    int _subrun;
    int _event;
    int _entry;

    int _vtxid;

    float _vtx_x;
    float _vtx_y;
    float _vtx_z;

    int _ntracks;
    int _npars;

    std::vector<int> _npx_v;
    std::vector<int> _npts_v;
    std::vector<int> _trk_type_v;
    std::vector<std::vector<int> > _trk_type_vv;

    void ClearVertex();
     
  };

}

#endif
