#ifndef __INTERTTREESPEC_H__
#define __INTERTTREESPEC_H__

#include "LLCVBase/llcv_err.h"
#include "TChain.h"
#include <unordered_map>

namespace llcv {

  class InterTTreeSpec {
  public:

    void LoadRSEV(TChain& tc);
    void LoadTree(TChain& tc);

    std::unordered_map<std::string, int> _imap;
    std::unordered_map<std::string, double> _dmap;
    std::unordered_map<std::string, float> _fmap;
    std::unordered_map<std::string, std::vector<float>* > _vmap;
    std::unordered_map<std::string, std::vector<std::vector<float> >* > _vvmap;

    int* _run_ptr;
    int* _subrun_ptr;
    int* _event_ptr;
    double* _vtxid_ptr;
    
  private:

    int _run;
    int _subrun;
    int _event;
    double _vtxid;
    
  };

}


#endif 
