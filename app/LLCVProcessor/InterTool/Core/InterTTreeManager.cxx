#ifndef __INTERTTREEMANAGER_CXX__
#define __INTERTTREEMANAGER_CXX__

#include "InterTTreeManager.h"
#include "InterTTreeManager.imp.h"

namespace llcv {
  
  void InterTTreeManager::Configure(const larcv::PSet& cfg) {
    set_verbosity((msg::Level_t)cfg.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
  }

  void InterTTreeManager::Initialize(const std::string& fname,
				     const std::string& tname) {
    LLCV_DEBUG() << "tname=" << tname << std::endl;
    LLCV_DEBUG() << "fname=" << fname << std::endl;
    _tchain = new TChain(tname.c_str());
    _tchain->Add(fname.c_str());
    _nentries = (int)_tchain->GetEntries();
    _spec.LoadRSEV(*_tchain);

    LLCV_INFO() << "Initializing RSEV map" << std::endl;
    for(size_t ientry=0; ientry<(size_t)_nentries; ++ientry) {
      _tchain->GetEntry(ientry);
      RSEVID rsev(this->Run(),this->SubRun(),this->Event(),this->Vertex());

      LLCV_DEBUG() << "insert=(" 
		   << this->Run()    << "," 
		   << this->SubRun() << "," 
		   << this->Event()  << "," 
		   << this->Vertex() << ") @ centry=" << ientry << std::endl;

      _rsev_m.insert(std::make_pair(rsev,ientry));
    }
    LLCV_INFO() << "initialized" << std::endl;
    _spec.LoadTree(*_tchain);
  }


  bool InterTTreeManager::GoTo(size_t run, size_t subrun, size_t event, size_t vtxid) {
    RSEVID rsev(run,subrun,event,vtxid);
    return _goto(rsev);
  }

  bool InterTTreeManager::_goto(const RSEVID& rsev) {
    LLCV_DEBUG() << "centry=" << _centry << std::endl;

    if(_nentries==-1) {
      LLCV_WARNING() << "No inter file detected" << std::endl;
      return true;
    }

    auto it_rsev = _rsev_m.find(rsev);
    if ( it_rsev==_rsev_m.end() ) {
      throw llcv_err("Could not find RSEV in the inter file");
    }

    _centry = (*it_rsev).second;

    if(_centry>=(size_t)_nentries) 
      throw llcv_err("Out of range centry");

    _tchain->GetEntry(_centry);

    LLCV_DEBUG() << "now centry=" << _centry << std::endl;
    return true;
  }

  template int    InterTTreeManager::Scalar<int   >(const std::string& str) const;
  template float  InterTTreeManager::Scalar<float >(const std::string& str) const;
  template double InterTTreeManager::Scalar<double>(const std::string& str) const;
  
  std::vector<float> InterTTreeManager::Vector(const std::string& str) const {
    
    static std::unordered_map<std::string, std::vector<float>*>::const_iterator search_v;
    
    search_v = _spec._vmap.find(str);
    if(search_v != _spec._vmap.end())
      return *(search_v->second);
    
    LLCV_CRITICAL() << "Could not find vector " << str << " in TTree" << std::endl;
    throw llcv_err("die");
    
    static std::vector<float> res;      
    return res;
    
  }
  
  std::vector<std::vector<float> > InterTTreeManager::VVector(const std::string& str) const {

    static std::unordered_map<std::string, std::vector<std::vector<float> >*>::const_iterator search_vv;
    
    search_vv = _spec._vvmap.find(str);
    if(search_vv != _spec._vvmap.end())
      return *(search_vv->second);
    
    LLCV_CRITICAL() << "Could not find vector vector " << str << " in TTree" << std::endl;
    throw llcv_err("die");
    
    static std::vector<std::vector<float> > res;      
    return res;
    
  }


}


#endif

//hi mom
