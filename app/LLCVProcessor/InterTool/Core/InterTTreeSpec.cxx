#ifndef __INTERTTREESPEC_CXX__
#define __INTERTTREESPEC_CXX__

#include "InterTTreeSpec.h"
#include "InterToolTypes.h"

#include "TLeafElement.h"

#include <iostream>

namespace llcv {

  void InterTTreeSpec::LoadRSEV(TChain& tc) {
    tc.SetBranchAddress("run"   ,&_run);
    tc.SetBranchAddress("subrun",&_subrun);
    tc.SetBranchAddress("event" ,&_event);
    tc.SetBranchAddress("vtxid" ,&_vtxid);

    _run_ptr    = &_run;
    _subrun_ptr = &_subrun;
    _event_ptr  = &_event;
    _vtxid_ptr  = &_vtxid;
  }

  void InterTTreeSpec::LoadTree(TChain& tc) {

    tc.ResetBranchAddresses();

    auto leaves = tc.GetListOfLeaves();

    for(size_t lid=0; lid<(size_t)leaves->GetEntries(); ++lid) {
      auto leaf = (TLeafElement*)leaves->At(lid);
      auto name  = std::string(leaf->GetName());
      auto type  = std::string(leaf->GetTypeName());
      auto stype = LeafToSpecType(type);
      switch (stype) {
      case kINT  :   tc.SetBranchAddress(name.c_str(),&_imap[name]);  break;
      case kFLOAT:   tc.SetBranchAddress(name.c_str(),&_fmap[name]);  break;
      case kDOUBLE:  tc.SetBranchAddress(name.c_str(),&_dmap[name]);  break;
      case kVFLOAT:  tc.SetBranchAddress(name.c_str(),&_vmap[name]);  break;
      case kVVFLOAT: tc.SetBranchAddress(name.c_str(),&_vvmap[name]); break;
      default: 
	std::cout << "@name=" << name << " @type=" << type << " is unknown" << std::endl;
	throw llcv_err("Unknown type");
      };
    }
    
    _run_ptr    = &_imap.at("run");
    _subrun_ptr = &_imap.at("subrun");
    _event_ptr  = &_imap.at("event");
    _vtxid_ptr  = &_dmap.at("vtxid");
  }




}


#endif 
