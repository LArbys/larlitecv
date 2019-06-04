#ifndef __INTERSELBASE_CXX__
#define __INTERSELBASE_CXX__

#include "InterSelBase.h"

namespace llcv {
  
  void InterSelBase::AttachRSEV(TTree* tree) {
    if(!tree) throw llcv_err("No tree specified");

    if(!_run)    throw llcv_err("No run specified");
    if(!_subrun) throw llcv_err("No subrun specified");
    if(!_event)  throw llcv_err("No event specified");
    if(!_vtxid)  throw llcv_err("No vtxid specified");

    tree->Branch("run"   ,_run   ,"run/I");
    tree->Branch("subrun",_subrun,"subrun/I");
    tree->Branch("event" ,_event ,"event/I");
    tree->Branch("vtxid" ,_vtxid ,"vtxid/I");
  }
  
  
}


#endif
