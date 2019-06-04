#ifndef LLCVANABASE_CXX
#define LLCVANABASE_CXX

#include "AnaBase.h"

namespace llcv {

  AnaBase::AnaBase(const std::string name)
    : llcv_base  ( name ),
      _name      ( name )
    , _fout      ( nullptr )
  {}

  void AnaBase::_configure_(const larcv::PSet& cfg)
  {
    set_verbosity((msg::Level_t)(cfg.get<unsigned short>("Verbosity",logger().level())));
    configure(cfg);
  }
  
  bool AnaBase::_process_(larcv::IOManager& mgr, larlite::storage_manager& sto)
  {
    bool status=false;
    status = this->process(mgr,sto);
    return status;
  }


}

#endif
