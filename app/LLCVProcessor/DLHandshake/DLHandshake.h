#ifndef __DLHANDSHAKE_H__
#define __DLHANDSHAKE_H__

#include "LLCVBase/AnaBase.h"
#include "HandShaker.h"

namespace llcv {

  class DLHandshake : public AnaBase {

  public:

  DLHandshake(const std::string name="DLHandshake") : AnaBase(name) {}
    ~DLHandshake() {}

    void configure(const larcv::PSet&);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();

  private:

    HandShaker _HandShaker;

    std::string _in_hit_prod;
    std::string _in_pgraph_prod;
    std::string _in_ctor_prod;
    std::string _out_prod;

    bool _use_ctor;

  };

}

#endif
