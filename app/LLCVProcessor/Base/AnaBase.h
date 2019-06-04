#ifndef LLCVANABASE_H
#define LLCVANABASE_H

#include "llcv_base.h"
#include "llcv_err.h"

#include "DataFormat/IOManager.h"
#include "DataFormat/storage_manager.h"
#include "Base/PSet.h"

namespace llcv {

  class Processor;
  
  class AnaBase : public llcv_base {

    friend class Processor;

  public:
    
    /// Default constructor
    AnaBase(const std::string name="AnaBase");
    
    /// Default destructor
    virtual ~AnaBase(){}

    virtual void configure(const larcv::PSet&) = 0;
    virtual void initialize() = 0;
    virtual bool process(larcv::IOManager& mgr, larlite::storage_manager& sto) = 0;
    virtual void finalize() = 0;

    bool has_ana_file() const
    { return _fout != nullptr; }
    
    TFile& ana_file()
      { if(!_fout) throw llcv_err("ana file does not exist"); return *_fout; }

    const std::string& name() const { return _name; }

  private:
    void _configure_(const larcv::PSet&);
    bool _process_(larcv::IOManager& mgr, larlite::storage_manager& sto);
    bool _event_creator;    ///< special flag to mark this algorithm an event creator
    std::string _name;
    TFile* _fout;           ///< output analysis file

  };
}

#endif
/** @} */ // end of doxygen group 

