#ifndef __INTERTTREEMANAGER_H__
#define __INTERTTREEMANAGER_H__

#include "LLCVBase/llcv_base.h"
#include "Base/PSet.h"
#include "InterTTreeSpec.h"

#include <map>

namespace llcv {


  class RSEVID {
  public:
  RSEVID(size_t run_val=0, size_t subrun_val=0, size_t event_val=0, size_t vertex_val=0)
    : run(run_val)
      , subrun(subrun_val)
      , event(event_val)
      , vertex(vertex_val)
    {}
    ~RSEVID(){}

    inline bool operator < (const RSEVID& rhs) const
    { if(run < rhs.run) return true;
      if(run > rhs.run) return false;
      if(subrun < rhs.subrun) return true;
      if(subrun > rhs.subrun) return false;
      if(event < rhs.event) return true;
      if(event > rhs.event) return false;
      if(vertex < rhs.vertex) return true;
      if(vertex > rhs.vertex) return false;
      return false;
    }

    inline bool operator == (const RSEVID& rhs) const
    { return (run == rhs.run && subrun == rhs.subrun && event == rhs.event && vertex == rhs.vertex); }

    inline bool operator != (const RSEVID& rhs) const
    { return !( (*this) == rhs ); }

    inline bool operator > (const RSEVID& rhs) const
    { return ( (*this) != rhs && !((*this) < rhs) ); }

    size_t run, subrun, event, vertex;
  };

  class InterModule;
  class InterMichel;
  class InterPMT;
  class InterDriver;

  class InterTTreeManager : public llcv_base {

    friend class InterModule;
    friend class InterMichel;
    friend class InterPMT;

    friend class InterDriver;

  public:
    
  InterTTreeManager(std::string name="InterTTreeManager")
    :  llcv_base(name)
      , _name(name)
      , _tchain(nullptr)
      , _nentries(-1)
      , _centry(kINVALID_SIZE)
    {}
    
  InterTTreeManager(const std::string& fname,
		    const std::string& tname,
		    std::string name="InterTTreeManager")
    :  llcv_base(name)
      , _name(name)
      , _tchain(nullptr)
      , _nentries(-1)
      , _centry(kINVALID_SIZE)
    {
      Initialize(fname,tname);
    }
    
    ~InterTTreeManager() {}

  public:
    
    template <class T> T Scalar(const std::string& str) const;
    
    std::vector<float> Vector(const std::string& str) const;

    std::vector<std::vector<float> > VVector(const std::string& str) const;
    
    
  private:

    std::string _name;

    TChain* _tchain;

    int _nentries;
    size_t _centry;

    InterTTreeSpec _spec;

    std::map<RSEVID,size_t> _rsev_m;

    void Configure(const larcv::PSet& cfg);
    void Initialize(const std::string& fname, const std::string& tname);

    bool GoTo(size_t run, size_t subrun, size_t event, size_t vtxid);
    bool _goto(const RSEVID& rsev);

    size_t Run()    const { return (size_t)(*(_spec._run_ptr)); }
    size_t SubRun() const { return (size_t)(*(_spec._subrun_ptr)); }
    size_t Event()  const { return (size_t)(*(_spec._event_ptr)); }
    size_t Vertex() const { return (size_t)(*(_spec._vtxid_ptr)); }


  };
}


#endif
