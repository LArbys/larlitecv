#ifndef __PROCESSOR_H__
#define __PROCESSOR_H__

// package
#include "llcv_base.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "AnaBase.h"

// larcv
#include "Processor/ProcessBase.h"

// larlite
#include "Analysis/ana_base.h"

namespace llcv {

  class Processor : public llcv_base {
  public:
    
  Processor(const std::string name="Processor") : 
    llcv_base(name), _name(name) {}

    ~Processor() {}
    
    void configure(const std::string& fname);
    void initialize();
    bool process_event(int run, int subrun, int event);
    void finalize();
    
    void add_ll_input_file(const std::string& fname);
    void add_lcv_input_file(const std::string& fname);

    void add_ll_input_files(const std::vector<std::string>& fname_v);
    void add_lcv_input_files(const std::vector<std::string>& fname_v);
    
    void set_output_ll_name(const std::string& fname);
    void set_output_lcv_name(const std::string& fname);

    bool batch_process_lcv(int start=0, int nentries=-1);
    bool batch_process_ll(int start=0, int nentries=-1);
 
    bool batch_process_lcv_reverse(int start=0, int nentries=-1);
    bool batch_process_ll_reverse(int start=0, int nentries=-1);

    larlitecv::DataCoordinator& dataco() { return _dataco; }
    const std::string& name() const { return _name; }
    
    size_t get_n_ll_entries()  { return (size_t)_dataco.get_nentries("larlite"); }
    size_t get_n_lcv_entries() { return (size_t)_dataco.get_nentries("larcv"); }

    void add_llcv_ana(AnaBase* ptr);
    void add_ll_ana(larlite::ana_base* ptr);
    void add_lc_proc(larcv::ProcessBase* ptr);

  private:
    //whoami
    std::string _name;

    //coordinator
    larlitecv::DataCoordinator _dataco;

    //llcv
    std::vector<AnaBase*> _llcv_ana_v;
    std::vector<bool> _llcv_status_v;
    bool _llcv_unit_status;

    //ll
    std::vector<larlite::ana_base*> _ll_ana_v;
    std::vector<bool> _ll_status_v;
    bool _ll_unit_status;
    
    //lcv
    std::vector<larcv::ProcessBase*> _lcv_proc_v;
    std::vector<bool> _lcv_status_v;
    bool _lcv_unit_status;

    // internal batch process
    bool _batch_process(const std::string& ftype,
			int start,int nentries);

    bool _batch_process_reverse(const std::string& ftype,
				int start,int nentries);

  };
}
#endif
