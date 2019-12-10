#ifndef __PROCESSOR_CXX__
#define __PROCESSOR_CXX__

#include "Processor.h"
#include "Base/LArCVBaseUtilFunc.h"
#include <cassert>

namespace llcv {

  void Processor::configure(const std::string& fname) 
  {
    const auto main_cfg = larcv::CreatePSetFromFile(fname);
    // get my pset
    const auto proc_cfg = main_cfg.get<larcv::PSet>(name());
    set_verbosity((msg::Level_t)(proc_cfg.get<unsigned short>("Verbosity",logger().level())));

    LLCV_INFO() << "start" << std::endl;

    // get dataco pset
    const auto dataco_config  = proc_cfg.get<larcv::PSet>("DataCoordinator");
    auto dataco_ll_cfg  = dataco_config.get<larcv::PSet>("StorageManager");
    auto dataco_lcv_cfg = dataco_config.get<larcv::PSet>("IOManager");

    _dataco.configure(dataco_lcv_cfg,dataco_ll_cfg);
    
    // get llcv & lcv algorithm configure blocks
    for(auto ptr : _llcv_ana_v) {
      LLCV_DEBUG() << "configure="<<ptr->name()<<std::endl;
      ptr->_configure_(proc_cfg.get<larcv::PSet>(ptr->name()));
    }

    for(auto ptr : _lcv_proc_v) {
      LLCV_DEBUG() << "configure="<<ptr->name()<<std::endl;
      ptr->configure(proc_cfg.get<larcv::PSet>(ptr->name()));
    }
    
    LLCV_INFO() << "end" << std::endl;
  }

  void Processor::initialize() {
    LLCV_INFO() << "start" << std::endl;
    _dataco.initialize();

    for(auto ptr : _llcv_ana_v)
      ptr->initialize();

    for(auto ptr : _ll_ana_v)
      ptr->initialize();

    for(auto ptr : _lcv_proc_v)
      ptr->initialize();

    LLCV_INFO() << "end" << std::endl;
  }
  
  bool Processor::batch_process_lcv(int start,int nentries)
  { 
    LLCV_INFO() << "start" << std::endl;
    return _batch_process("larcv",start,nentries);
    LLCV_INFO() << "end" << std::endl;
  }

  bool Processor::batch_process_ll(int start,int nentries)
  { 
    LLCV_INFO() << "start" << std::endl;
    return _batch_process("larlite",start,nentries);
    LLCV_INFO() << "end" << std::endl;
  }

  bool Processor::batch_process_lcv_reverse(int start,int nentries)
  { 
    LLCV_INFO() << "start" << std::endl;
    return _batch_process_reverse("larcv",start,nentries);
    LLCV_INFO() << "end" << std::endl;
  }

  bool Processor::batch_process_ll_reverse(int start,int nentries)
  { 
    LLCV_INFO() << "start" << std::endl;
    return _batch_process_reverse("larlite",start,nentries);
    LLCV_INFO() << "end" << std::endl;
  }
  
  bool Processor::_batch_process(const std::string& ftype,int start, int nentries) {
    int type_nentries = (int)_dataco.get_nentries(ftype);

    if (start > type_nentries) 
      throw llcv_err("requested index out of range");
    
    if (nentries<0) 
      nentries = type_nentries;
    
    nentries += start;
    
    if (nentries>(start+type_nentries))
      nentries = type_nentries;

    LLCV_INFO() << "Processing " << type_nentries << " entries @ ftype=" << ftype  << std::endl;

    for(int entry=start; entry<nentries; ++entry) {
      
      LLCV_DEBUG() << "@entry=" << entry << std::endl;
      _dataco.goto_entry(entry,ftype);
      
      // ll & llcv
      LLCV_DEBUG() << "processing ll+lcv" << std::endl;
      for(size_t aid=0; aid < _llcv_ana_v.size(); ++aid) {
	LLCV_DEBUG() << "@id=" << aid << " name=" << _llcv_ana_v[aid]->name() << " (" << _llcv_ana_v[aid] << ")"  << std::endl;
	try {
	_llcv_status_v[aid] = _llcv_ana_v[aid]->_process_(_dataco.get_larcv_io(),_dataco.get_larlite_io());
	}
	catch ( std::exception& e ) {
	  LLCV_CRITICAL() << "@id=" << aid << " name=" << _llcv_ana_v[aid]->name() << " (" << _llcv_ana_v[aid] << ")"  << std::endl
			  << "ERROR: " << e.what() << std::endl;
	}
	_llcv_unit_status = _llcv_unit_status && _llcv_status_v[aid];
      }

      // ll
      LLCV_DEBUG() << "processing ll" << std::endl;
      for(size_t aid=0; aid < _ll_ana_v.size(); ++aid) {
	LLCV_DEBUG() << "@id=" << aid << " name=" << _ll_ana_v[aid]->class_name() << " (" << _ll_ana_v[aid] << ")"  << std::endl;
	_ll_status_v[aid] = _ll_ana_v[aid]->analyze(&_dataco.get_larlite_io());
	_ll_unit_status = _ll_unit_status && _ll_status_v[aid];
      }
    
      //lcv
      LLCV_DEBUG() << "processing lcv" << std::endl;
      for(size_t aid=0; aid < _lcv_proc_v.size(); ++aid) {
	LLCV_DEBUG() << "@id=" << aid << " name=" << _lcv_proc_v[aid]->name() << " (" << _lcv_proc_v[aid] << ")"  << std::endl;
	_lcv_status_v[aid] = _lcv_proc_v[aid]->process(_dataco.get_larcv_io());
	_lcv_unit_status = _lcv_unit_status && _lcv_status_v[aid];
      }

      // (?)
      _dataco.save_entry();
    }
    
    return true;
  }

  bool Processor::_batch_process_reverse(const std::string& ftype,int start, int nentries) {
    int type_nentries = (int)_dataco.get_nentries(ftype);

    if (start > type_nentries) 
      throw llcv_err("requested index out of range");
    
    if (nentries<0) 
      nentries = type_nentries;
    
    nentries += start;
    
    if (nentries>(start+type_nentries))
      nentries = type_nentries;

    LLCV_INFO() << "Processing " << type_nentries << " entries @ ftype=" << ftype  << std::endl;

    for(int entry=start; entry<nentries; ++entry) {
      
      LLCV_DEBUG() << "@entry=" << entry << std::endl;
      _dataco.goto_entry(entry,ftype);
      
      _ll_unit_status   = true;
      _lcv_unit_status  = true;
      _llcv_unit_status = true;

      // ll
      LLCV_DEBUG() << "processing ll" << std::endl;
      for(size_t aid=0; aid < _ll_ana_v.size(); ++aid) {
	LLCV_DEBUG() << "@id=" << aid << " name=" << _ll_ana_v[aid]->class_name() << " (" << _ll_ana_v[aid] << ")"  << std::endl;
	try {
	  _ll_status_v[aid] = _ll_ana_v[aid]->analyze(&_dataco.get_larlite_io());
	}
	catch ( std::exception& e ) {
	  LLCV_CRITICAL() << "@id=" << aid << " name=" << _llcv_ana_v[aid]->name() << " (" << _llcv_ana_v[aid] << ")"  << std::endl
			  << "ERROR: " << e.what() << std::endl;
	}	  
	_ll_unit_status = _ll_unit_status && _ll_status_v[aid];
	if(!_ll_unit_status) break;
      }
    
      //lcv
      LLCV_DEBUG() << "processing lcv" << std::endl;
      for(size_t aid=0; aid < _lcv_proc_v.size(); ++aid) {
	LLCV_DEBUG() << "@id=" << aid << " name=" << _lcv_proc_v[aid]->name() << " (" << _lcv_proc_v[aid] << ")"  << std::endl;
	_lcv_status_v[aid] = _lcv_proc_v[aid]->process(_dataco.get_larcv_io());
	_lcv_unit_status = _lcv_unit_status && _lcv_status_v[aid];
	if(!_lcv_unit_status) break;
      }
      
      // ll & llcv
      LLCV_DEBUG() << "processing ll+lcv" << std::endl;
      for(size_t aid=0; aid < _llcv_ana_v.size() and _lcv_unit_status and _ll_unit_status; ++aid) {
	LLCV_DEBUG() << "@id=" << aid << " name=" << _llcv_ana_v[aid]->name() << " (" << _llcv_ana_v[aid] << ")"  << std::endl;
	_llcv_status_v[aid] = _llcv_ana_v[aid]->_process_(_dataco.get_larcv_io(),_dataco.get_larlite_io());
	_llcv_unit_status = _llcv_unit_status && _llcv_status_v[aid];
	if(!_llcv_unit_status) break;
      }

      // (?)
      _dataco.save_entry();
    }
    
    return true;
  }

  void Processor::finalize() {
    LLCV_INFO() << "start" << std::endl;
    
    for(auto ptr : _llcv_ana_v)
      ptr->finalize();

    for(auto ptr : _ll_ana_v)
      ptr->finalize();

    for(auto ptr : _lcv_proc_v)
      ptr->finalize();

    _dataco.finalize();
    LLCV_INFO() << "end" << std::endl;
  }

  bool Processor::process_event(int run, int subrun, int event)
  { return true; }
  
  void Processor::set_output_ll_name(const std::string& fname) 
  {
    LLCV_DEBUG() << "set ll output name=" << fname << std::endl;
    _dataco.get_larlite_io().set_out_filename(fname);
  }

  void Processor::set_output_lcv_name(const std::string& fname) 
  {
    LLCV_DEBUG() << "set lcv output name=" << fname << std::endl;
    _dataco.get_larcv_io().set_out_file(fname);
  }

  void Processor::add_ll_input_files(const std::vector<std::string>& fname_v)  
  {
    for(const auto& fname : fname_v) 
      add_ll_input_file(fname);
  }
  
  void Processor::add_lcv_input_files(const std::vector<std::string>& fname_v) 
  {
    for(const auto& fname : fname_v) 
      add_lcv_input_file(fname);
  }  
  
  void Processor::add_ll_input_file(const std::string& fname) 
  {
    LLCV_INFO() << "insert ll fname=" << fname << std::endl;
    _dataco.add_inputfile(fname, "larlite");
  }
  
  void Processor::add_lcv_input_file(const std::string& fname) 
  {
    LLCV_INFO() << "insert lcv fname=" << fname << std::endl;
    _dataco.add_inputfile(fname, "larcv");
  }
  
  void Processor::add_llcv_ana(AnaBase* ptr) {
    _llcv_ana_v.push_back(ptr);
    _llcv_status_v.push_back(true);
    assert(_llcv_ana_v.size() == _llcv_status_v.size());
  }

  void Processor::add_ll_ana(larlite::ana_base* ptr) {
    _ll_ana_v.push_back(ptr);
    _ll_status_v.push_back(true);
    assert(_ll_ana_v.size() == _ll_status_v.size());
  }

  void Processor::add_lc_proc(larcv::ProcessBase* ptr) {
    _lcv_proc_v.push_back(ptr);
    _lcv_status_v.push_back(true);
    assert(_lcv_proc_v.size() == _lcv_status_v.size());
  }
  
}

#endif
