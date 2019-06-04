#ifndef __INTERPMT_CXX__
#define __INTERPMT_CXX__

#include "InterPMT.h"

// ll
#include "DataFormat/vertex.h"
#include "DataFormat/opdetwaveform.h"

namespace llcv {

  void InterPMT::configure(const larcv::PSet& cfg) {
    LLCV_DEBUG() << "start" << std::endl;

    _adc_img_prod = cfg.get<std::string>("ADCImageProducer");
    _trk_img_prod = cfg.get<std::string>("TrackImageProducer");
    _shr_img_prod = cfg.get<std::string>("ShowerImageProducer");

    _track_vertex_prod   = cfg.get<std::string>("TrackVertexProducer");
    _opdigit_prod  = cfg.get<std::string>("OpDigitProducer","saturation");

    LLCV_DEBUG() << "adc_img_prod........." << _adc_img_prod << std::endl;
    LLCV_DEBUG() << "trk_img_prod........." << _trk_img_prod << std::endl;
    LLCV_DEBUG() << "shr_img_prod........." << _shr_img_prod << std::endl;
    LLCV_DEBUG() << "opdigit_prod............." << _opdigit_prod     << std::endl;
    
    _epsilon = cfg.get<float>("EPS",1e-5);

    _driver.Configure(cfg.get<larcv::PSet>(_driver.name()));

    _driver._tree_mgr.Configure(cfg.get<larcv::PSet>(_driver._tree_mgr.name()));

    for(auto selptr : _driver._sel_base_v)
      selptr->Configure(cfg.get<larcv::PSet>(selptr->name()));

    LLCV_DEBUG() << "end" << std::endl;
  }

  void InterPMT::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    _driver.Initialize();
    LLCV_DEBUG() << "end" << std::endl;
  }

  bool InterPMT::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "@sto (r,s,e,e)=(" 
		 << sto.run_id()    << "," 
		 << sto.subrun_id() << "," 
		 << sto.event_id()  << "," 
		 << sto.get_index() << ")" <<std::endl;

    //
    // larcv data products
    //
    larcv::EventImage2D* ev_adc_img = nullptr;
    if (!_adc_img_prod.empty()) 
      ev_adc_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_adc_img_prod);
    else 
      throw llcv_err("Must specify ADC image");


    LLCV_DEBUG() << "@mgr (r,s,e,e)=(" 
		 << ev_adc_img->run()    << "," 
		 << ev_adc_img->subrun() << "," 
		 << ev_adc_img->event()  << "," 
		 << mgr.current_entry()  << ")" << std::endl;

    _driver._run    = (int)ev_adc_img->run();
    _driver._subrun = (int)ev_adc_img->subrun();
    _driver._event  = (int)ev_adc_img->event();

    _driver.AttachImage(ev_adc_img->Image2DArray(),kImageADC);

    larcv::EventImage2D* ev_trk_img = nullptr;
    if(!_trk_img_prod.empty()){
      ev_trk_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_trk_img_prod);
      _driver.AttachImage(ev_trk_img->Image2DArray(),kImageTrack);
    }
    
    larcv::EventImage2D* ev_shr_img = nullptr;
    if(!_shr_img_prod.empty()) {
      ev_shr_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_shr_img_prod);
      _driver.AttachImage(ev_shr_img->Image2DArray(),kImageShower);
    }
    
    //
    // larlite data products
    //
    larlite::event_vertex* ev_track_vertex = nullptr;
    if(_track_vertex_prod.empty()) 
      throw llcv_err("Empty track vertex producer specified");
      
    ev_track_vertex = (larlite::event_vertex*)sto.get_data(larlite::data::kVertex,_track_vertex_prod);



    larlite::event_opdetwaveform* ev_opdigit = nullptr;
    if(!_opdigit_prod.empty()) 
      ev_opdigit = (larlite::event_opdetwaveform*)sto.get_data(larlite::data::kOpDetWaveform,_opdigit_prod);
    
    //
    // configure the driver
    //
    size_t num_vertex = ev_track_vertex->size();
    for(size_t vtxid=0; vtxid < num_vertex; ++vtxid) {
      const auto track_vertex_ptr = &((*ev_track_vertex)[vtxid]);
      auto vid  = _driver.AttachVertex(track_vertex_ptr);
      for(const auto& opdigit : (*ev_opdigit))
	_driver.AttachOpDigit(vid,&opdigit);
    }
    
    _driver.Process();
    _driver.Reset();
    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  
  void InterPMT::finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _driver.Finalize();
    LLCV_DEBUG() << "end" << std::endl;
  }

}

#endif
	
