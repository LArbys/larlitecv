#ifndef __DLHANDSHAKE_CXX__
#define __DLHANDSHAKE_CXX__

#include "DLHandshake.h"

namespace llcv {

  void DLHandshake::configure(const larcv::PSet& cfg) {
    LLCV_DEBUG() << "start" << std::endl;

    _in_hit_prod     = cfg.get<std::string>("InputHitProducer", "gaushit");
    _in_pgraph_prod  = cfg.get<std::string>("InputPGraphProducer", "test");
    _in_ctor_prod    = cfg.get<std::string>("InputCtorProducer", "test_ctor");
    _out_prod        = cfg.get<std::string>("OutputProducer","dl");
    _use_ctor        = cfg.get<bool>("UseContour",true);

    LLCV_DEBUG() << "end" << std::endl;
  }

  void DLHandshake::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    _HandShaker.reset();
    LLCV_DEBUG() << "end" << std::endl;
  }

  bool DLHandshake::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {
    LLCV_DEBUG() << "start" << std::endl;

    _HandShaker.reset();

    auto ev_pfpart  = (larlite::event_pfpart*)  sto.get_data(larlite::data::kPFParticle, _out_prod);
    auto ev_vertex  = (larlite::event_vertex*)  sto.get_data(larlite::data::kVertex,     _out_prod);
    auto ev_shower  = (larlite::event_shower*)  sto.get_data(larlite::data::kShower,     _out_prod);
    auto ev_track   = (larlite::event_track*)   sto.get_data(larlite::data::kTrack,      _out_prod);
    auto ev_cluster = (larlite::event_cluster*) sto.get_data(larlite::data::kCluster,    _out_prod);
    auto ev_hit     = (larlite::event_hit*)     sto.get_data(larlite::data::kHit,        _out_prod);
    auto ev_ass     = (larlite::event_ass*)     sto.get_data(larlite::data::kAssociation,_out_prod);
    
    if(!ev_pfpart->empty())  throw llcv_err("LL dataproduct not empty");
    if(!ev_vertex->empty())  throw llcv_err("LL dataproduct not empty");
    if(!ev_shower->empty())  throw llcv_err("LL dataproduct not empty");
    if(!ev_track->empty())   throw llcv_err("LL dataproduct not empty");
    if(!ev_cluster->empty()) throw llcv_err("LL dataproduct not empty");
    if(!ev_hit->empty())     throw llcv_err("LL dataproduct not empty");

    auto ev_hit_in  = (larlite::event_hit*)     sto.get_data(larlite::data::kHit, _in_hit_prod);
    LLCV_DEBUG() << "GOT: " << ev_hit_in->size() << " gaushit" << std::endl;
        
    auto ev_pgraph  = (larcv::EventPGraph*) mgr.get_data(larcv::kProductPGraph, _in_pgraph_prod);
    auto ev_pixel2d = (larcv::EventPixel2D*)mgr.get_data(larcv::kProductPixel2D,_in_ctor_prod);

    LLCV_DEBUG() << "GOT: " << ev_pgraph->PGraphArray().size() << " vertices" << std::endl;
    LLCV_DEBUG() << "GOT: " << ev_pixel2d->Pixel2DArray().size() << " pixel array" << std::endl;
    LLCV_DEBUG() << "GOT: " << ev_pixel2d->Pixel2DClusterArray().size() << " pixel cluster array" << std::endl;


    _HandShaker.pixel_distance_threshold(1.);
    _HandShaker.set_larlite_pointers(ev_pfpart, ev_vertex,
				     ev_shower, ev_track,
				     ev_cluster, ev_hit,
				     ev_ass);
      
    if (_use_ctor) { 
      LLCV_DEBUG() << "Using contours" << std::endl;
      _HandShaker.construct_contour(*ev_pgraph, *ev_pixel2d, ev_hit_in);
    }
    else {
      LLCV_DEBUG() << "Using image" << std::endl;
      _HandShaker.construct_image(*ev_pgraph, *ev_pixel2d, ev_hit_in);
    }

    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  
  void DLHandshake::finalize() {

    LLCV_DEBUG() << "start" << std::endl;
      
    LLCV_DEBUG() << "end" << std::endl;
  }
}

#endif
	
