#ifndef SHOWERTRUTHMATCH_CXX
#define SHOWERTRUTHMATCH_CXX

// llcv
#include "ShowerTruthMatch.h"
#include "MatchUtils.h"

// larcv
#include "DataFormat/EventImage2D.h"

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/shower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include <array>
#include <cassert>

namespace llcv {

  void ShowerTruthMatch::configure(const larcv::PSet& cfg) {
    LLCV_DEBUG() << "start" << std::endl;

    _shr_reco_prod = cfg.get<std::string>("ShowerRecoProducer");
    _adc_img_prod  = cfg.get<std::string>("ADCImageProducer");
    _seg_img_prod  = cfg.get<std::string>("TrueImageProducer");


    LLCV_INFO() << "ShowerRecoProducer: " << _shr_reco_prod << std::endl;
    LLCV_INFO() << "ADCImageProducer:  " << _adc_img_prod << std::endl;
    LLCV_INFO() << "TrueImageProducer: " << _seg_img_prod << std::endl;

    LLCV_DEBUG() << "end" << std::endl;
  }

  void ShowerTruthMatch::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    _tree = new TTree("ShowerTruthMatch","");
    _tree->Branch("run"     , &_run     , "run/I");
    _tree->Branch("subrun"  , &_subrun  , "subrun/I");
    _tree->Branch("event"   , &_event   , "event/I");
    _tree->Branch("entry"   , &_entry   , "entry/I");
    _tree->Branch("vtxid"   , &_vtxid   , "vtxid/I");
    _tree->Branch("nshowers", &_nshowers, "nshowers/I");

    _tree->Branch("unknownfrac_v" , &_unknownfrac_v);
    _tree->Branch("electronfrac_v", &_electronfrac_v);
    _tree->Branch("gammafrac_v"   , &_gammafrac_v);
    _tree->Branch("pizerofrac_v"  , &_pizerofrac_v);
    _tree->Branch("muonfrac_v"    , &_muonfrac_v);
    _tree->Branch("kminusfrac_v"  , &_kminusfrac_v);
    _tree->Branch("piminusfrac_v" , &_piminusfrac_v);
    _tree->Branch("protonfrac_v"  , &_protonfrac_v);

    _tree->Branch("npx_v"      , &_npx_v);
    _tree->Branch("shr_type_v" , &_shr_type_v);
    _tree->Branch("shr_type_vv", &_shr_type_vv);

    LLCV_DEBUG() << "end" << std::endl;
  }

  bool ShowerTruthMatch::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {

    if (_seg_img_prod.empty()) return true;

    assert(!_shr_reco_prod.empty());
    assert(!_adc_img_prod.empty());
    
    const auto ev_vertex  = (larlite::event_vertex*)sto.get_data(larlite::data::kVertex,_shr_reco_prod);
    const auto ev_adc_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_adc_img_prod);
    const auto ev_seg_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_seg_img_prod);

    LLCV_DEBUG() << "LC RSE=(" << ev_seg_img->run() << "," << ev_seg_img->subrun() << "," << ev_seg_img->event() << ") " 
		 << "& LL RSE=(" << sto.run_id()  << "," << sto.subrun_id() << "," << sto.event_id() << ")" << std::endl;
    
    LLCV_DEBUG() << "GOT: " << ev_vertex->size() << " vertices" << std::endl;
    
    _run    = (int) ev_adc_img->run();
    _subrun = (int) ev_adc_img->subrun();
    _event  = (int) ev_adc_img->event();
    _entry  = (int) mgr.current_entry();
    
    // check data coordinator...
    assert(_run    == (int)sto.run_id());
    assert(_subrun == (int)sto.subrun_id());
    assert(_event  == (int)sto.event_id());

    if (ev_vertex->empty()) {
      LLCV_INFO() << "No vertex, next" << std::endl;
      return true;
    }

    // get the associated pf particles
    larlite::event_pfpart *ev_pfpart = nullptr;
    auto const& ass_pfpart_vv = sto.find_one_ass(ev_vertex->id(), ev_pfpart, ev_vertex->name());
    if (!ev_pfpart or ev_pfpart->empty()) {
      LLCV_INFO()  << "No pfpart, next" << std::endl;
      return true;
    }

    LLCV_INFO() << "... got " << ev_pfpart->size() << " pfparts" << std::endl;

    // get the associated pf clusters
    larlite::event_shower *ev_shower = nullptr;
    auto const& ass_shower_vv = sto.find_one_ass(ev_vertex->id(), ev_shower, "showerreco");
    //if (!ev_shower or ev_shower->empty()) {
    //  LLCV_INFO()  << "No shower, next" << std::endl;
    //  return true;
    //}
    
    larlite::event_cluster * ev_cluster = nullptr;
    larlite::event_hit * ev_hit = nullptr;
    auto ass_cluster_vv = ass_shower_vv; // i forget the type
    ass_cluster_vv.clear();
    auto ass_hit_vv = ass_cluster_vv;
 
    if (ev_shower) { 
      LLCV_INFO() << "... got " << ev_shower->size() << " showers" << std::endl;
  
      // get the associated cluster(s)
      ass_cluster_vv = sto.find_one_ass(ev_shower->id(), ev_cluster, ev_shower->name());
      if (!ev_cluster or ev_cluster->empty()) {
        LLCV_CRITICAL() << "NO associated cluster to shower" << std::endl;
        throw llcv_err("die");
      }

      LLCV_INFO() << "... got " << ev_cluster->size() << " clusters" << std::endl;

      // get the associated hit(s)
      ass_hit_vv = sto.find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());
      if (!ev_hit or ev_hit->empty()) {
        LLCV_CRITICAL() << "NO associated hits to cluster" << std::endl;
        throw llcv_err("die");
      }

      LLCV_INFO() << "... got " << ev_hit->size() << " hits" << std::endl;
    } 

    std::array<larcv::ImageMeta,3> meta_v;
    for(size_t plane=0; plane<3; ++plane) 
      meta_v[plane] = ev_seg_img->Image2DArray()[plane].meta();
    
    double xpixel = kINVALID_DOUBLE;
    double ypixel = kINVALID_DOUBLE;

    for( size_t vtx_id = 0; vtx_id < ass_pfpart_vv.size(); ++vtx_id) {
      ClearVertex();
      const auto& vertex = ev_vertex->at(vtx_id);
      const auto& ass_pfpart_v = ass_pfpart_vv.at(vtx_id);

      LLCV_DEBUG() << "@vtx_id=" << vtx_id << " : (" << vertex.X() << "," << vertex.Y() << "," << vertex.Z() << ")" << std::endl;
      LLCV_DEBUG() << "..." << ass_pfpart_v.size() << " associated pfparts" << std::endl;

      _vtxid = vtx_id;

      _nshowers = 0;
      

      if(!ev_shower) {
        _tree->Fill(); 
        continue;
      }
          
      std::vector<const larlite::shower* > shower_v;
      std::vector<std::array<std::vector<const larlite::hit*>, 3> > hit_vvv;

      hit_vvv.reserve(ev_shower->size());

      for( size_t pfp_id=0; pfp_id < ass_pfpart_v.size(); ++ pfp_id) {

	const auto pfpart_id = ass_pfpart_v.at(pfp_id);
	
	if (pfpart_id >= ass_shower_vv.size()) {
	  std::cout << std::endl;
	  std::cout << "---------------------------" << std::endl;
	  std::cout << "Edge case detected!" << std::endl;
	  std::cout << "Vertex found, pfparticle exists, but no shower found?" << std::endl;
	  std::cout << "Skip for now..." << std::endl;
	  std::cout << "---------------------------" << std::endl << std::endl;
	  continue;
	}

	const auto& ass_shower_v = ass_shower_vv.at(pfpart_id);

	for(size_t shr_id=0; shr_id < ass_shower_v.size(); ++ shr_id) {

	  const auto shower_id = ass_shower_v.at(shr_id);
	  const auto& ass_cluster_v = ass_cluster_vv.at(shower_id);

	  shower_v.push_back(&ev_shower->at(shower_id));
	  hit_vvv.resize(hit_vvv.size()+1);
	  auto& hit_vv = hit_vvv.back();

	  for(size_t clu_id=0; clu_id < ass_cluster_v.size(); ++ clu_id) {
	    const auto cluster_id = ass_cluster_v.at(clu_id);
	    const auto& ass_hit_v = ass_hit_vv.at(cluster_id);
	    
	    for(size_t h_id=0; h_id < ass_hit_v.size(); ++ h_id) {
	      const auto hit_id = ass_hit_v.at(h_id);
	      const auto& hit = ev_hit->at(hit_id);
	      hit_vv.at((size_t)hit.View()).push_back(&hit);
	    }
	  }
	}
      }

      LLCV_DEBUG() << shower_v.size() << ".. showers" << std::endl;
      for(size_t shr_id=0; shr_id<shower_v.size(); ++shr_id) {
	LLCV_DEBUG() << "@shrid=" << shr_id << std::endl;
	for(auto plane_v : hit_vvv.at(shr_id)) {
	  LLCV_DEBUG() << "..." << plane_v.size() << std::endl;
	}
      }

      _nshowers = (int)shower_v.size();
      ResizeFrac(_nshowers);
      LLCV_DEBUG() << "loop" << std::endl;

      for(size_t shr_id=0; shr_id < shower_v.size(); ++shr_id) {

	_npx_v[shr_id] = 0;
	auto& shr_type = _shr_type_v[shr_id];
	auto& shr_type_v = _shr_type_vv[shr_id];

	for(size_t plane=0; plane<3; ++plane) {
	  const auto& adc_img = ev_adc_img->Image2DArray().at(plane);
	  const auto& seg_img = ev_seg_img->Image2DArray().at(plane);
	  const auto& hit_v  = hit_vvv.at(shr_id).at(plane);
	  const auto& meta = meta_v.at(plane);

	  for(const auto hit : hit_v) {
	    assert(hit);
	    auto time = hit->PeakTime() + 2400;
	    auto wire = hit->WireID().Wire;
	    
	    auto xpixel = (wire - meta.min_x()) / meta.pixel_width();
	    auto ypixel = (meta.max_y() - time) / meta.pixel_height();

	    int xx = (int)(xpixel+0.5);
	    int yy = (int)(ypixel+0.5);

	    // float pixel_value = adc_img.pixel(yy,xx);
	    // assert(pixel_value);
	    float pixel_type  = seg_img.pixel(yy,xx);
	    shr_type_v.at((size_t)pixel_type) += 1;
	    _npx_v[shr_id] += 1;
	    LLCV_DEBUG() << "p:" << plane << "(" << yy << "," << xx << ")=" << (size_t)pixel_type << std::endl;
	  }
	} // end this plane 



	auto res_iter = std::max_element(std::begin(shr_type_v), std::end(shr_type_v));
	shr_type  = (int)std::distance(std::begin(shr_type_v), res_iter);
      
	auto& unknownfrac  = _unknownfrac_v[shr_id];
	auto& electronfrac = _electronfrac_v[shr_id];
	auto& gammafrac    = _gammafrac_v[shr_id];
	auto& pizerofrac   = _pizerofrac_v[shr_id];
	auto& muonfrac     = _muonfrac_v[shr_id];
	auto& kminusfrac   = _kminusfrac_v[shr_id];
	auto& piminusfrac  = _piminusfrac_v[shr_id];
	auto& protonfrac   = _protonfrac_v[shr_id];
      
	unknownfrac  = 0.0;
	electronfrac = 0.0;
	gammafrac    = 0.0;
	pizerofrac   = 0.0;
	muonfrac     = 0.0;
	kminusfrac   = 0.0;
	piminusfrac  = 0.0;
	protonfrac   = 0.0;

	float npixels = _npx_v[shr_id];
	if (npixels) {
	  unknownfrac  = (float) shr_type_v.at(0) / npixels;
	  electronfrac = (float) shr_type_v.at(3) / npixels;
	  gammafrac    = (float) shr_type_v.at(4) / npixels;
	  pizerofrac   = (float) shr_type_v.at(5) / npixels;
	  muonfrac     = (float) shr_type_v.at(6) / npixels;
	  kminusfrac   = (float) shr_type_v.at(7) / npixels;
	  piminusfrac  = (float) shr_type_v.at(8) / npixels;
	  protonfrac   = (float) shr_type_v.at(9) / npixels;
	}
      
      } // end this shower
      
      _tree->Fill();
    } // end vertex
    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  
  void ShowerTruthMatch::finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _tree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void ShowerTruthMatch::ClearVertex() {
    _vtxid    = -1.0*larcv::kINVALID_INT;
    _nshowers = -1.0*larcv::kINVALID_INT;

    _npx_v.clear();

    _shr_type_v.clear();
    _shr_type_vv.clear();

    _unknownfrac_v.clear(); 
    _electronfrac_v.clear();
    _gammafrac_v.clear(); 
    _pizerofrac_v.clear(); 
    _muonfrac_v.clear(); 
    _kminusfrac_v.clear(); 
    _piminusfrac_v.clear();
    _protonfrac_v.clear(); 
  }

  void ShowerTruthMatch::ResizeFrac(int nshowers) {
    _unknownfrac_v.resize(nshowers);
    _electronfrac_v.resize(nshowers);
    _gammafrac_v.resize(nshowers);
    _pizerofrac_v.resize(nshowers);
    _muonfrac_v.resize(nshowers);
    _kminusfrac_v.resize(nshowers);
    _piminusfrac_v.resize(nshowers);
    _protonfrac_v.resize(nshowers);
    _npx_v.resize(nshowers);      
    _shr_type_v.resize(nshowers);
    _shr_type_vv.resize(nshowers);
    for(auto& v : _shr_type_vv) 
      v.resize(((int)larcv::kROITypeMax)+1,0);
  }

}

#endif
	
