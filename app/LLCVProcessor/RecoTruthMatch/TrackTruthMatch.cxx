#ifndef TRACKTRUTHMATCH_CXX
#define TRACKTRUTHMATCH_CXX

// llcv
#include "TrackTruthMatch.h"
#include "MatchUtils.h"

// larcv
#include "DataFormat/EventImage2D.h"

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"

#include <array>
#include <cassert>

namespace llcv {

  void TrackTruthMatch::configure(const larcv::PSet& cfg) {
    LLCV_DEBUG() << "start" << std::endl;

    _trk_reco_prod = cfg.get<std::string>("TrackRecoProducer");
    _adc_img_prod  = cfg.get<std::string>("ADCImageProducer");
    _seg_img_prod  = cfg.get<std::string>("TrueImageProducer");


    LLCV_INFO() << "TrackRecoProducer: " << _trk_reco_prod << std::endl;
    LLCV_INFO() << "ADCImageProducer:  " << _adc_img_prod << std::endl;
    LLCV_INFO() << "TrueImageProducer: " << _seg_img_prod << std::endl;

    LLCV_DEBUG() << "end" << std::endl;
  }

  void TrackTruthMatch::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    _tree = new TTree("TrackTruthMatch","");
    _tree->Branch("run"    , &_run    , "run/I");
    _tree->Branch("subrun" , &_subrun , "subrun/I");
    _tree->Branch("event"  , &_event  , "event/I");
    _tree->Branch("entry"  , &_entry  , "entry/I");
    _tree->Branch("vtxid"  , &_vtxid  , "vtxid/I");
    _tree->Branch("ntracks", &_ntracks, "ntracks/I");

    _tree->Branch("unknownfrac_v" , &_unknownfrac_v);
    _tree->Branch("electronfrac_v", &_electronfrac_v);
    _tree->Branch("gammafrac_v"   , &_gammafrac_v);
    _tree->Branch("pizerofrac_v"  , &_pizerofrac_v);
    _tree->Branch("muonfrac_v"    , &_muonfrac_v);
    _tree->Branch("kminusfrac_v"  , &_kminusfrac_v);
    _tree->Branch("piminusfrac_v" , &_piminusfrac_v);
    _tree->Branch("protonfrac_v"  , &_protonfrac_v);

    _tree->Branch("npx_v"  , &_npx_v);
    _tree->Branch("npts_v" , &_npts_v);
    _tree->Branch("trk_type_v", &_trk_type_v);
    _tree->Branch("trk_type_vv", &_trk_type_vv);

    LLCV_DEBUG() << "end" << std::endl;
  }

  bool TrackTruthMatch::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {

    if (_seg_img_prod.empty()) return true;

    assert(!_trk_reco_prod.empty());
    assert(!_adc_img_prod.empty());
    
    const auto ev_vertex  = (larlite::event_vertex*)sto.get_data(larlite::data::kVertex,_trk_reco_prod);
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

    if (ev_vertex->empty()) return true;

    larlite::event_track *ev_track = nullptr;
    auto const& ass_track_vv = sto.find_one_ass(ev_vertex->id(), ev_track, ev_vertex->name());
    if (!ev_track) return true;
    
    std::array<larcv::ImageMeta,3> meta_v;
    for(size_t plane=0; plane<3; ++plane) 
      meta_v[plane] = ev_seg_img->Image2DArray()[plane].meta();
    
    double xpixel=kINVALID_DOUBLE;
    double ypixel=kINVALID_DOUBLE;

    LLCV_DEBUG() << "FOUND: " << ass_track_vv.size() << " associations" << std::endl;
    for( size_t vtx_id = 0; vtx_id < ass_track_vv.size(); ++vtx_id) {
      ClearVertex();
      const auto& vertex = ev_vertex->at(vtx_id);
      const auto& ass_track_v = ass_track_vv[vtx_id];

      _vtxid = vtx_id;
      _ntracks = (int)ass_track_v.size();

      _unknownfrac_v.resize(_ntracks);
      _electronfrac_v.resize(_ntracks);
      _gammafrac_v.resize(_ntracks);
      _pizerofrac_v.resize(_ntracks);
      _muonfrac_v.resize(_ntracks);
      _kminusfrac_v.resize(_ntracks);
      _piminusfrac_v.resize(_ntracks);
      _protonfrac_v.resize(_ntracks);
      _npx_v.resize(_ntracks);      
      _npts_v.resize(_ntracks);      
      _trk_type_v.resize(_ntracks);
      _trk_type_vv.resize(_ntracks);
      


      LLCV_DEBUG() << "@vtx_id=" << vtx_id << std::endl;
      LLCV_DEBUG() << "..." << ass_track_vv[vtx_id].size() << " associated tracks" << std::endl;

      for(size_t aid=0; aid<ass_track_v.size(); ++aid) {
	auto assid = ass_track_v[aid];
	const auto& track = ev_track->at(assid);

	_npts_v[aid] = track.NumberTrajectoryPoints();
	_npx_v[aid] = 0;

	auto& pt_type_v = _trk_type_vv[aid];
	pt_type_v.clear();
	pt_type_v.resize(((int)larcv::kROITypeMax)+1,0);

	LLCV_DEBUG() << "@track=" << aid << " sz=" << track.NumberTrajectoryPoints() << std::endl;

	for(size_t pid=0; pid< track.NumberTrajectoryPoints(); ++pid) {
	  const auto& pt = track.LocationAtPoint(pid);

	  LLCV_DEBUG() << "@pid: "<<pid<<"=(" << pt.X() << ","  << pt.Y() << "," << pt.Z() << ")" << std::endl;

	  for(size_t plane=0; plane<3; ++plane) {
	    const auto& adc_img = ev_adc_img->Image2DArray()[plane];
	    const auto& seg_img = ev_seg_img->Image2DArray()[plane];
	    xpixel = kINVALID_DOUBLE;
	    ypixel = kINVALID_DOUBLE;
	    Project3D(meta_v[plane],pt.X(),pt.Y(),pt.Z(),0.0,plane,xpixel,ypixel);
	    int xx = (int)(xpixel+0.5);
	    int yy = (int)(ypixel+0.5);
	    yy = seg_img.meta().rows() - yy - 1;

	    if (yy >= (int)adc_img.meta().rows() or yy < 0) {
	      LLCV_WARNING() << "3D location outside of image! Skip this plane" << std::endl;
	      LLCV_WARNING() << "(rx,ry,rz)=(" << pt.X() << "," << pt.Y() << "," << pt.Z() << ")" << std::endl;
	      LLCV_WARNING() << "(x,y)=(" << xx << "," << yy << ")" << std::endl;
	      continue;
	    }
	    
	    if (xx >= (int)adc_img.meta().cols() or xx < 0) {
	      LLCV_WARNING() << "3D location outside of image! Skip this plane" << std::endl;
	      LLCV_WARNING() << "(rx,ry,rz)=(" << pt.X() << "," << pt.Y() << "," << pt.Z() << ")" << std::endl;
	      LLCV_WARNING() << "(x,y)=(" << xx << "," << yy << ")" << std::endl;
	      continue;
	    }

	    // float pixel_value = adc_img.pixel(yy,xx);
	    // if(pixel_value==0) continue;
	    // float pixel_type  = seg_img.pixel(yy,xx);

	    float pixel_type = TestPixelType(yy,xx,adc_img,seg_img,true);
	    if (pixel_type < -1) continue;	    	    
	    if (pixel_type < 0) pixel_type = 0;

	    pt_type_v.at((size_t)pixel_type) += 1;
	    _npx_v[aid] += 1;
	    LLCV_DEBUG() << "p:" << plane << "(" << yy << "," << xx << ")=" << (size_t)pixel_type << std::endl;
	  } // end plane
	} // end 3D point

	auto res_iter = std::max_element(std::begin(pt_type_v), std::end(pt_type_v));
	auto res_loc  = std::distance(std::begin(pt_type_v), res_iter);

	auto& unknownfrac  = _unknownfrac_v[aid];
	auto& electronfrac = _electronfrac_v[aid];
	auto& gammafrac    = _gammafrac_v[aid];
	auto& pizerofrac   = _pizerofrac_v[aid];
	auto& muonfrac     = _muonfrac_v[aid];
	auto& kminusfrac   = _kminusfrac_v[aid];
	auto& piminusfrac  = _piminusfrac_v[aid];
	auto& protonfrac   = _protonfrac_v[aid];

	unknownfrac  = 0.0;
	electronfrac = 0.0;
	gammafrac    = 0.0;
	pizerofrac   = 0.0;
	muonfrac     = 0.0;
	kminusfrac   = 0.0;
	piminusfrac  = 0.0;
	protonfrac   = 0.0;
	
	float npixels = (float) _npx_v[aid];
	if (npixels) {
	  unknownfrac  = (float) pt_type_v.at(0) / npixels;
	  electronfrac = (float) pt_type_v.at(3) / npixels;
	  gammafrac    = (float) pt_type_v.at(4) / npixels;
	  pizerofrac   = (float) pt_type_v.at(5) / npixels;
	  muonfrac     = (float) pt_type_v.at(6) / npixels;
	  kminusfrac   = (float) pt_type_v.at(7) / npixels;
	  piminusfrac  = (float) pt_type_v.at(8) / npixels;
	  protonfrac   = (float) pt_type_v.at(9) / npixels;
	}
	_trk_type_v[aid] = (int)res_loc;

      } // End track

      _tree->Fill();
    } // end vertex
    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  
  void TrackTruthMatch::finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _tree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void TrackTruthMatch::ClearVertex() {
    _vtxid   = -1.0*larcv::kINVALID_INT;
    _ntracks = -1.0*larcv::kINVALID_INT;

    _npx_v.clear();
    _npts_v.clear();
    _trk_type_v.clear();
    _trk_type_vv.clear();

    _unknownfrac_v.clear(); 
    _electronfrac_v.clear();
    _gammafrac_v.clear(); 
    _pizerofrac_v.clear(); 
    _muonfrac_v.clear(); 
    _kminusfrac_v.clear(); 
    _piminusfrac_v.clear();
    _protonfrac_v.clear(); 
  }

}

#endif
	
