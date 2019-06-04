#ifndef TRACKPGRAPHMATCH_CXX
#define TRACKPGRAPHMATCH_CXX

// llcv
#include "TrackPGraphMatch.h"
#include "MatchUtils.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPGraph.h"
#include "DataFormat/EventPixel2D.h"

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"

#include <array>
#include <cassert>

namespace llcv {

  void TrackPGraphMatch::configure(const larcv::PSet& cfg) {
    LLCV_DEBUG() << "start" << std::endl;

    _trk_reco_prod = cfg.get<std::string>("TrackRecoProducer");
    _adc_img_prod  = cfg.get<std::string>("ADCImageProducer");
    _pgraph_prod   = cfg.get<std::string>("PGraphProducer");
    _pixel_prod    = cfg.get<std::string>("PixelProducer");

    LLCV_INFO() << "TrackRecoProducer: " << _trk_reco_prod << std::endl;
    LLCV_INFO() << "ADCImageProducer:  " << _adc_img_prod << std::endl;
    LLCV_INFO() << "PGraphProducer:    " << _pgraph_prod << std::endl;
    LLCV_INFO() << "PixelProducer:     " << _pixel_prod << std::endl;

    LLCV_DEBUG() << "end" << std::endl;
  }

  void TrackPGraphMatch::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    _tree = new TTree("TrackPGraphMatch","");
    _tree->Branch("run"    , &_run    , "run/I");
    _tree->Branch("subrun" , &_subrun , "subrun/I");
    _tree->Branch("event"  , &_event  , "event/I");
    _tree->Branch("entry"  , &_entry  , "entry/I");
    _tree->Branch("vtxid"  , &_vtxid  , "vtxid/I");

    _tree->Branch("vtx_x"  , &_vtx_x  , "vtx_x/F");
    _tree->Branch("vtx_y"  , &_vtx_y  , "vtx_y/F");
    _tree->Branch("vtx_z"  , &_vtx_z  , "vtx_z/F");

    _tree->Branch("ntracks", &_ntracks, "ntracks/I");
    _tree->Branch("npars"  , &_npars  , "npars/I");

    _tree->Branch("npx_v"  , &_npx_v);
    _tree->Branch("npts_v" , &_npts_v);
    _tree->Branch("trk_type_v", &_trk_type_v);
    _tree->Branch("trk_type_vv", &_trk_type_vv);

    LLCV_DEBUG() << "end" << std::endl;
  }

  bool TrackPGraphMatch::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {

    assert(!_trk_reco_prod.empty());
    assert(!_adc_img_prod.empty());
    assert(!_pgraph_prod.empty());
    assert(!_pixel_prod.empty());
    
    const auto ev_vertex  = (larlite::event_vertex*)sto.get_data(larlite::data::kVertex,_trk_reco_prod);
    const auto ev_adc_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_adc_img_prod);
    const auto ev_pgraph  = (larcv::EventPGraph*)mgr.get_data(larcv::kProductPGraph,_pgraph_prod);
    const auto ev_pixel   = (larcv::EventPixel2D*)mgr.get_data(larcv::kProductPixel2D,_pixel_prod);
    
    LLCV_DEBUG() << "LC RSE=(" << ev_adc_img->run() << "," << ev_adc_img->subrun() << "," << ev_adc_img->event() << ") " 
		 << "& LL RSE=(" << sto.run_id()  << "," << sto.subrun_id() << "," << sto.event_id() << ")" << std::endl;
    
    LLCV_DEBUG() << "GOT: " << ev_vertex->size() << " vertices" << std::endl;

    assert(ev_vertex->size() == ev_pgraph->PGraphArray().size());
    
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
      meta_v[plane] = ev_adc_img->Image2DArray()[plane].meta();
    
    double xpixel = kINVALID_DOUBLE;
    double ypixel = kINVALID_DOUBLE;

    const auto& pix_m       = ev_pixel->Pixel2DClusterArray();
    const auto& pix_meta_m  = ev_pixel->ClusterMetaArray();

    LLCV_DEBUG() << "FOUND: " << ass_track_vv.size() << " associations" << std::endl;

    std::vector<larcv::Image2D> pgraph_img_v;

    for(size_t img_idx=0; img_idx < meta_v.size(); ++img_idx) {
      const auto& meta = meta_v[img_idx];

      if(pgraph_img_v.size() <= img_idx)
	pgraph_img_v.emplace_back(larcv::Image2D(meta));

      if(pgraph_img_v[img_idx].meta() != meta)
	pgraph_img_v[img_idx] = larcv::Image2D(meta);

      auto& pgraph_img = pgraph_img_v[img_idx];
      pgraph_img.paint(-1);
    }

    for( size_t vtx_id = 0; vtx_id < ass_track_vv.size(); ++vtx_id) {
      ClearVertex();
      const auto& vertex      = ev_vertex->at(vtx_id);
      const auto& ass_track_v = ass_track_vv[vtx_id];

      const auto& pgraph        = ev_pgraph->PGraphArray().at(vtx_id);
      const auto& roi_v         = pgraph.ParticleArray();
      const auto& cluster_idx_v = pgraph.ClusterIndexArray();

      _vtx_x = vertex.X();
      _vtx_y = vertex.Y();
      _vtx_z = vertex.Z();
      
      auto npars = roi_v.size();
      
      for(size_t plane=0; plane<3; ++plane) {
	auto& pgraph_img = pgraph_img_v[plane];
	pgraph_img.paint(-1);
      }
      
      for(size_t roid=0; roid < npars; ++roid) {
	const auto& roi = roi_v[roid];
	auto cidx = cluster_idx_v.at(roid);

	float px_val = (float) roid;

	for(size_t plane=0; plane<3; ++plane) {
	  
	  auto iter_pix = pix_m.find(plane);
	  if(iter_pix == pix_m.end())
	    continue;

	  auto iter_pix_meta = pix_meta_m.find(plane);
	  if(iter_pix_meta == pix_meta_m.end())
	    continue;
	  
	  auto const& pix_v      = (*iter_pix).second;
	  auto const& pix_meta_v = (*iter_pix_meta).second;
	  
	  auto const& pix      = pix_v.at(cidx);
	  auto const& pix_meta = pix_meta_v.at(cidx);
	  
	  auto& pgraph_img = pgraph_img_v[plane];

	  for(const auto& px : pix) {
	    auto posx = pix_meta.pos_x(px.Y());
	    auto posy = pix_meta.pos_y(px.X());
	    auto row  = pgraph_img.meta().row(posy);
	    auto col  = pgraph_img.meta().col(posx);
	    pgraph_img.set_pixel(row,col,px_val);
	  } // end pixel
	} // end plane
      }	// end this ROI

      
      _vtxid   = (int)vtx_id;
      _ntracks = (int)ass_track_v.size();
      _npars   = (int)npars;

      LLCV_DEBUG() << "@vtx_id=" << vtx_id << std::endl;
      LLCV_DEBUG() << "..." << ass_track_vv[vtx_id].size() << " associated tracks" << std::endl;
      LLCV_DEBUG() << "..." << npars << " particles" << std::endl;

      _npx_v.resize(_ntracks);
      _npts_v.resize(_ntracks);
      _trk_type_v.resize(_ntracks);
      _trk_type_vv.resize(_ntracks);

      for(size_t aid=0; aid<ass_track_v.size(); ++aid) {
	auto assid = ass_track_v[aid];
	const auto& track = ev_track->at(assid);

	_npts_v[aid] = track.NumberTrajectoryPoints();
	_npx_v[aid] = 0;

	auto& trk_type_v = _trk_type_vv[aid];
	trk_type_v.clear();
	trk_type_v.resize(npars+1,0);

	LLCV_DEBUG() << "@track=" << aid << " sz=" << track.NumberTrajectoryPoints() << std::endl;

	for(size_t pid=0; pid< track.NumberTrajectoryPoints(); ++pid) {
	  const auto& pt = track.LocationAtPoint(pid);

	  LLCV_DEBUG() << "@pid: "<<pid<<"=(" << pt.X() << ","  << pt.Y() << "," << pt.Z() << ")" << std::endl;

	  for(size_t plane=0; plane<3; ++plane) {
	    LLCV_DEBUG() << "@plane=" << plane << std::endl;
	    const auto& adc_img = ev_adc_img->Image2DArray()[plane];
	    const auto& pgraph_img = pgraph_img_v[plane];
	    xpixel = kINVALID_DOUBLE;
	    ypixel = kINVALID_DOUBLE;

	    Project3D(meta_v[plane],pt.X(),pt.Y(),pt.Z(),0.0,plane,xpixel,ypixel);

	    int xx = (int)(xpixel+0.5);
	    int yy = (int)(ypixel+0.5);
	    yy = pgraph_img.meta().rows() - yy - 1;
	    
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
	    // float pixel_type  = pgraph_img.pixel(yy,xx);

	    float pixel_type = TestPixelType(yy,xx,adc_img,pgraph_img);
	    if (pixel_type < -1) continue;	    

	    LLCV_DEBUG() << "(xx,yy) = " << "(" << yy << "," << xx << ")=" << pixel_type << std::endl;

	    if (pixel_type<0) trk_type_v.back() += 1;
	    else trk_type_v.at((size_t)pixel_type) += 1;

	    _npx_v[aid] += 1;
	  } // end plane
	} // end 3D point

	auto res_itr  = std::max_element(std::begin(trk_type_v), std::end(trk_type_v)-1);
	auto res_loc  = std::distance(std::begin(trk_type_v), res_itr);
	
	_trk_type_v[aid] = (int)res_loc;
	
	LLCV_DEBUG() << "...next track" << std::endl;
	LLCV_DEBUG() << std::endl;
      } // End track
      
	LLCV_DEBUG() << "...next vertex" << std::endl;
	LLCV_DEBUG() << std::endl;
      _tree->Fill();
    } // end vertex
    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  
  void TrackPGraphMatch::finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _tree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }
  
  void TrackPGraphMatch::ClearVertex() {
    _vtxid   = -1.0*larcv::kINVALID_INT;
    _ntracks = -1.0*larcv::kINVALID_INT;
        
    _vtx_x = -1.0*kINVALID_FLOAT;
    _vtx_y = -1.0*kINVALID_FLOAT;
    _vtx_z = -1.0*kINVALID_FLOAT;

    _npx_v.clear();
    _npts_v.clear();
    _trk_type_v.clear();
    _trk_type_vv.clear();
  }



}

#endif
	
