#ifndef SEARCHDETACHED_CXX
#define SEARCHDETACHED_CXX

// llcv
#include "SearchDetached.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPGraph.h"
#include "DataFormat/EventPixel2D.h"

// laropencvhandl
#include "RecoTruthMatch/MatchUtils.h"

// locv

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/shower.h"

#include <cassert>

namespace llcv {

  void SearchDetached::configure(const larcv::PSet& cfg) {
    LLCV_DEBUG() << "start" << std::endl;

    _adc_img_prod = cfg.get<std::string>("ADCImageProducer");
    _shr_img_prod = cfg.get<std::string>("ShowerImageProducer");

    _pgraph_prod  = cfg.get<std::string>("PGraphProducer");
    _pixel_prod   = cfg.get<std::string>("Pixel2DProducer");

    _out_pgraph_prod = cfg.get<std::string>("OutPGraphProducer");
    _out_pixel_prod  = cfg.get<std::string>("OutPixel2DProducer");

    _search_distance = cfg.get<float>("SearchDistance");
    
    _mask_particles = cfg.get<bool>("MaskParticles");

    _larmkr.Configure(cfg.get<larcv::PSet>("LArbysImageMaker"));

    if(!_algo) throw llcv_err("No algo specified");

    _algo->Configure(cfg.get<larcv::PSet>(_algo->Name()));

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SearchDetached::initialize() {
    LLCV_DEBUG() << "start" << std::endl;
    LLCV_DEBUG() << "end" << std::endl;
  }

  bool SearchDetached::process(larcv::IOManager& mgr, larlite::storage_manager& sto) {
    LLCV_DEBUG() << "start" << std::endl;

    auto ev_adc_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_adc_img_prod);
    auto ev_shr_img = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_shr_img_prod);
    
    auto ev_pgraph  = (larcv::EventPGraph*) mgr.get_data(larcv::kProductPGraph ,_pgraph_prod);
    auto ev_pixel   = (larcv::EventPixel2D*)mgr.get_data(larcv::kProductPixel2D,_pixel_prod);
    
    auto ev_out_pgraph = (larcv::EventPGraph*) mgr.get_data(larcv::kProductPGraph ,_out_pgraph_prod);
    auto ev_out_pixel  = (larcv::EventPixel2D*)mgr.get_data(larcv::kProductPixel2D,_out_pixel_prod);

    auto adc_img_v = ev_adc_img->Image2DArray();
    auto shr_img_v = ev_shr_img->Image2DArray();
    const auto& pgraph_v   = ev_pgraph->PGraphArray();
    const auto& pix_m      = ev_pixel->Pixel2DClusterArray();
    const auto& pix_meta_m = ev_pixel->ClusterMetaArray();

    //
    // mask the particles around the vertex out of the image
    //
    if (_mask_particles) {
      MaskImage(pgraph_v,pix_m,pix_meta_m,adc_img_v);
      MaskImage(pgraph_v,pix_m,pix_meta_m,shr_img_v);
    }

    // 
    // search for a detached candidate per vertex
    //
    std::vector<larcv::Image2D> adc_crop_img_v(3);
    std::vector<larcv::Image2D> shr_crop_img_v(3);

    std::vector<std::pair<int,int> > vtx_pixel_v;

    for(size_t vtxid=0; vtxid < pgraph_v.size(); ++vtxid) {
      
      const auto& pgraph = pgraph_v[vtxid];
      const auto& par = pgraph.ParticleArray().front();
      auto vtx_X = par.X();
      auto vtx_Y = par.Y();
      auto vtx_Z = par.Z();

      LLCV_DEBUG() << "@vtxid=" << vtxid 
		   << " (" << vtx_X << "," << vtx_Y << "," << vtx_Z << ")" << std::endl;
      
      // project 3d vertex into plane & crop to specified dimension
      double xpixel = kINVALID_DOUBLE;
      double ypixel = kINVALID_DOUBLE;

      vtx_pixel_v.clear();
      vtx_pixel_v.resize(3);

      for(size_t plane=0; plane<3; ++plane) {
	xpixel = ypixel = kINVALID_DOUBLE;
	const auto& meta = shr_img_v[plane].meta();
	Project3D(meta,vtx_X,vtx_Y,vtx_Z,0.0,plane,xpixel,ypixel);
	int xx = (int)(xpixel+0.5);
	int yy = (int)(ypixel+0.5);
	yy = meta.rows() - yy - 1;
	vtx_pixel_v[plane] = std::make_pair(yy,xx);
	LLCV_DEBUG() << "@plane=" << plane << " (" << yy << "," << xx << ")" << std::endl;
      }

      // crop the image by defined number of pixels left and right in user config
      for(size_t plane=0; plane<3; ++plane) {
	auto row = vtx_pixel_v[plane].first;
	auto col = vtx_pixel_v[plane].second;

	const auto& adc_img  = adc_img_v[plane];
	const auto& shr_img  = shr_img_v[plane];
	const auto& meta = shr_img.meta();

	LLCV_DEBUG() << "@plane=" << plane << std::endl;
	LLCV_DEBUG() << "(row,col)="<< row << "," << col << ")" << std::endl;
	LLCV_DEBUG() << "search dist=" << _search_distance << std::endl;
	LLCV_DEBUG() << "pixel (width,height)=(" << meta.pixel_width() << "," << meta.pixel_height() << ")" << std::endl;

	double width  = _search_distance*meta.pixel_width();
	double height = _search_distance*meta.pixel_height();

	LLCV_DEBUG() << "(width,height)=("<< width << "," << height << ")" << std::endl;

	size_t row_count = (size_t) _search_distance;
	size_t col_count = (size_t) _search_distance;

	LLCV_DEBUG() << "(rows,cols)=("<< row_count << "," << col_count << ")" << std::endl;
	LLCV_DEBUG() << "origin (x,y)=( " << meta.tl().x << "," << meta.tl().y << ")" << std::endl;

	double origin_x = meta.tl().x + meta.pixel_width() * ((double)col);
	double origin_y = meta.tl().y - meta.pixel_height() * ((double)row);

	LLCV_DEBUG() << "0 (" << origin_x << "," << origin_y << ")" << std::endl;

	origin_x -= width/2.0;
	origin_y += height/2.0;

	LLCV_DEBUG() << "1 (" << origin_x << "," << origin_y << ")" << std::endl;

	LLCV_DEBUG() << "tl: " << meta.tl().x << "," << meta.tl().y << std::endl;
	LLCV_DEBUG() << "tr: " << meta.tr().x << "," << meta.tr().y << std::endl;
	LLCV_DEBUG() << "bl: " << meta.bl().x << "," << meta.bl().y << std::endl;
	LLCV_DEBUG() << "br: " << meta.br().x << "," << meta.br().y << std::endl;
	
	// check if origin is on the edge, move it to the edge if needed
	if (origin_x < meta.tl().x) origin_x = meta.tl().x;
	if (origin_x > meta.tr().x) origin_x = meta.tr().x;

	if (origin_y > meta.tl().y) origin_y = meta.tl().y;
	if (origin_y < meta.br().y) origin_y = meta.br().y;

	// check if origin on the other edge, move if on the edge
	LLCV_DEBUG() << "2 (" << origin_x << "," << origin_y << ")" << std::endl;

	auto max_x = origin_x + width;
	auto min_y = origin_y - height;

	if (max_x > meta.max_x()) {
	  auto dist = meta.max_x() - max_x;
	  origin_x += dist;
	}

	if (min_y < meta.min_y()) {
	  auto dist = meta.min_y() - min_y;
	  origin_y += dist;
	}

	LLCV_DEBUG() << "3 (" << origin_x << "," << origin_y << ")" << std::endl;

	larcv::ImageMeta crop_meta(width,height,
				   row_count,col_count,
				   origin_x,origin_y,
				   plane);
	
	LLCV_DEBUG() << meta.dump();
	LLCV_DEBUG() << crop_meta.dump();

	adc_crop_img_v[plane] = adc_img.crop(crop_meta);
	shr_crop_img_v[plane] = shr_img.crop(crop_meta);
	
	LLCV_DEBUG() << std::endl;
      }
      
      LLCV_DEBUG() << "extract image" << std::endl;
      // convert cropped image2d to cv::Mat
      auto adc_mat_meta_v = _larmkr.ExtractImage(adc_crop_img_v);
      auto shr_mat_meta_v = _larmkr.ExtractImage(shr_crop_img_v);
      LLCV_DEBUG() << "... extracted" << std::endl;

      larocv::data::Vertex3D vtx3d;
      vtx3d.x = vtx_X;
      vtx3d.y = vtx_Y;
      vtx3d.z = vtx_Z;
      
      // search for detached candidates
      LLCV_DEBUG() << "search!" << std::endl;

      auto ret_v = _algo->Search(vtx3d,adc_mat_meta_v,shr_mat_meta_v);

      LLCV_DEBUG() << "ret_v sz=" << ret_v.size() << std::endl;

      // no detached candidates
      if (ret_v.empty()) continue;

      // convert to larcv data products
      FillOutput(ret_v,ev_out_pgraph,ev_out_pixel,adc_crop_img_v);

      LLCV_DEBUG() << "...next vertex" << std::endl;
    }
    
    LLCV_DEBUG() << "end" << std::endl;
    return true;
  }
  

  void SearchDetached::FillOutput(const std::vector<DetachedCandidate>& detached_v,
				  larcv::EventPGraph* ev_pgraph,
				  larcv::EventPixel2D* ev_pixel,
				  const std::vector<larcv::Image2D>& adc_img_v) {
    larcv::PGraph out_pg;

    for(size_t pid=0;pid<detached_v.size();++pid) {

      const auto& detached = detached_v[pid];
      
      larcv::ROI proi;
      
      proi.Position(detached.origin.x,detached.origin.y,detached.origin.z,kINVALID_DOUBLE);
      proi.EndPosition(detached.start.x,detached.start.y,detached.start.z,kINVALID_DOUBLE);
      proi.Shape(larcv::kShapeShower);

      for(size_t plane=0; plane<3; ++plane)
	proi.AppendBB(adc_img_v.at(plane).meta());

      out_pg.Emplace(std::move(proi),pid);

      std::vector<larcv::Pixel2D> ctor_v;

      for(size_t plane=0; plane<3; ++plane) {
	
	const auto& pmeta = adc_img_v[plane].meta();
	
	ctor_v.clear();

	const auto& pctor = detached.CandidateCluster(plane).ctor;
	const auto& img2d = adc_img_v.at(plane);
	
	if(!pctor.empty()) {
	  ctor_v.reserve(pctor.size());
	  
	  // Store contour
	  for(const auto& pt : pctor)  {
	    auto col  = pmeta.rows() - pt.x - 1;
	    auto row  = pt.y;
	    auto gray = 1.0;
	    ctor_v.emplace_back(row,col);
	    ctor_v.back().Intensity(gray);
	  }
	} // empty particle contour on this plane
	  
	larcv::Pixel2DCluster pixctor(std::move(ctor_v));
	ev_pixel->Emplace(plane,std::move(pixctor),pmeta);
      } // 3 planes
    } // end this detached shower

    ev_pgraph->Emplace(std::move(out_pg));
        
    return; 
  }

  void SearchDetached::MaskImage(const std::vector<larcv::PGraph>& pgraph_v,
				 const std::map<larcv::PlaneID_t, std::vector<larcv::Pixel2DCluster> >& pix_m,
				 const std::map<larcv::PlaneID_t, std::vector<larcv::ImageMeta> >&   pix_meta_m,
				 std::vector<larcv::Image2D>& img_v) {
    
    for(const auto& pgraph : pgraph_v) {
      auto const& roi_v         = pgraph.ParticleArray();
      auto const& cluster_idx_v = pgraph.ClusterIndexArray();

      for(size_t roid=0; roid < roi_v.size(); ++roid) {
	
	const auto& roi = roi_v[roid];
	auto cidx = cluster_idx_v.at(roid);

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
	  
	  auto& plane_img = img_v.at(plane);

	  for(const auto& px : pix) {
	    auto posx = pix_meta.pos_x(px.Y());
	    auto posy = pix_meta.pos_y(px.X());
	    auto row = plane_img.meta().row(posy);
	    auto col = plane_img.meta().col(posx);
	    plane_img.set_pixel(row,col,0);
	  }
	} // end plane
      } // end roi
    } // end pgraph
    
  }
    

  
  void SearchDetached::finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    
    LLCV_DEBUG() << "end" << std::endl;
  }

}

#endif
	
