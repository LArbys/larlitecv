#include "StopMuFoxTrot.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "UBWireTool/UBWireTool.h"

#include "TaggerTypes/BoundaryEndPt.h"

namespace larlitecv {

  StopMuFoxTrot::StopMuFoxTrot( const StopMuFoxTrotConfig& config )
    : m_config(config)  {
    m_algo = new FoxTrotTrackerAlgo( m_config.foxtrotalgo_cfg );
  }

  StopMuFoxTrot::~StopMuFoxTrot() {
    delete m_algo;
  }

  std::vector< larlitecv::BMTrackCluster3D> StopMuFoxTrot::findStopMuTracks( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
									     const std::vector<larcv::Image2D>& thrumu_v, const std::vector<BoundarySpacePoint>& startpts_v ) {
    std::vector< larlitecv::BMTrackCluster3D > track3d_v;
    const float cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    const float cm_per_wire = 0.3; // 0.3 cm per wire on each plane

    // incorporating tagged information
    // we mask out the charge information, but we also put that information
    // into a badchannel image as well (later if needed)
    if ( m_config.verbosity>0 ) {
      std::cout << __FILE__ << ":" << __LINE__ << " searching for StopMu tracks using " << startpts_v.size() << " starting points." << std::endl;
    }

    for ( auto const& startpt : startpts_v ) {

      if ( m_config.verbosity>0 ) {
        std::vector<int> startcoords = larcv::UBWireTool::getProjectedImagePixel( startpt.pos(), img_v.front().meta(), img_v.size() );
        std::cout << __FILE__ << ":" << __LINE__ << " starting point"
                  << " (" << startpt.pos()[0] << "," << startpt.pos()[1] << "," << startpt.pos()[2] << ")"
                  << " imgcoords=(" << startcoords[0] << "," << startcoords[1] << "," << startcoords[2] << "," << startcoords[3] << ")"
                  << std::endl;
      }

      // by request from janet/jarrett, we skip tracking for anode/cathode boundary pts
      if ( m_config.SkipAnodeCathodeStartPts && (startpt.type()==larlitecv::kAnode || startpt.type()==larlitecv::kCathode) ) {
	std::cout << __FILE__ << ":" << __LINE__ << " by request, skipping anode and cathode start points" << std::endl;
	continue;
      }
      
      FoxTrack ft = m_algo->followTrack( img_v, badch_v, thrumu_v, startpt );
      if ( (int)ft.size()<m_config.min_num_steps ) {
        if ( m_config.verbosity>0)
          std::cout << __FILE__ << ":" << __LINE__ << " track too short " << ft.size() << " < " << m_config.min_num_steps << std::endl;
        continue;
      }

      // if long enough, we continue the end of the track to find the end (not done)

      // make a boundary endpt for the end position
      std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( ft.back().pos(), img_v.front().meta(), img_v.size() );

      // correct rounding error for top of image (x=292 or so)
      if ( imgcoords[0]==-1 ) {
	if ( fabs(3200.0+ft.back().pos()[0]/cm_per_tick-img_v.front().meta().max_y())<img_v.front().meta().pixel_height() ) {
	  imgcoords[0] = 0;
	  std::vector<float> backend = ft.back().pos();
	  backend[0] = img_v.front().meta().pos_y(0);
	  std::vector<float> backdir = ft.back().dir();
	  ft.pop_back();
	  FoxStep fixedback( backend, backdir );
	  ft.push_back( fixedback );
	}
	else {
	  std::cout << "row=-1 and dtick=" << fabs(3200.0+ft.back().pos()[0]/cm_per_tick-img_v.front().meta().max_y()) << std::endl;
	  continue;
	}
	
      }

      // Set this variable up here so that it is recognized when it is referenced.
      bool pixel_outside_of_image = false;
      
      // correct the rounding error for the columns of the image.
      // This is the analog of what is done above.
      for (size_t p = 0; p < img_v.size(); p++) {
	if (imgcoords[p+1] == -1) {

	  // Make the value that you need to feed into the WireCoordinate function.
	  Double_t xyz[3] = { ft.back().pos()[0], ft.back().pos()[1], ft.back().pos()[2] };
	  // Convert the y and z coordinates into a wire on this plane.
	  float wid = larutil::Geometry::GetME()->WireCoordinate( xyz, p );

	  // The compression factor is 0 here and the minimum x value is 0, so you can just compare to the pixel width.
	  if ( wid<img_v.front().meta().pixel_width() ) {

	    std::cout << "Wire ID: " << wid << std::endl;
	    std::cout << "Pixel Width: " << img_v.front().meta().pixel_width() << std::endl;
	    
	    imgcoords[p+1] = 0;
	    std::vector<float> backend = ft.back().pos();
	    backend[p] = img_v.front().meta().pos_x(p);
	    std::vector<float> backdir = ft.back().dir();
	    ft.pop_back();
	    FoxStep fixedback( backend, backdir );
	    ft.push_back( fixedback );

	  }
	  else {
	    std::cout << "col=-1 and dwire=" << wid << std::endl;
	    pixel_outside_of_image = true;

	    break;
	  }

	}
      }

      // If the pixel is outside of the image, then continue.
      // I do this here so that we do not interfere with that loop.
      if (pixel_outside_of_image) continue;
      
      std::vector<BoundaryEndPt> bendpt_v;
      for ( int p=0; p<(int)img_v.size(); p++ ) {
        //const larcv::ImageMeta& meta = img_v[p].meta();
        BoundaryEndPt bendpt( imgcoords[0], imgcoords[p+1], larlitecv::kUndefined ); // for stopmu end
        bendpt_v.emplace_back( std::move(bendpt) );
      }
      BoundarySpacePoint endpt( larlitecv::kUndefined, std::move(bendpt_v), ft.back().pos()[0], ft.back().pos()[1], ft.back().pos()[2] );
      if ( m_config.verbosity>0) {
          std::cout << __FILE__ << ":" << __LINE__ << " stopmu endpoint "
                    << " (" << endpt.pos()[0] << "," << endpt.pos()[1] << "," << endpt.pos()[2] << ")" << std::endl;
      }

      // pop the end
      //ft.pop_back();

      std::vector< std::vector<double> > path3d;
      path3d.reserve( ft.size() );

      for ( auto& foxstep : ft ) {
        if ( !foxstep.isgood() )
          continue;
        std::vector<double> posd(foxstep.pos().size(),0);
        //std::cout << "posd: (";
        for (int i=0; i<(int)foxstep.pos().size(); i++) {
          posd[i] = foxstep.pos()[i];
          //std::cout << posd[i] << " ";
        }
        //std::cout << ")" << std::endl;
        path3d.push_back( posd );
      }

      BMTrackCluster3D track3d( startpt, endpt, path3d );
      //std::cout << "Stopmu foxtrot track has " << path3d.size() << " steps" << std::endl;
      track3d_v.emplace_back( std::move(track3d) );
    }

    return track3d_v;
  }

}
