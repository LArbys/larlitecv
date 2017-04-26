#include "StopMuFoxTrot.h"

#include "UBWireTool/UBWireTool.h"

#include "ThruMu/BoundaryEndPt.h"

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

    // incorporating tagged information
    // we mask out the charge information, but we also put that information
    // into a badchannel image as well (later if needed)
    
    for ( auto const& startpt : startpts_v ) {

      FoxTrack ft = m_algo->followTrack( img_v, badch_v, thrumu_v, startpt );
      if ( (int)ft.size()<m_config.min_num_steps )
	continue;

      // if long enough, we continue the end of the track to find the end (not done)

      // make a boundary endpt for the end position
      std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( ft.back().pos(), img_v.front().meta(), img_v.size() );
      std::vector<BoundaryEndPt> bendpt_v;
      for ( int p=0; p<(int)img_v.size(); p++ ) {
	//const larcv::ImageMeta& meta = img_v[p].meta();
	BoundaryEndPt bendpt( imgcoords[0], imgcoords[p+1], larlitecv::kUndefined ); // for stopmu end
	bendpt_v.emplace_back( std::move(bendpt) );
      }
      BoundarySpacePoint endpt( larlitecv::kUndefined, std::move(bendpt_v), ft.back().pos()[0], ft.back().pos()[1], ft.back().pos()[2] );

      // pop the end
      ft.pop_back();

      std::vector< std::vector<double> > path3d;
      path3d.reserve( ft.size() );
      
      for ( auto& foxstep : ft ) {
	if ( !foxstep.isgood() )
	  continue;
	std::vector<double> posd(foxstep.pos().size(),0);
	std::cout << "posd: (";
	for (int i=0; i<(int)foxstep.pos().size(); i++) {
	  posd[i] = foxstep.pos()[i];
	  std::cout << posd[i] << " ";
	}
	std::cout << ")" << std::endl;
	path3d.push_back( posd );
      }
      
      BMTrackCluster3D track3d( startpt, endpt, path3d );
      std::cout << "Stopmu foxtrot track has " << path3d.size() << " steps" << std::endl;
      track3d_v.emplace_back( std::move(track3d) );
    }

    return track3d_v;
  }
  
}
