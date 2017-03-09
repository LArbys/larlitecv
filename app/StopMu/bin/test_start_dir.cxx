#include <iostream>
#include <cmath>
#include <cstdlib>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventChStatus.h"
#include "UBWireTool/UBWireTool.h"
#include "dbscan/DBSCANAlgo.h"

// larlitecv
#include "GapChs/EmptyChannelAlgo.h"

#include "StopMuAlgo.h"

int main( int nargs, char** argv ) {

  std::string cfg_file = "kfstop.cfg";

  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( "output_larcv_testextbnb_001.root", "larcv" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "KFStopMu" );
  dataco.initialize();

  dataco.goto_entry(0, "larcv");

  // PARAMETERS
  float fThreshold = 10.0;
  int rneighbor = 10;
  int cneighbor = 10;

  // for the test, we target a top-passing, stop muon
  // at time around tick 3880
  larcv::EventImage2D* imgs            = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
  larcv::EventPixel2D* top_spacepoints = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "topspacepts" );
  
  // make the bad channel image
  larlitecv::EmptyChannelAlgo emptyalgo;
  larcv::EventChStatus* ev_status      = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
  std::vector< larcv::Image2D > badchimgs = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status ); //< fill this out!
  std::cout << "number of bad ch imgs: " << badchimgs.size() << std::endl;
  
  const std::vector<larcv::Image2D>& img_v = imgs->Image2DArray();
  const larcv::ImageMeta& meta = img_v.at(0).meta();

  int nendpts = top_spacepoints->Pixel2DArray(0).size();
  std::vector<larcv::Pixel2D> start;
  for (int i=0; i<nendpts; i++) {
    const larcv::Pixel2D& pix = top_spacepoints->Pixel2DArray(0).at(i);
    float tick = meta.pos_y( pix.Y() );
    std::cout << "top endpoint #" << i << ": tick=" << tick << std::endl;
    //if ( tick>3850 && tick<3950 ) {
    if ( tick>6900 && tick<6910 ) {
      std::cout << "Found test start point:tick= " << tick << std::endl;
      for (int p=0; p<3; p++) {
	larcv::Pixel2D copy( top_spacepoints->Pixel2DArray(p).at(i) );
	start.emplace_back( copy );
      }
    }
  }
  
  std::cout << "start point: " << start.size() << std::endl;
  
  larlitecv::StopMuAlgo algo;
  algo.setVerbose(1);

  //std::vector< std::array<float,2> > start_dir2d;
  //std::vector< float > start_dir3d;
  //algo.getStartDirection( img_v, badchimgs, start, rneighbor, cneighbor, fThreshold, start_dir2d, start_dir3d );
  
  larlitecv::Vec3DList_t spacepoints;
  std::vector< std::vector<larcv::Pixel2D> > pixellist;
  algo.runTimeStepTracker( img_v, badchimgs, start, pixellist, spacepoints );
  return 0;
}
