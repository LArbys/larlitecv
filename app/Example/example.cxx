#include <iostream>
#include <cmath>
#include <utility>
#include <ctime>

// config/storage: from LArCV
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h" // example of data product

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"

// ROOT includes
#include "TRandom3.h"


int main( int nargs, char** argv ) {
  
  std::cout << "[BOUNDARY MUON TAGGER]" << std::endl;
  
  // configuration
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );
  larcv::PSet exconfig = cfg.get<larcv::PSet>("ExampleConfigurationFile");
  std::string larcv_image_producer = exconfig.get<std::string>("InputLArCVImages");
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  // larlite
  //dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_mcinfo_0000.root", "larlite" );
  //dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_wire_0000.root", "larlite" );
  //dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opdigit_0000.root", "larlite" );
  //dataco.add_inputfile( "data/data_samples/v05/spoon/larlite/larlite_opreco_0000.root", "larlite" );

  // larcv
  dataco.add_inputfile( "data/data_samples/v05/spoon/larcv/spoon_larcv_out_0000.root", "larcv" );

  // configure
  dataco.configure( "config.cfg", "StorageManager", "IOManager", "ExampleConfigurationFile" );
  
  // initialize
  dataco.initialize();


  // Start Event Loop
  int nentries = dataco.get_nentries("larcv");
  nentries = 3; // to shorten the loop

  TRandom3 rand(time(NULL));
  
  for (int ientry=0; ientry<nentries; ientry++) {
    
    dataco.goto_entry(ientry,"larcv");

    // get images (from larcv)
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, larcv_image_producer );
    std::cout << "get data: number of images=" << event_imgs->Image2DArray().size() << std::endl;
    
    // make a new image that's the same size as the original ones
    std::vector< larcv::Image2D > output_container;
    // we iterate through the images
    for ( auto &img : event_imgs->Image2DArray() ) {
      const larcv::ImageMeta& meta = img.meta(); // get a constant reference to the Image's meta data
      larcv::Image2D newimg( meta );
      for (int irow=0; irow<meta.rows(); irow++) {
	for (int icol=0; icol<meta.cols(); icol++) {
	  newimg.set_pixel( irow, icol, rand.Uniform() );
	}
      }
      output_container.emplace_back( std::move(newimg) ); // we use emplace and move so that we avoid making a copy
    }

    
    // make the output event container
    larcv::EventImage2D* out_event_imgs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "rando" ); // rando matches output list
    out_event_imgs->Emplace( std::move( output_container ) );
    
    // go to tree
    dataco.save_entry();
  }
  
  dataco.finalize();

  return 0;
}
