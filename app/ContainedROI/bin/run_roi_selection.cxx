#include <iostream>
#include <string>
#include <vector>

// larcv
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventImage2D.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/chstatus.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "ThruMu/EmptyChannelAlgo.h"
#include "ContainedROIConfig.h"
#include "ContainedROI.h"


int main( int nargs, char** argv ) {

  std::cout << "[[ ROI SELECTION ]]" << std::endl;

  // we need to have a data coordinator for each stage because the number of entries could be different.
  // we'll coordinate by using event,subrun,run information
  std::string data_folder = "~/working/data/larbys/cosmic_tagger_dev/";

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( data_folder+"/output_larcv.root", "larcv" ); // segment image/original image
  dataco_source.add_inputfile( data_folder+"/output_larlite.root", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
  dataco_thrumu.add_inputfile( data_folder+"/output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.add_inputfile( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larcv_p1.root", "larcv" ); //stopmu-tagger output


  dataco_source.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_thrumu.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_stopmu.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;

  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( "containedroi.cfg" );
  larcv::PSet pset = cfg.get<larcv::PSet>("ContainedROI");
  larcv::PSet contained_pset = pset.get<larcv::PSet>("ContainedROIConfig");
  std::vector<std::string> flashproducers = pset.get< std::vector<std::string> >("OpFlashProducers");

  larlitecv::ContainedROIConfig contained_cfg = larlitecv::CreateContainedROIConfig( contained_pset );

  // contained ROI selection algorithm
  larlitecv::ContainedROI algo( contained_cfg );

  // event loop

  int nentries = dataco_stopmu.get_nentries("larcv");
  int start_entry = 0;
  //int end_entry = nentries;
  int end_entry = 10;

  for ( int ientry=start_entry; ientry<end_entry; ientry++ ) {

    // load same event in all three data sets
    dataco_stopmu.goto_entry(ientry,"larcv");
    int run,subrun,event;
    dataco_stopmu.get_id(run,subrun,event);
    std::cout << "entry " << ientry << std::endl;
    std::cout << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    dataco_thrumu.goto_event(run,subrun,event,"larcv");
    dataco_source.goto_event(run,subrun,event,"larcv");
    algo.m_run = run;
    algo.m_subrun = subrun;
    algo.m_event = event;
    algo.m_entry = ientry;

    // get images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"modimgs");
    larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"marked3d");
    larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D,"stopmu");  	

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
    const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();

    // use larlite chstatus to get badch
    larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco_source.get_larlite_data( larlite::data::kChStatus, "chstatus" );
    larlitecv::EmptyChannelAlgo emptyalgo;
    std::cout << "ch status planes: " << ev_status->size() << std::endl;
    std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    std::cout << "number of bad ch imgs: " << badch_v.size() << std::endl;

    // larlite opflash data
    std::vector< larlite::event_opflash* > opflash_containers;
    for ( auto &flashproducer : flashproducers ) {
      std::cout << "search for flash hits from " << flashproducer << ". ";
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco_source.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << " has " << opdata->size() << " flashes" << std::endl;
      opflash_containers.push_back( opdata );
    }

    std::vector<larcv::ROI> untagged_rois = algo.SelectROIs( imgs_v, thrumu_v, stopmu_v, badch_v, opflash_containers );

  }//end of event loop

  return 0;
}
