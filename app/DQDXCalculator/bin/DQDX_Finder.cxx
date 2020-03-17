#ifndef DQDX_FINDER_CXX
#define DQDX_FINDER_CXX
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

// ROOT
#include "TH2F.h"
// #include "TGraph.h"
// #include "TCanvas.h"
// #include "TStyle.h"

// // larutil
// #include "LArUtil/LArProperties.h"
// #include "LArUtil/DetectorProperties.h"
// #include "LArUtil/Geometry.h"
// #include "LArUtil/ClockConstants.h"
// #include "LArUtil/SpaceChargeMicroBooNE.h"

// larcv
#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPGraph.h"

// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/track.h"

//DQDXCalculator Package
#include "DQDXCalculator/UtilFunctions.h"
#include "DQDXCalculator/Visualize_Functions.h"
#include "DQDXCalculator/DQDXBuilder.h"

// std::vector<double> calc_dqdx_vector(const larlite::track& reco_track, std::vector<larcv::Image2D> img_v);

// int main(int nargs, char** argv){
int main(){

	std::cout << "Hello world " << "\n";

	std::vector<std::string> files = {
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_00/out/v08_00_00_29/3315480_683/merged_dlreco_11fe5236-0e2b-4644-9092-403daeae8eb1.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_00/out/v08_00_00_29/3315586_772/merged_dlreco_ab3df3ea-80e7-49a4-a2a9-55fd7f42ad98.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_02/out/v08_00_00_29/27789778_225/merged_dlreco_12c2ee92-7818-49ee-b320-20f3e42ec79c.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_02/out/v08_00_00_29/27789967_314/merged_dlreco_0a635ea5-61b3-49ad-a016-e05eb2652657.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_02/out/v08_00_00_29/27791451_805/merged_dlreco_ee3c9895-ca73-4b0c-b00e-4df0e78939cd.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_03/out/v08_00_00_29/3322582_257/merged_dlreco_ea6192de-db64-41cb-b3df-bcc1e7f2a88c.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_03/out/v08_00_00_29/3322813_465/merged_dlreco_f7c4ceb8-9d25-49f6-87ab-cd9fa7812b2c.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_03/out/v08_00_00_29/3323018_648/merged_dlreco_9df8853a-7939-4b00-9937-2fc57f85b7f5.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_03/out/v08_00_00_29/3323036_665/merged_dlreco_672dce3e-76d5-4d49-86f1-6a3a48c3bf8b.root ",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_04/out/v08_00_00_29/28147902_267/merged_dlreco_49d7b5ab-d489-4ab8-8397-34490087dc6c.root",
		"/cluster/tufts/wongjiradlab/jmills09/files_merged/5e19_wctagger_ptest_04/out/v08_00_00_29/28148481_766/merged_dlreco_2fe23124-ef32-49a8-8639-5e1663f33223.root"
	};
	std::vector<int> specific_entry_number = {
		2,
		8,
		7,
		13,
		8,
		4,
		15,
		8,
		7,
		6,
		8
	};
	std::string producer_name_img = "wire";
	std::string producer_name_vtx = "test";
	std::string producer_name_in_track = "trackReco";
	std::string producer_name_out_track = "dqdxTrack";

	std::string output_dir = "p_dqdx_images/";

	larlitecv::DQDXBuilder dqdxbuilder;


	for (unsigned int file_index = 0;file_index<files.size(); file_index++){

	  std::string input_file = files[file_index];
		if (file_index !=5){
			continue;
		}

		int specific_entry = specific_entry_number[file_index];
		// int specific_entry = -1;


		gStyle->SetOptStat(0);
		larcv::IOManager* io_larcv  = new larcv::IOManager(larcv::IOManager::kREAD,"IOManager_Tagger", larcv::IOManager::kTickBackward);
		io_larcv->add_in_file(input_file);
		io_larcv->initialize();

		larlite::storage_manager* io_larlite  = new larlite::storage_manager(larlite::storage_manager::kREAD);
		io_larlite->add_in_filename(input_file);
		io_larlite->open();


		larlite::storage_manager* io_out_ll  = new larlite::storage_manager(larlite::storage_manager::kWRITE);
		io_out_ll->set_out_filename("testout_ll.root");
		io_out_ll->open();


		int start_entry = 0;
		int nentries_mc_cv = io_larcv->get_n_entries();
		std::cout << "Entries in File:  " << nentries_mc_cv << "\n";
		if (specific_entry != -1){
			start_entry = specific_entry;
			nentries_mc_cv = start_entry+1;
		}
		for (int entry=start_entry; entry < nentries_mc_cv; entry++){
		  std::cout << "Entry : " << entry << "\n\n";
			io_larcv->read_entry(entry);
			io_larlite->go_to(entry);
			// LArCV Imports
			std::cout << "\n";
		  larcv::EventImage2D* ev_in_adc_dlreco  = (larcv::EventImage2D*)(io_larcv->get_data(larcv::kProductImage2D, producer_name_img.c_str()));
			larlite::event_track& ev_dl_track       = *((larlite::event_track*)io_larlite->get_data(larlite::data::kTrack,  producer_name_in_track.c_str() ));

			larlite::event_track& ev_output_track       = *((larlite::event_track*)io_out_ll->get_data(larlite::data::kTrack,  producer_name_out_track.c_str() ));

			int run = ev_in_adc_dlreco->run();
			int subrun = ev_in_adc_dlreco->subrun();
			int event = ev_in_adc_dlreco->event();
			io_out_ll->set_id(run,subrun,event);
			std::cout << run << " " << subrun << " " << event << "\n";
		  std::vector< larcv::Image2D > img_v = ev_in_adc_dlreco->Image2DArray();
			larlitecv::make_evdisp_single(img_v[2], "ev_disp_raw_y");
			dqdxbuilder.set_img_v(img_v);
			// Try to instantiate DQDXBuilder class
			int idx1 = 0;
			for (larlite::track reco3d_track : ev_dl_track){
				if (idx1 !=0){continue;}
				std::vector<double> reco_vertex = {(float)reco3d_track.Vertex().X(), (float)reco3d_track.Vertex().Y(), (float)reco3d_track.Vertex().Z()};
				std::vector<int> trk_vtx_rc = larlitecv::getProjectedPixel(reco_vertex, img_v[0].meta(), 3);
				// larlite::geo::View_t views[3] = {larlite::geo::kU,larlite::geo::kV,larlite::geo::kZ};
				larlitecv::make_dqdx_curve(reco3d_track, 2, Form("dqdx_adrien_%d_%d", event,2));

				larlite::track dqdx_track = dqdxbuilder.calc_dqdx_track_revamp(reco3d_track);
				if (dqdx_track.NumberTrajectoryPoints() > 0){
					larlitecv::make_dqdx_curve(dqdx_track, 2, Form("dqdx_minetrack_%d_%d", event,2));
				}
				ev_output_track.push_back(dqdx_track);


				idx1++;

			}//end loop through tracks

			io_out_ll->next_event();
		} //End of entry loop
		io_larcv->finalize();
		io_larlite->close();
		io_out_ll->close();

		delete io_larcv;
		delete io_larlite;
		delete io_out_ll;
	}
	// larlitecv::make_evdisp(dqdxbuilder.residual_dqdx, "Residual_dqdx", "Distance From End (cm)", "dQdX");
	larlitecv::print_signal();


	return 0;
	}//End of main
#endif
