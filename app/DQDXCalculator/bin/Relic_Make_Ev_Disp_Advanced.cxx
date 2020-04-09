#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

// ROOT
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

// larutil
#include "LArUtil/LArProperties.h"
#include "LArUtil/DetectorProperties.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/ClockConstants.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

// larcv
#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPGraph.h"

// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/track.h"
void print_signal();
std::vector<int> getProjectedPixel(const std::vector<double>& pos3d,
				   const larcv::ImageMeta& meta,
				   const int nplanes,
				   const float fracpixborder=1.5 );
bool is_close_enough(std::vector<double> pt1, std::vector<double> compare_pt, double thresh=0.1);
double distance_between_pt(const std::vector<double>& pt1,const std::vector<double>& pt2);
double distance_between_pt2d(const int x1, const int y1, const int x2, const int y2);
// std::vector<double> calc_dqdx_vector(const larlite::track& reco_track, std::vector<larcv::Image2D> img_v);
double calc_dqdx_for_track(const larlite::track& reco_track, std::vector<larcv::Image2D> img_v);

int main(int nargs, char** argv){
	std::cout << "Hello world " << "\n";

	bool draw_vtx = true;
	bool draw_track = true;



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
	std::vector<double> xvertexvector = {
		194.142258,
		94.242149,
		23.969749,
		205.965363,
		46.043491,
		96.187157,
		175.152603,
		112.214577,
		49.811749,
		36.981663,
		113.615570
	};
	std::vector<double> yvertexvector = {
		-17.649546,
		-10.008907,
		26.537363,
		-86.864594,
		-100.528175,
		69.279793,
		-82.081833,
		-99.075546,
		71.224274,
		-92.300934,
		72.763504
	};
	std::vector<double> zvertexvector = {
		494.549988,
		748.750000,
		818.950012,
		869.299988,
		710.400024,
		674.750000,
		667.950012,
	  775.549988,
		200.149994,
		262.549988,
		656.650024
	};
	std::string producer_name = "wire";
	std::string output_dir = "p_dqdx_images/";
	std::string producer_name_vtx = "test";
	TH2D residual_dqdx =TH2D("Residual_dQdX","Residual_dQdX ",1000,0,50.,20,0,100.);
	TH1D step_distance_h =TH1D("Step_Distance","Step_Distance ",100,0,5.);


	for (int file_index = 0;file_index<files.size(); file_index++){

	  std::string input_file = files[file_index];
		double x_vtx = xvertexvector[file_index];
		double y_vtx = yvertexvector[file_index];
		double z_vtx = zvertexvector[file_index];

		int specific_entry = specific_entry_number[file_index];

		std::string output_custom = "_trackprint_";

		gStyle->SetOptStat(0);
		larcv::IOManager* io_larcv  = new larcv::IOManager(larcv::IOManager::kREAD,"IOManager_Tagger", larcv::IOManager::kTickBackward);
		// io_larcv->reverse_all_products();
		io_larcv->add_in_file(input_file);
		io_larcv->initialize();

		larlite::storage_manager* io_larlite  = new larlite::storage_manager(larlite::storage_manager::kREAD);
		io_larlite->add_in_filename(input_file);
		io_larlite->open();
		int start_entry = 0;
		int nentries_mc_cv = io_larcv->get_n_entries();
		std::cout << "Entries in File:  " << nentries_mc_cv << "\n";

		if (specific_entry != -1){
			start_entry = specific_entry;
			nentries_mc_cv = start_entry+1;
		}
		for (int entry=start_entry; entry < nentries_mc_cv; entry++){
 			TH2D ev_disp_raw_v =TH2D("ev_disp_raw_v","ev_disp_raw_v ",3456,0,3456.,1008,0,1008.);
			TH2D ev_disp_raw_y =TH2D("ev_disp_raw_y","ev_disp_raw_y ",3456,0,3456.,1008,0,1008.);
			TH2D vtx_u =TH2D("vtx_u","vtx_u ",3456,0,3456.,1008,0,1008.);
			TH2D vtx_v =TH2D("vtx_v","vtx_v ",3456,0,3456.,1008,0,1008.);
			TH2D vtx_y =TH2D("vtx_y","vtx_y ",3456,0,3456.,1008,0,1008.);
			// TH2D trk_u =TH2D("trk_u","trk_u ",3456,0,3456.,1008,0,1008.);
			// TH2D trk_v =TH2D("trk_v","trk_v ",3456,0,3456.,1008,0,1008.);
			// TH2D trk_y =TH2D("trk_y","trk_y ",3456,0,3456.,1008,0,1008.);
			int num_points=0;
			std::vector<int> vec_col_u;
			std::vector<int> vec_col_v;
			std::vector<int> vec_col_y;
			std::vector<int> vec_row;

		  std::cout << "Entry : " << entry << "\n\n";
			io_larcv->read_entry(entry);
			io_larlite->go_to(entry);
			// LArCV Imports
			std::cout << "\n";



		  larcv::EventImage2D* ev_in_adc_dlreco  = (larcv::EventImage2D*)(io_larcv->get_data(larcv::kProductImage2D, producer_name.c_str()));
			int run = ev_in_adc_dlreco->run();
			int subrun = ev_in_adc_dlreco->subrun();
			int event = ev_in_adc_dlreco->event();
			std::cout << run << " " << subrun << " " << event << "\n";
		  std::vector< larcv::Image2D > img_2d_in_v = ev_in_adc_dlreco->Image2DArray();
		  for (int col = 0; col<3456;col++){
		    for (int row = 0; row<1008;row++){
		      double val_u = img_2d_in_v[0].pixel(row,col);
		      double val_v = img_2d_in_v[1].pixel(row,col);
		      double val_y = img_2d_in_v[2].pixel(row,col);
					if (val_u > 100.) val_u = 100.0;
					if (val_v > 100.) val_v = 100.0;
					if (val_y > 100.) val_y = 100.0;

		      ev_disp_raw_u.SetBinContent(col,row,val_u);
		      ev_disp_raw_v.SetBinContent(col,row,val_v);
		      ev_disp_raw_y.SetBinContent(col,row,val_y);
		    }
		  }
			if (draw_vtx ){
				vtx_u.SetXTitle("Column (Wire)");
				vtx_u.SetYTitle("Row Stacked (6 Ticks)");
				vtx_u.SetOption("COLZ");
				vtx_u.SetMarkerStyle(kCircle);
				vtx_u.SetMarkerColor(2);
				vtx_u.SetMarkerSize(3);
				vtx_v.SetXTitle("Column (Wire)");
				vtx_v.SetYTitle("Row Stacked (6 Ticks)");
				vtx_v.SetOption("COLZ");
				vtx_v.SetMarkerStyle(kCircle);
				vtx_v.SetMarkerColor(2);
				vtx_v.SetMarkerSize(3);
				vtx_y.SetXTitle("Column (Wire)");
				vtx_y.SetYTitle("Row (6 Ticks)");
				vtx_y.SetOption("COLZ");
				vtx_y.SetMarkerStyle(kCircle);
				vtx_y.SetMarkerColor(2);
				vtx_y.SetMarkerSize(3);
				larcv::EventPGraph* ev_test_pgraph  = (larcv::EventPGraph*)(io_larcv->get_data(larcv::kProductPGraph,producer_name_vtx.c_str()));
				std::vector<larcv::PGraph> test_pgraph_v = ev_test_pgraph->PGraphArray();
				for (auto pgraph : test_pgraph_v) {
					for (larcv::ROI reco_roi : pgraph.ParticleArray()){
						// std::cout << reco_roi.X() << " " << reco_roi.Y() << " " << reco_roi.Z() << " " << reco_roi.T() <<"\n";
						std::vector<double> reco_vertex = {reco_roi.X(), reco_roi.Y(), reco_roi.Z()};
						std::vector<double> target_vertex = {x_vtx,y_vtx,z_vtx};
						// std::cout << "This:           " << reco_vertex[0] << " " << reco_vertex[1] << " " << reco_vertex[2] << "\n";
						// std::cout << "Looking for:    " << x_vtx << " " << y_vtx<< " " << z_vtx << "\n";

						if (false == is_close_enough(reco_vertex,target_vertex)) continue;
						std::cout << "Projecting!\n";
						std::vector<int> vtx_rc = getProjectedPixel(reco_vertex, img_2d_in_v[0].meta(), 3);
						vtx_u.Fill(vtx_rc[1], vtx_rc[0]);
						vtx_v.Fill(vtx_rc[2], vtx_rc[0]);
						vtx_y.Fill(vtx_rc[3], vtx_rc[0]);
					}
				}
			} //End draw vtx
			if (draw_track){
				// trk_u.SetXTitle("Column (Wire)");
				// trk_u.SetYTitle("Row Stacked (6 Ticks)");
				// trk_u.SetOption("COLZ");
				// trk_u.SetMarkerStyle(kCircle);
				// trk_u.SetMarkerColor(6);
				// trk_u.SetMarkerSize(3);
				// trk_v.SetXTitle("Column (Wire)");
				// trk_v.SetYTitle("Row Stacked (6 Ticks)");
				// trk_v.SetOption("COLZ");
				// trk_v.SetMarkerStyle(kCircle);
				// trk_v.SetMarkerColor(6);
				// trk_v.SetMarkerSize(3);
				// trk_y.SetXTitle("Column (Wire)");
				// trk_y.SetYTitle("Row (6 Ticks)");
				// trk_y.SetOption("COLZ");
				// trk_y.SetMarkerStyle(kCircle);
				// trk_y.SetMarkerColor(6);
				// trk_y.SetMarkerSize(3);
				const auto& ev_dl_track       = *((larlite::event_track*)io_larlite->get_data(larlite::data::kTrack,  "trackReco" ));
				int proton_idx = -1;
				int this_idx = 0;
				double shortest_length = 99999;
				for (larlite::track track : ev_dl_track){
					double distance_along_track = 0;
					for (int pt_idx = 0; pt_idx <track.NumberTrajectoryPoints(); pt_idx++){
						TVector3 locpt = track.LocationAtPoint(pt_idx);
						std::vector<double> pt_trck = {locpt.X(), locpt.Y() , locpt.Z() };

						if (pt_idx!=0){
							TVector3 locpt_last = track.LocationAtPoint(pt_idx-1);
							std::vector<double> pt_trck_last = {locpt_last.X(), locpt_last.Y() , locpt_last.Z() };
							distance_along_track += distance_between_pt(pt_trck, pt_trck_last);
						}

					}
					if (distance_along_track < shortest_length) {
						proton_idx = this_idx;
						shortest_length=distance_along_track;
					}
					this_idx++;
				}
				this_idx = 0 ;
				for (larlite::track track : ev_dl_track){
					if (this_idx != proton_idx) {
						this_idx ++;
						continue;
					}
					std::vector<double> reco_vertex = {(float)track.Vertex().X(), (float)track.Vertex().Y(), (float)track.Vertex().Z()};
					std::vector<double> target_vertex = {x_vtx,y_vtx,z_vtx};
					std::vector<int> vtx_rc = getProjectedPixel(target_vertex, img_2d_in_v[0].meta(), 3);

					// std::cout << "This:           " << reco_vertex[0] << " " << reco_vertex[1] << " " << reco_vertex[2] << "\n";
					// std::cout << "Looking for:    " << x_vtx << " " << y_vtx<< " " << z_vtx << "\n";

					if (false == is_close_enough(reco_vertex,target_vertex)) continue;
					std::cout << "TrackPrinting!\n";
					double distance_along_track_u = 0;
					double distance_from_vertex_u = 0;
					double distance_along_track_v = 0;
					double distance_from_vertex_v = 0;
					double distance_along_track_y = 0;
					double distance_from_vertex_y = 0;
					double distance_along_track = 0;
					double distance_from_vertex = 0;
					double num_nodes = 0;
					double tot_dqdx = 0;
					double dist_this_step=0;
					larlite::geo::View_t views[3] = {larlite::geo::kU,larlite::geo::kV,larlite::geo::kZ};

					for (int pt_idx = track.NumberTrajectoryPoints()-1; pt_idx >= 0 ; pt_idx--){
						TVector3 locpt = track.LocationAtPoint(pt_idx);
						std::vector<double> pt_trck = {locpt.X(), locpt.Y() , locpt.Z() };
						std::vector<int> pt_rc = getProjectedPixel(pt_trck, img_2d_in_v[0].meta(), 3);
						vec_col_u.push_back(pt_rc[1]);
						vec_col_v.push_back(pt_rc[2]);
						vec_col_y.push_back(pt_rc[3]);
						vec_row.push_back(pt_rc[0]);
						double this_dqdx_sum_adrien = 0;
						int num_ok_planes =0;
						for (int a=0;a<3;a++){
							double this_dqdx = track.DQdxAtPoint(pt_idx,views[a]);
							if (this_dqdx != 0) {
								this_dqdx_sum_adrien+=this_dqdx;
								num_ok_planes++;
							}
						}
						if (this_dqdx_sum_adrien !=0){
							if (num_ok_planes!=0) this_dqdx_sum_adrien=this_dqdx_sum_adrien*3./num_ok_planes;
							tot_dqdx += this_dqdx_sum_adrien;
							num_nodes++;
						}

						if (pt_idx==track.NumberTrajectoryPoints()-1){
							residual_dqdx.Fill(0.,track.DQdxAtPoint(pt_idx,larlite::geo::kW));
						}
						else {
							TVector3 locpt_last = track.LocationAtPoint(pt_idx+1);
							std::vector<double> pt_trck_last = {locpt_last.X(), locpt_last.Y() , locpt_last.Z() };
							dist_this_step = distance_between_pt(pt_trck, pt_trck_last);
							step_distance_h.Fill(dist_this_step);
							distance_along_track += distance_between_pt(pt_trck, pt_trck_last);
							residual_dqdx.Fill(distance_along_track,track.DQdxAtPoint(pt_idx,larlite::geo::kW));
						}
						if (num_nodes == 0) num_nodes=1;
						std::cout << "Tot DQDX	" << tot_dqdx << "	Nodes:	" << num_nodes << "	Division:	" << tot_dqdx/num_nodes << "	Dist this step:	" << dist_this_step <<  "\n";
						// trk_u.Fill(pt_rc[1], pt_rc[0]);
						// trk_v.Fill(pt_rc[2], pt_rc[0]);
						// trk_y.Fill(pt_rc[3], pt_rc[0]);
					}
					double track_dqdx = calc_dqdx_for_track(track, img_2d_in_v);
					this_idx++;
					std::cout << distance_along_track << " Selected Track Length\n";
				}//end loop through tracks
			}//end draw trk
			// int colu[vec_col_u.size()];
			// int colv[vec_col_v.size()];
			// int coly[vec_col_y.size()];
			// int row[vec_col_y.size()];
			TGraph graph_u(vec_row.size(), vec_col_u.data(), vec_row.data());
			TGraph graph_v(vec_row.size(), vec_col_v.data(), vec_row.data());
			TGraph graph_y(vec_row.size(), vec_col_y.data(), vec_row.data());


		  TCanvas can("can", "histograms ", 3456, 1008);
		  can.cd();
		  ev_disp_raw_u.SetTitle(Form("Image Raw U Plane Run: %d Subrun: %d Event: %d",run,subrun,event));
		  ev_disp_raw_u.SetXTitle("Column (Wire)");
		  ev_disp_raw_u.SetYTitle("Row (6 Ticks)");
		  ev_disp_raw_u.SetOption("COLZ");
		  ev_disp_raw_u.Draw("");
			if (draw_vtx ) vtx_u.Draw("SAME");
			graph_u.Draw("SAMEL");
			can.SaveAs(Form("%soutput_%d_%d_%d_u%s.png",output_dir.c_str(),run,subrun,event,output_custom.c_str()));
		  ev_disp_raw_v.SetTitle(Form("Image Raw V Plane Run: %d Subrun: %d Event: %d",run,subrun,event));
		  ev_disp_raw_v.SetXTitle("Column (Wire)");
		  ev_disp_raw_v.SetYTitle("Row (6 Ticks)");
		  ev_disp_raw_v.SetOption("COLZ");
		  ev_disp_raw_v.Draw("");
			if (draw_vtx ) vtx_v.Draw("SAME");
			graph_v.Draw("SAMEL");

			can.SaveAs(Form("%soutput_%d_%d_%d_v%s.png",output_dir.c_str(),run,subrun,event,output_custom.c_str()));
		  ev_disp_raw_y.SetTitle(Form("Image Raw Y Plane Run: %d Subrun: %d Event: %d",run,subrun,event));
		  ev_disp_raw_y.SetXTitle("Column (Wire)");
		  ev_disp_raw_y.SetYTitle("Row (6 Ticks)");
		  ev_disp_raw_y.SetOption("COLZ");
		  ev_disp_raw_y.Draw("");
			if (draw_vtx ) vtx_y.Draw("SAME");
			graph_y.Draw("SAMEL");

		  can.SaveAs(Form("%soutput_%d_%d_%d_y%s.png",output_dir.c_str(),run,subrun,event,output_custom.c_str()));

		} //End of entry loop
		io_larcv->finalize();
		delete io_larcv;
	}
	print_signal();
	// TCanvas can2("can2", "histograms2 ", 3000, 1000);
	// residual_dqdx.SetTitle("Residual DQdx");
	// residual_dqdx.SetXTitle("Distance from End");
	// residual_dqdx.SetYTitle("dQdX Yplane");
	// residual_dqdx.SetOption("COLZ");
	// residual_dqdx.Draw();
	// can2.SaveAs(Form("%sResidual_End_Thin_dQdX_yplane.png",output_dir.c_str()));
	// residual_dqdx.ProjectionX()->Draw();
	// can2.SaveAs(Form("%sResidual_End_Thin_Projection_X_yplane.png",output_dir.c_str()));
	// step_distance_h.SetTitle("Step Distance");
	// step_distance_h.SetXTitle("Step Distance 3D (cm)");
	// step_distance_h.SetYTitle("Count");
	// step_distance_h.Draw();
	// can2.SaveAs(Form("%sStep_Distance.png",output_dir.c_str()));

	return 0;
	}//End of main

bool is_close_enough(std::vector<double> pt1, std::vector<double> compare_pt, double thresh){
	bool isclose=true;
	for (int p=0;p<3;p++){
		if ((pt1[p] > compare_pt[p]+thresh) || (pt1[p] < compare_pt[p]-thresh)) isclose=false;
	}
	return isclose;
}
std::vector<int> getProjectedPixel( const std::vector<double>& pos3d,
				    const larcv::ImageMeta& meta,
				    const int nplanes,
				    const float fracpixborder ) {
  std::vector<int> img_coords( nplanes+1, -1 );
  float row_border = fabs(fracpixborder)*meta.pixel_height();
  float col_border = fabs(fracpixborder)*meta.pixel_width();

  // tick/row
  float tick = pos3d[0]/(::larutil::LArProperties::GetME()->DriftVelocity()*::larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3) + 3200.0;
  if ( tick<meta.min_y() ) {
    if ( tick>meta.min_y()-row_border )
      // below min_y-border, out of image
      img_coords[0] = meta.rows()-1; // note that tick axis and row indicies are in inverse order (same order in larcv2)
    else
      // outside of image and border
      img_coords[0] = -1;
  }
  else if ( tick>meta.max_y() ) {
    if ( tick<meta.max_y()+row_border )
      // within upper border
      img_coords[0] = 0;
    else
      // outside of image and border
      img_coords[0] = -1;
  }
  else {
    // within the image
    img_coords[0] = meta.row( tick );
  }

  // Columns
  Double_t xyz[3] = { pos3d[0], pos3d[1], pos3d[2] };
  // there is a corner where the V plane wire number causes an error
  if ( (pos3d[1]>-117.0 && pos3d[1]<-116.0) && pos3d[2]<2.0 ) {
    xyz[1] = -116.0;
  }
  for (int p=0; p<nplanes; p++) {
    float wire = larutil::Geometry::GetME()->WireCoordinate( xyz, p );

    // get image coordinates
    if ( wire<meta.min_x() ) {
      if ( wire>meta.min_x()-col_border ) {
	// within lower border
	img_coords[p+1] = 0;
      }
      else
	img_coords[p+1] = -1;
    }
    else if ( wire>=meta.max_x() ) {
      if ( wire<meta.max_x()+col_border ) {
	// within border
	img_coords[p+1] = meta.cols()-1;
      }
      else
	// outside border
	img_coords[p+1] = -1;
    }
    else
      // inside image
      img_coords[p+1] = meta.col( wire );
  }//end of plane loop

  // there is a corner where the V plane wire number causes an error
  if ( pos3d[1]<-116.3 && pos3d[2]<2.0 && img_coords[1+1]==-1 ) {
    img_coords[1+1] = 0;
  }
  return img_coords;
}

double distance_between_pt(const std::vector<double>& pt1,const std::vector<double>& pt2){
	/*
	This function finds the distance between pt1 and pt2
	which should both be in x,y,z
	*/

  double distance = -1;
  distance = std::sqrt(std::pow(pt1[0]-pt2[0],2)+std::pow(pt1[1]-pt2[1],2)+std::pow(pt1[2]-pt2[2],2));
  return distance ;
}
double distance_between_pt2d(const int x1, const int y1, const int x2, const int y2){
	/*
	This function finds the distance between pt1 and pt2
	which should both be in x,y
	*/
	double xd1 = (double)x1;
	double yd1 = (double)y1;
	double xd2 = (double)x2;
	double yd2 = (double)y2;
	double distance = -1;
	distance = std::sqrt(std::pow(xd1-xd2,2)+std::pow(yd1-yd2,2));
	return distance ;

}

void print_signal(){
	std::cout << "\n\n";
	std::cout << "/////////////////////////////////////////////////////////////\n";
	std::cout << "/////////////////////////////////////////////////////////////\n";
	std::cout << "\n\n";

	std::cout << "         _==/            i     i           \\==_ \n";
	std::cout << "        /XX/             |\\___/|            \\XX\\    \n";
	std::cout << "       /XXXX\\            |XXXXX|            /XXXX\\   \n";
	std::cout << "      |XXXXXX\\_         _XXXXXXX_         _/XXXXXX|   \n";
	std::cout << "     XXXXXXXXXXXxxxxxxxXXXXXXXXXXXxxxxxxxXXXXXXXXXXX   \n";
	std::cout << "    |XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX|  \n";
	std::cout << "    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   \n";
	std::cout << "    |XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX|   \n";
	std::cout << "     XXXXXX/^^^^^\\XXXXXXXXXXXXXXXXXXXXX/^^^^^\\XXXXXX    \n";
	std::cout << "      |XXX|       \\XXX/^^\\XXXXX/^^\\XXX/       |XXX|    \n";
	std::cout << "       \\XX\\        \\X/    \\XXX/    \\X/       /XX/    \n";
	std::cout << "           \\        |      \\X/      |       /     \n";
	std::cout << "\n\n";
	std::cout << "/////////////////////////////////////////////////////////////\n";
	std::cout << "/////////////////////////////////////////////////////////////\n";
	std::cout << "\n\n";
	return;
	}

	// std::vector<double> calc_dqdx_vector(const larlite::track& reco_track, std::vector<larcv::Image2D> img_v){
	// 	std::cout << "Only doing Plane y for now\n";
	// 	std::vector<double> dqdx_vector(reco_track.NumberTrajectoryPoints(), -1);
	// 	double total_q_charge = 0; //Keep track of entire track's total charge
	//
	// 	for (int trkpt_idx = 0; trkpt_idx < reco_track.NumberTrajectoryPoints()-1; trkpt_idx++){
	// 		TVector3 pt1 = reco_track.LocationAtPoint(trkpt_idx);
	// 		TVector3 pt2 = reco_track.LocationAtPoint(trkpt_idx+1);
	//
	// 		std::vector<double> xyz_pt1 = {pt1.X(),pt1.Y(),pt1.Z()};
	// 		std::vector<double> xyz_pt2 = {pt2.X(),pt2.Y(),pt2.Z()};
	// 		double this_dqdx_dx = distance_between_pt(xyz_pt1,xyz_pt2);
	// 		double this_dqdx_dq = 0;
	// 		std::vector<int> rc_pt1 = getProjectedPixel(xyz_pt1, img_v[0].meta(), 3);
	// 		std::vector<int> rc_pt2 = getProjectedPixel(xyz_pt2, img_v[0].meta(), 3)
	// 		double slope_dy = rc_pt2[0] - rc_pt1[0];
	// 		double slope_dx = rc_pt2[3] - rc_pt1[3];
	//
	// 		if (slope_dx == 0){
	// 			// Handle Perfect Vertical Step
	// 			for (int y_coord = rc_pt1[0]; y_coord < rc_pt2[0]; y_coord++){
	// 				this_dqdx_dq
	// 			} //End loop along y_coords
	// 		} // End Edge Case of Perfect Vertical Step
	// 		else {
	// 			double slope = slope_dy/slope_dx;
	// 			double y_cept = rc_pt2[0] - slope*rc_pt2[3];
	// 			for (int x_coord = rc_pt1[3]; x_coord < rc_pt2[3]; x_coord++){
	//
	// 			} // End loop along x coord of line segment
	// 		} // End normal case of line_segment
	// 	} // End loop through track steps
	//
	// 	return dqdx_vector;
	// }

double calc_dqdx_for_track(const larlite::track& reco_track, std::vector<larcv::Image2D> img_v){
		std::cout << "Only doing Plane y for now\n";
		double total_q_charge = 0; //Keep track of entire track's total charge
		double total_length = 0;
		for (int trkpt_idx = 0; trkpt_idx < reco_track.NumberTrajectoryPoints()-1; trkpt_idx++){
					TVector3 pt1 = reco_track.LocationAtPoint(trkpt_idx);
					TVector3 pt2 = reco_track.LocationAtPoint(trkpt_idx+1);

					std::vector<double> xyz_pt1 = {pt1.X(),pt1.Y(),pt1.Z()};
					std::vector<double> xyz_pt2 = {pt2.X(),pt2.Y(),pt2.Z()};
			double this_dqdx_dx = distance_between_pt(xyz_pt1,xyz_pt2);
			total_length += this_dqdx_dx;
			double this_dqdx_dq = 0;
			std::vector<int> rc_pt1 = getProjectedPixel(xyz_pt1, img_v[0].meta(), 3);
			std::vector<int> rc_pt2 = getProjectedPixel(xyz_pt2, img_v[0].meta(), 3);
			double slope_dy = rc_pt2[0] - rc_pt1[0];
			double slope_dx = rc_pt2[3] - rc_pt1[3];

			if (slope_dx == 0){
				// Handle Perfect Vertical Step
				for (int y_coord = rc_pt1[0]; y_coord <= rc_pt2[0]; y_coord++){
					for(int x_coord = rc_pt2[3]-3; x_coord <= rc_pt2[3]+3; x_coord++){
						total_q_charge += img_v[2].pixel(y_coord,x_coord);
						this_dqdx_dq += img_v[2].pixel(y_coord,x_coord);
						img_v[2].set_pixel(y_coord,x_coord,0.0) ; //Ensure no double counting
					}// End loop along x_coords
				} //End loop along y_coords
			} // End Edge Case of Perfect Vertical Step
			else {
				double slope = slope_dy/slope_dx;
				double y_cept = rc_pt2[0] - slope*rc_pt2[3];
				for (int x_coord = rc_pt1[3]-3; x_coord <= rc_pt2[3]+3; x_coord++){
					int y_step = int( x_coord*slope + y_cept + 0.5); //The 0.5 is to help round from double to int. Casting to int always rounds down
					for (int y_coord = y_step-3; y_coord <= y_step+3; y_coord++){
						total_q_charge += img_v[2].pixel(y_coord,x_coord);
						this_dqdx_dq += img_v[2].pixel(y_coord,x_coord);
						img_v[2].set_pixel(y_coord,x_coord,0.0) ; //Ensure no double counting
					}
				} // End loop along x coord of line segment
			} // End normal case of line_segment
			double this_dqdx = this_dqdx_dq/this_dqdx_dx;
			std::cout << "nStep	"<< trkpt_idx << "	Length " << this_dqdx_dx << "	My dQ	" << this_dqdx_dq<< "	My dQdX	:	" << this_dqdx << "	Stored dQdX	" << reco_track.DQdxAtPoint(trkpt_idx,larlite::geo::kZ) << "\n";
		} // End loop through track steps

		double track_dqdx = total_q_charge/total_length;
		std::cout << "Track dQdX	" << track_dqdx << "	Total Charge	" << total_q_charge << "	Total Length	" << total_length << "\n";
		assert (1==2);
		return track_dqdx;
	}
