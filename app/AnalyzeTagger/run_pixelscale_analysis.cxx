#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"

// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/hit.h"
#include "LArUtil/LArProperties.h"

// larcv
#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"

int main( int nargs, char** argv ) {

  gStyle->SetOptStat(0);

  std::string input_larcv = argv[1];
  std::string input_reco2d = argv[2];
  std::string output_hist = argv[3];
  std::string mccversion = argv[4];
  bool make_image = false;

  larcv::IOManager iolarcv( larcv::IOManager::kREAD );
  iolarcv.add_in_file( input_larcv );
  iolarcv.initialize();

  larlite::storage_manager iolarlite( larlite::storage_manager::kREAD );
  iolarlite.add_in_filename( input_reco2d );
  iolarlite.open();

  TFile fout( output_hist.c_str(), "new" ); // do not rewrite

  // output is a histogram of the pixel distribution
  // not saving a tree, because presumably, the number of pixels is too large
  // maybe its not

  TTree tpix("pixamp","PixelAmpTree"); // entry per hit
  float pixamp;
  float hitamp;
  int planeid;
  int wireid;
  tpix.Branch("planeid",&planeid,"planeid/I");
  tpix.Branch("wireid",&wireid,"wireid/I");
  tpix.Branch("pixamp",&pixamp,"pixamp/F");
  tpix.Branch("hitamp",&hitamp,"hitamp/F");

  const float driftv = larutil::LArProperties::GetME()->DriftVelocity();

  int nentries = iolarcv.get_n_entries();
  for (int i=0; i<nentries; i++) {
    iolarcv.read_entry(i);
    iolarlite.go_to(i);

    std::cout << "entry " << i << std::endl;
    
    auto ev_hit = (larlite::event_hit*)iolarlite.get_data( larlite::data::kHit, "gaushit" );
    auto ev_img = (larcv::EventImage2D*)iolarcv.get_data( larcv::kProductImage2D, "wire" );

    std::cout << "  number of hits = " << ev_hit->size() << std::endl;
    std::cout << "  number of images = " << ev_img->Image2DArray().size() << std::endl;

    
    auto const& p2meta = ev_img->Image2DArray().at(2).meta();
    TH2D* hp2 = nullptr;
    if ( make_image ) {
      hp2 = new TH2D("p2","",p2meta.cols(),p2meta.min_x(),p2meta.max_x(),p2meta.rows(),p2meta.min_y(),p2meta.max_y());

      for (int r=0; r<p2meta.rows(); r++) {
	for (int c=0; c<p2meta.cols(); c++) {
	  hp2->SetBinContent( c+1, p2meta.rows()-1-r+1, ev_img->Image2DArray().at(2).pixel(r,c) );
	}
      }
      
      TCanvas c("c","c",1200,600);
      hp2->Draw("colz");
      c.Update();
      c.SaveAs("test_nohits.png");
    }

    TGraph ghit( ev_hit->size() );
    ghit.SetMarkerStyle(20);
    ghit.SetMarkerSize(0.3);
    int ihit=0;
    for ( auto const& hit : *ev_hit ) {
      // for each hit, we get its pixel value at the peak amplitude
      planeid = (int)hit.WireID().Plane;
      wireid  = (int)hit.WireID().Wire;

      hitamp = hit.PeakAmplitude();
      float peak_us = hit.PeakTime(); // usec from trigger?
      float tick = peak_us;
      if ( mccversion=="mcc9" )
	tick += 4800;
      else if ( mccversion=="mcc8" )
	tick += 2400;

      
      auto const& img = ev_img->Image2DArray().at( planeid );
      auto const& meta = img.meta();
      std::cout << "planeid=" << planeid << " peak_us=" << peak_us << " -> tick=" << tick << " meta=[" << meta.min_y() << "," << meta.max_y() << ")" << std::endl;
      
      if ( tick <= meta.min_y() || tick >=meta.max_y() )
	continue;
      
      int row = meta.row( tick );
      int col = meta.col( wireid );
      pixamp = img.pixel( row, col );

      if ( make_image ) {
	if ( planeid==2 ) {
	  ghit.SetPoint( ihit, wireid, tick );
	  ihit++;
	}
      }
      
      tpix.Fill();
    }
    if ( make_image )
      ghit.Set(ihit);

    if ( make_image ) {
      TCanvas c("c","c",1200,600);
      hp2->Draw("colz");
      ghit.Draw("P");
      c.Update();
      c.SaveAs("test.png");
      delete hp2;
      break;
    }
    
  }//end of entry loop

  
  
  fout.Write();
  fout.Close();

  
  return 0;
}
