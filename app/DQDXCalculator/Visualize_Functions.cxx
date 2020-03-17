#ifndef VISUALIZE_FUNCTIONS_CXX
#define VISUALIZE_FUNCTIONS_CXX
#include "Visualize_Functions.h"

namespace larlitecv {

void print_tvec(const TVector3 tvec){
  std::cout << tvec.x() << "  "<< tvec.y() << " " << tvec.z() << "\n";
}

void make_evdisp(TH2D& hist, std::string title, std::string x_title, std::string y_title, double maxval, TH2D* vtx_h){
  int xdim = hist.GetNbinsX();
  int ydim = hist.GetNbinsY();
  while ((xdim < 500) || (ydim < 500)){
    xdim *=2;
    ydim *=2;
  }
  TCanvas can("can", "histograms ", xdim, ydim);
  can.cd();
  hist.SetTitle(title.c_str());
  hist.SetXTitle(x_title.c_str());
  hist.SetYTitle(y_title.c_str());
  hist.GetXaxis()->SetNdivisions(hist.GetNbinsX()/4);
  hist.SetOption("COLZ");
  if (maxval > 0){
    hist.SetMaximum(maxval);
  }
  hist.Draw("");
  if (vtx_h != NULL) vtx_h->Draw("SAME");
  can.SaveAs(Form("%s.png", title.c_str() ));
  return;
}
void make_evdisp_single(const larcv::Image2D& img, std::string title, std::string x_title, std::string y_title, double maxval, TH2D* vtx_h){
  TH2D hist = TH2D("hist","hist ",3456,0,3456.,1008,0,1008.);
  for (int r=0;r<1008;r++){
    for (int c=0;c<3456;c++){
      hist.SetBinContent(c,r,img.pixel(r,c));
    }
  }
  make_evdisp(hist, title,x_title,y_title,maxval,vtx_h);
  return;
}

void make_evdisp_triple(const std::vector<larcv::Image2D>& img_v, std::string title, std::string x_title, std::string y_title, double maxval, TH2D* vtx_h){
  TH2D hist = TH2D("hist","hist ",3456,0,3456.,3024,0,3024.);
  for (int r=0;r<1008;r++){
    for (int c=0;c<3456;c++){
      hist.SetBinContent(c,r+2016,img_v[0].pixel(r,c));
      hist.SetBinContent(c,r+1008,img_v[1].pixel(r,c));
      hist.SetBinContent(c,r,img_v[2].pixel(r,c));
    }
  }
  make_evdisp(hist, title,x_title,y_title,maxval,vtx_h);
  return;
}

void make_dqdx_curve(const std::vector<double>& dqdx_v, const std::vector<double>& distance_v, std::string filename){
  if (dqdx_v.size() != distance_v.size()){
    std::cout << "Input Vectors not the same size!\n";
    return;
  }
  TGraph dqdx_curve_g = TGraph(dqdx_v.size(), distance_v.data(), dqdx_v.data());
  TCanvas can("can", "histograms ", 2000, 1000);
  can.cd();
  std::string title = filename+";Distance from End (cm);dQdX";
  dqdx_curve_g.SetTitle(title.c_str());

  dqdx_curve_g.Draw("");
  can.SaveAs(Form("%s.png", filename.c_str() ));
  return;
}
void make_dqdx_curve(const larlite::track& track, int plane, std::string filename){
  larlite::geo::View_t views[3] = {larlite::geo::kU,larlite::geo::kV,larlite::geo::kZ};
  std::vector<double> dqdx_v;
  std::vector<double> dist_v;
  dqdx_v.reserve(track.NumberTrajectoryPoints());
  dist_v.reserve(track.NumberTrajectoryPoints());
  double distance_sum = 0;
  // Iterate backwards so that 0 length is the track end
  for ( int pt_idx = (int)track.NumberTrajectoryPoints()-1;pt_idx>=0;pt_idx--){
    if (pt_idx != (int)track.NumberTrajectoryPoints()-1){
      TVector3 this_pt = track.LocationAtPoint(pt_idx);
      TVector3 last_pt = track.LocationAtPoint(pt_idx+1);
      distance_sum += larlitecv::distance_between_pt(this_pt,last_pt);
      dist_v.push_back(distance_sum);
      dqdx_v.push_back(track.DQdxAtPoint(pt_idx,views[plane]));
    }
  }
  larlitecv::make_dqdx_curve(dqdx_v,dist_v,filename);
  return;
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
}
#endif
