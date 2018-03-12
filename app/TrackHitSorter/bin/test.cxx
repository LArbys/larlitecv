#include <iostream>
#include <string>

// ROOT
#include "TH2D.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TStyle.h"

// larlitecv
#include "Base/DataCoordinator.h"

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

//larcv
#include "DataFormat/EventImage2D.h"

#include "TrackHitSorter.h"
#include "HandScanTable.h"

int main( int nargs, char** argv ) {

  std::string hitfile    = "/home/twongjirad/working/data/larbys/mcc8v6/goodreco1mu1p/goodreco_1mu1p_shower.root";
  std::string trackfile  = "/home/twongjirad/working/data/larbys/mcc8v6/goodreco1mu1p/goodreco_1mu1p_track.root";
  std::string oprecofile = "/home/twongjirad/working/data/larbys/mcc8v6/goodreco1mu1p/goodreco_1mu1p_opreco.root";
  std::string ssnetfile  = "/home/twongjirad/working/data/larbys/mcc8v6/goodreco1mu1p/goodreco_1mu1p_ssnet.root";


  gStyle->SetOptStat(0);
  
  thsort::HandScanTable handscaninfo("goodrecohandscan_1mu1p.tab");
  //int vtxid = handscaninfo.GetVertexID( 5128, 34, 1729 );
  //std::cout << "table test: vtx="  << vtxid << " should equal 2" << std::endl;
  
  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( hitfile,    "larlite" );
  dataco.add_inputfile( trackfile,  "larlite" );
  dataco.add_inputfile( oprecofile, "larlite" );
  dataco.add_inputfile( ssnetfile,  "larcv" );
  dataco.initialize();

  int nentries = dataco.get_nentries("larlite");
  std::cout << "nentries: " << nentries << std::endl;

  // parameters
  const float max_radius = 0.65;
  const float dedx_bin_step  = 0.5;
  const float dedx_bin_width = 1.0;
  bool dump_images     = true;   // dump image of all hits with markers for track locations
  bool draw_track_hits = true;  // dump hits only for those associated to certain vertex. image for each vertex and each plane (i.e. 3 per vertex)
  bool draw_hit_graph  = true;  // draw hits connected in path- and dist- order
  bool plot_dedx = true;
  int colors[4] = { kCyan, kMagenta, kRed };

  
  //for ( int ientry=0; ientry<nentries; ientry++) {
  for ( int ientry=0; ientry<10; ientry++) {
    dataco.goto_entry( ientry, "larlite" );

    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << "[entry " << ientry << "]" << std::endl;

    int run = 0;
    int subrun = 0;
    int event = 0;
    dataco.get_id( run, subrun, event );

    std::cout << "Run, Subrun, Event: " << run << ", " << subrun << ", "  << event << std::endl;

    int goodvtxid = handscaninfo.GetVertexID( run, subrun, event );

    larcv::EventImage2D* ev_img      = (larcv::EventImage2D*) dataco.get_larcv_data( larcv::kProductImage2D, "modimg" );
    larlite::event_vertex* ev_vertex = (larlite::event_vertex*)dataco.get_larlite_data( larlite::data::kVertex, "trackReco" );
    larlite::event_track*  ev_track  = (larlite::event_track*) dataco.get_larlite_data( larlite::data::kTrack,  "trackReco" );
    larlite::event_hit*    ev_hit    = (larlite::event_hit*)   dataco.get_larlite_data( larlite::data::kHit,    "gaushit" );
    const std::vector<larlite::vertex>& vertex_v = *ev_vertex;
    const std::vector<larlite::track>&  track_v  = *ev_track;
    const std::vector<larlite::hit>&    hit_v    = *ev_hit;
    const std::vector<larcv::Image2D>&  img_v    = ev_img->Image2DArray();

    std::cout << "number of vertices: " << vertex_v.size() << std::endl;
    std::cout << "number of tracks: " << track_v.size() << std::endl;
    std::cout << "number of hits: " << hit_v.size() << std::endl;

    int vtxid = handscaninfo.GetVertexID( 5128, 34, 1729 );    

    // start loop over vertex
    int ivertex=-1;
    for (auto const& vtx : vertex_v ) {
      ivertex++;

      //if ( ivertex!=goodvtxid )
      //continue;

      std::vector< std::vector<float> > dedxplots[3];
      std::vector<int> hitmask( hit_v.size(), 1 );
      
      // loop over tracks, keep only those that start from this vertex
      int itrack=-1;
      for ( auto const& track : track_v ) {
	itrack++;
	
	const TVector3& start = track.LocationAtPoint(0);
	float dist2start = (start.X()-vtx.X())*(start.X()-vtx.X());
	dist2start += (start.Y()-vtx.Y())*(start.Y()-vtx.Y());
	dist2start += (start.Z()-vtx.Z())*(start.Z()-vtx.Z());
	dist2start = sqrt(dist2start);

	if ( dist2start>0.1 )
	  continue;

	// associated track. collect hits for it
	thsort::TrackHitSorter algo;
	algo.buildSortedHitList( vtx, track, hit_v, max_radius, hitmask );

	std::vector< std::vector<float> > dedx_per_plane(3);
	algo.getPathBinneddEdx( dedx_bin_step, dedx_bin_width, dedx_per_plane );
	dedxplots[2].push_back( dedx_per_plane[2] ); // saved per track
	
	std::cout << "Hit Sorter" << std::endl;
	std::cout << " hits: p0=" << algo.pathordered[0].size() << " p1=" << algo.pathordered[1].size() << " p2=" << algo.pathordered[2].size() << std::endl;
	algo.dump();

	if ( draw_track_hits ) {
	  const larcv::ImageMeta& meta = img_v.front().meta();
	  const larutil::Geometry*      geo  = larutil::Geometry::GetME();
	  const larutil::LArProperties* larp = larutil::LArProperties::GetME();
	  TH2D* hist2d[3];      
	  for (int plane=0; plane<3; plane++) {
	    // loop through the event images: one per TPC wireplane
	    
	    // dimensions/coordinates of image (all images have the same dimensions)
	    char histname[50];
	    sprintf( histname, "trackhits_p%d", plane );
	    hist2d[plane] = new TH2D( histname, "TPC data;wire number;time tick number", meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );

	    // draw hit.
	    for ( auto const& ho : algo.pathordered[plane] ) {
	      const larlite::hit& hohit = *ho.phit;
	      int wirecol = meta.col( hohit.WireID().Wire );
	      int tickrow = meta.row( 2400+hohit.PeakTime() );
	      hist2d[plane]->SetBinContent( wirecol+1, meta.rows()-tickrow, hohit.Integral() );
	    }
	    
	    TCanvas c("c","c",1800,1200);
	    c.Draw();
	    hist2d[plane]->Draw("COLZ");

	    float vtxtick = vtx.X()/(larp->DriftVelocity()*0.5)+3200;
	    Double_t vtx_xyz[3];
	    vtx.XYZ( vtx_xyz );
	    float vtxwire = geo->NearestWire( vtx_xyz, plane );
	    
	    TMarker mvtx( vtxwire, vtxtick, 24 );
	    //mvtx.SetMarkerStyle(24);
	    mvtx.SetMarkerColor(kRed);
	    mvtx.SetMarkerSize(1.0);
	    mvtx.Draw();
	    
	    c.Update();
	    char zname[100];
	    sprintf( zname, "tracks_entry%03d_vertex%d_track%d_plane%d.png", ientry, ivertex, itrack, plane );
	    c.SaveAs(zname);

	    delete hist2d[plane];
	  }//end of loop over drawn planes
	}//if draw track hits

	if ( draw_hit_graph ) {
	  const larcv::ImageMeta& meta = img_v.front().meta();
	  const larutil::Geometry*      geo  = larutil::Geometry::GetME();
	  const larutil::LArProperties* larp = larutil::LArProperties::GetME();
	  TGraph* gpath[3];
	  TGraph* gdist[3];      
	  for (int plane=0; plane<3; plane++) {
	    // loop through the event images: one per TPC wireplane
	    
	    // dimensions/coordinates of image (all images have the same dimensions)
	    gpath[plane] = new TGraph( algo.pathordered[plane].size() );
	    gdist[plane] = new TGraph( algo.distordered[plane].size() );	    

	    // draw hit. set value by s
	    int ihit=-1;
	    for ( auto const& ho : algo.pathordered[plane] ) {
	      ihit++;
	      const larlite::hit& hohit = *ho.phit;
	      std::cout << "[ihit " << ihit << "] w=" << hohit.WireID().Wire << " tick=" << 2400+hohit.PeakTime() << std::endl;
	      gpath[plane]->SetPoint(ihit, hohit.WireID().Wire, 2400+hohit.PeakTime() );
	    }
	    
	    ihit=-1;
	    for ( auto const& ho : algo.distordered[plane] ) {
	      ihit++;
	      const larlite::hit& hohit = *ho.phit;
	      gdist[plane]->SetPoint(ihit, hohit.WireID().Wire, 2400+hohit.PeakTime() );
	    }

	    gdist[plane]->SetLineColor(kRed);
	    gdist[plane]->SetMarkerStyle(24);
	    gpath[plane]->SetMarkerStyle(24);
	    
	    TCanvas c("c","c",1800,1200);
	    c.Draw();
	    gpath[plane]->Draw("ALP");
	    gdist[plane]->Draw("LPsame");

	    float vtxtick = vtx.X()/(larp->DriftVelocity()*0.5)+3200;
	    Double_t vtx_xyz[3];
	    vtx.XYZ( vtx_xyz );
	    float vtxwire = geo->NearestWire( vtx_xyz, plane );
	    
	    TMarker mvtx( vtxwire, vtxtick, 24 );
	    //mvtx.SetMarkerStyle(24);
	    mvtx.SetMarkerColor(kRed);
	    mvtx.SetMarkerSize(1.0);
	    mvtx.Draw();
	    
	    c.Update();
	    char zname[100];
	    sprintf( zname, "tgraphs_entry%03d_vertex%d_track%d_plane%d.png", ientry, ivertex, itrack, plane );
	    c.SaveAs(zname);
	    
	    delete gpath[plane];
	    delete gdist[plane];	    
	  }//end of loop over drawn planes
	}//if draw track graphs

	if ( plot_dedx ) {

	  TGraph* gdedx[ dedxplots[2].size() ] = {NULL};
	  int itrack=-1;
	  int maxtrack = 0;
	  int maxlength = 0;
	  for (auto const& dedx_v : dedxplots[2]) {
	    itrack++;
	    gdedx[itrack] = new TGraph( dedx_v.size() );
	    if ( maxlength<dedx_v.size() ) {
	      maxtrack = itrack;
	      maxlength = dedx_v.size();
	    }
	    for (int ipt=0; ipt<dedx_v.size(); ipt++) {
	      gdedx[itrack]->SetPoint( ipt, ipt*dedx_bin_step, dedx_v.at(ipt) );
	    }
	    gdedx[itrack]->SetLineColor( colors[itrack%4] );
	    gdedx[itrack]->SetMarkerColor( colors[itrack%4] );	    
	  }

	  TCanvas c("c","c",800,600);
	  gdedx[maxtrack]->Draw("ALP");
	  //gdedx[maxtrack]->GetYaxis()->SetRangeUser(0,1200);
	  gdedx[maxtrack]->GetYaxis()->SetRangeUser(0,10.0);
	  for (int i=0; i<dedxplots[2].size(); i++) {
	    gdedx[i]->SetMarkerStyle(24);	    
	    gdedx[i]->Draw("LP");
	  }
	  char zname[100];
	  sprintf( zname, "dedx_entry%03d_vertex%d.png", ientry, ivertex, 2 );
	  c.SaveAs(zname);

	}//end of if plot dedx
	
      }//end of track loop
    }//end of vertex loop

    if ( dump_images ) {
      const larcv::ImageMeta& meta = img_v.front().meta();
      const larutil::Geometry*      geo  = larutil::Geometry::GetME();
      const larutil::LArProperties* larp = larutil::LArProperties::GetME();	
      
      TH2D* hist2d[3];
      for (int plane=0; plane<3; plane++) {
	// loop through the event images: one per TPC wireplane
	
	// dimensions/coordinates of image (all images have the same dimensions)
	char histname[50];
	sprintf( histname, "hist2d_p%d", plane );
	hist2d[plane] = new TH2D( histname, "TPC data;wire number;time tick number", meta.cols(), meta.min_x(), meta.max_x(), meta.rows(), meta.min_y(), meta.max_y() );
      }
      
      for (auto const& hit : hit_v ) {
	
	int wireno = hit.WireID().Wire;
	int plane  = hit.WireID().Plane;

	// formula for translating x-position into time
	float tick = hit.PeakTime()+2400;
	
	//std::cout << "(" << plane << "," << wireno << "," << tick << "," << hit.PeakTime() << ")" << std::endl;
	
	if ( wireno>meta.min_x() && wireno<meta.max_x()
	     && tick>meta.min_y() && tick<meta.max_y() ) {
	  // translate wire into col and tick into row
	  int wirecol = meta.col( wireno );
	  int tickrow = meta.row( tick );
	  hist2d[plane]->SetBinContent( wirecol+1, meta.rows()-tickrow, hit.Integral() );
	  //hist2d[plane]->SetBinContent( wirecol+1, tickrow, hit.Integral() );
	}
      }//end of loop over hits

      // int itrack=-1;
      // std::vector< TMarker* > pmarkers[3];
      // for (auto const& track : track_v ) {
      // 	itrack++;
	
      // 	int npts = track.NumberTrajectoryPoints();
      // 	for (int ipt=0; ipt<npts; ipt++) {
	  
      // 	  const TVector3& pt = track.LocationAtPoint( ipt );
      // 	  float tick = pt.X()/(larp->DriftVelocity()*0.5)+3200;
      // 	  for (int p=0; p<3; p++) {
      // 	    float wireno = geo->NearestWire( pt, p );
      // 	    TMarker* pmarker = new TMarker( wireno, tick, 20 );
      // 	    pmarker->SetMarkerColor(colors[itrack%3]);
      // 	    pmarker->SetMarkerStyle(24 + itrack%4 );
      // 	    pmarker->SetMarkerSize(0.3);
      // 	    pmarkers[p].push_back( pmarker );
      // 	  }
      // 	}
      // }

      for (int p=0; p<3; p++) {
	char zname[100];
	sprintf( zname, "img_entry%03d_plane%d.png", ientry, p );
	TCanvas canvas("c", "c", 2400,1800);

	canvas.Draw();
	hist2d[p]->SetMaximum(100);
	hist2d[p]->SetMinimum(0);
	hist2d[p]->Draw("COLZ");

	// for ( auto& pmarker : pmarkers[p] )
	//   pmarker->Draw();
	
	canvas.SaveAs( zname );
      }

      for (int p=0; p<3; p++) {
	// for (int i=0; i<pmarkers[p].size(); i++)
	//   delete pmarkers[p][i];
	delete hist2d[p];
      }
    }

    //if ( ientry>=3 )
    //break;
  }//end of loop over entries
  
  return 0;
}
