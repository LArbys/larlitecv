#include "NumuFilter.h"

#include "TTreeReader.h"
#include "TFile.h"

namespace larlitecv {
namespace dllee {

  void NumuFilter::ReturnFilteredDictionary(std::string inputfile,
                                            std::map<std::tuple<int,int,int>,bool>& rse,
                                            std::map<std::tuple<int,int,int,int>,bool>& rsev,
                                            std::vector< std::vector<float> >& cutvars,                                            
                                            int verbosity ) {

    TFile* infile  = new TFile(inputfile.c_str(),"read");
  
    TTreeReader vertex_reader("ShapeAnalysis",infile);
    TTreeReader track_reader("_recoTree",infile);
  
    TTreeReaderValue<int> run(vertex_reader,"run");
    TTreeReaderValue<int> subrun(vertex_reader,"subrun");
    TTreeReaderValue<int> event(vertex_reader,"event");
    TTreeReaderValue<int> vtxid(vertex_reader,"vtxid");
    TTreeReaderValue<double> vtx_x(vertex_reader,"x");
    TTreeReaderValue<double> vtx_y(vertex_reader,"y");
    TTreeReaderValue<double> vtx_z(vertex_reader,"z");
    TTreeReaderValue<std::vector<float>> shfrac(vertex_reader,"shower_frac_Y_v");

      
    TTreeReaderValue<int> runT(track_reader,"run._run");
    TTreeReaderValue<int> subrunT(track_reader,"subrun._subrun");
    TTreeReaderValue<int> eventT(track_reader,"event._event");
    TTreeReaderValue<int> vtxidT(track_reader,"vtx_id");
    TTreeReaderValue<std::vector<double>> trackv(track_reader,"Length_v");
    TTreeReaderValue<std::vector<double>> walldist(track_reader,"closestWall");

    rse.clear();
    rsev.clear();
    cutvars.clear();

    int idx = 0;
    while(vertex_reader.Next()) {

      bool keepMe = true;
    
      int runNum = *run;
      int subNum = *subrun;
      int evtNum = *event;
      int vtxId  = *vtxid;
      float x    = *vtx_x;
      float y    = *vtx_y;
      float z    = *vtx_z;
      std::vector<float>  shfracs = *shfrac;

    
      // Check if maximum shower fraction is too large
      float maxfrac = -1;
      if (shfracs.size() > 0) maxfrac = *max_element(shfracs.begin(), shfracs.end());
      if ( maxfrac<0 ) maxfrac = 0.;

      if (maxfrac > 0.5) keepMe = false;
      // ------------------------------------------- //

      std::vector<float> vtx_cutvars = { (float)runNum, (float)subNum, (float)evtNum, (float)vtxId, x, y, z, maxfrac };
      cutvars.push_back( vtx_cutvars );
      
      // Check if vertex position is outside the fiducial volume
      if (x < 15)         keepMe = false;
      if (x > 256-15)     keepMe = false;
      if (y < -116.5+15)  keepMe = false;
      if (y > 116.5-25)   keepMe = false;
      if (z < 15)         keepMe = false;
      if (z > 1000)       keepMe = false;

      if (z > 700 && z < 740) keepMe = false;
      if (z > 50  && z < 58)  keepMe = false;
      if (z > 90  && z < 98)  keepMe = false;
      if (z > 118 && z < 125) keepMe = false;
      if (z > 244 && z < 250) keepMe = false;
      if (z > 285 && z < 292) keepMe = false;
      if (z > 397 && z < 402) keepMe = false;
      if (z > 415 && z < 421) keepMe = false;
      if (z > 807 && z < 813) keepMe = false;
      if (z > 818 && z < 824) keepMe = false;
      if (z > 872 && z < 880) keepMe = false;

      float m = 1/sqrt(3);

      if (y <  m*z - 180 && y >  m*z - 200) keepMe = false;
      if (y <  m*z - 370 && y >  m*z - 380) keepMe = false;
      if (y <  m*z - 345 && y >  m*z - 350) keepMe = false;
      if (y <  m*z - 549 && y >  m*z - 555) keepMe = false;
      if (y <  m*z - 605 && y >  m*z - 600) keepMe = false;
      if (y <  m*z - 630 && y >  m*z - 625) keepMe = false;
      if (y < -m*z + 435 && y > -m*z + 415) keepMe = false;
      if (y < -m*z + 615 && y > -m*z + 605) keepMe = false;
      if (y < -m*z + 160 && y > -m*z + 155) keepMe = false;
      if (y < -m*z + 235 && y > -m*z + 231) keepMe = false;
      if (y > m*z -117) keepMe = false;
      // ------------------------------------------- //

      if ( verbosity>1 )
        std::cout << "[vtx (" << runNum << "," << subNum << "," << evtNum << "," << vtxId << ")] "
                  << "xyz=(" << x << "," << y << "," << z << ") "
                  << "maxfrac=" << maxfrac
                  << " keepMe=" << keepMe
                  << std::endl;

      
      
      rsev[std::make_tuple(runNum,subNum,evtNum,vtxId)] = keepMe;
    
    }

    if ( verbosity>0 ) {
      std::cout << std::endl;
      std::cout << "[ RSE-VERTEXID ] ================" << std::endl;  
      for ( auto result=rsev.begin(); result!=rsev.end(); result++ ) {
        std::cout << "[" << std::get<0>(result->first) << "," << std::get<1>(result->first) << ","
                  << std::get<2>(result->first) << "," << std::get<3>(result->first) << "]: "
                  << "result=" << result->second << std::endl;
      }
    }

    idx = 0;
    while(track_reader.Next()) {

      bool keepMe = true;
    
      int runNum = *runT;
      int subNum = *subrunT;
      int evtNum = *eventT;
      int vtxId  = *vtxidT;
      std::vector<double> trklen  = *trackv;
      std::vector<double> walld   = *walldist;

      // Check if there are 2 5cm long tracks //
      int n5trks = 0;
      int ntrks  = 0;
      for(float i : trklen) {
        ntrks++;
        if (i >= 5) n5trks++;
      }
      if (ntrks != 2 || n5trks !=2) keepMe = false;
      // ----------------------------------- //

      // Check if tracks are contained  //
      float minwalld = 1.0e9;
      for(float i : walld) {
        if (i < 15) keepMe = false;
        if ( i<minwalld ) minwalld = i;
      }
      if ( (int)cutvars[idx][0]==runNum
           && (int)cutvars[idx][1]==subNum
           && (int)cutvars[idx][2]==evtNum
           && (int)cutvars[idx][3]==vtxId ) {
        cutvars[idx].push_back( ntrks );
        cutvars[idx].push_back( n5trks );
        cutvars[idx].push_back( minwalld );
      }
      else {
        throw std::runtime_error("(RSEV) mismatch in tracker tree");
      }

      if ( verbosity>1 )  {
        std::cout << "[track (" << runNum << "," << subNum << "," << evtNum << "," << vtxId << ")] "
                  << "ntrks=" << ntrks << " n5trks=" << n5trks;
        std::cout << " walld=(";
        for ( float i : walld ) std::cout << i << ",";
        std::cout << ")"
                  << " keep=" << keepMe
                  << std::endl;
      }
      

      if (!keepMe) rsev[std::make_tuple(runNum,subNum,evtNum,vtxId)] = keepMe;
      //std::cout << "track filter: (" << runNum << "," << subNum << "," << evtNum << "," << vtxId << "): " << keepMe << std::endl;

      std::tuple<int,int,int> thisrse = std::make_tuple( runNum, subNum, evtNum );
      auto it = rse.find( thisrse );
      if ( it==rse.end() ) {
        rse[thisrse] = false;
      }
      if (keepMe)
        rse[thisrse] = true;

      idx++;
    }
    infile->Close();

    if ( verbosity>0 ) {
      std::cout << std::endl;
      std::cout << "[ RSE-VERTEXID ] ================" << std::endl;  
      for ( auto result=rsev.begin(); result!=rsev.end(); result++ ) {
        std::cout << "[" << std::get<0>(result->first) << "," << std::get<1>(result->first) << ","
                  << std::get<2>(result->first) << "," << std::get<3>(result->first) << "]: "
                  << "result=" << result->second << std::endl;
      }

      std::cout << "[ RSE ] ================" << std::endl;
      for ( auto result=rse.begin(); result!=rse.end(); result++ ) {
        std::cout << "[" << std::get<0>(result->first) << "," << std::get<1>(result->first) << ","
                  << std::get<2>(result->first) << "]: "
                  << "result=" << result->second << std::endl;
      }
      std::cout << std::endl;
      
      std::cout<< "All Done" <<std::endl;
    }
  
  }

}
}
