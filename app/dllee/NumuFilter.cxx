#include "NumuFilter.h"

#include "TTreeReader.h"
#include "TFile.h"

namespace larlitecv {
namespace dllee {

  void NumuFilter::ReturnFilteredDictionary(std::string inputfile,
                                            std::map<std::tuple<int,int,int>,bool>& rse,
                                            std::map<std::tuple<int,int,int,int>,bool>& rsev,
                                            int verbosity ) {
    /*
       for 1m1p, let's do:
       PassPMTPrecut == 1
       PassSimpleCuts == 1
       MaxShrFrac < 0.2
       OpenAng > 0.5
       ChargeNearTrunk > 0
       FailedBoost != 1
       Lepton_EdgeDist > 15
       Proton_EdgeDist > 15
       CosmicBkgScore > 0
       NuBkgScore > 0
    */
    
    TFile* infile  = new TFile(inputfile.c_str(),"read");
  
    TTreeReader vertex_reader("dlana/FinalVertexVariables",infile);
  
    TTreeReaderValue<int> run(vertex_reader,        "run");
    TTreeReaderValue<int> subrun(vertex_reader,     "subrun");
    TTreeReaderValue<int> event(vertex_reader,      "event");
    TTreeReaderValue<int> vtxid(vertex_reader,      "vtxid");
    TTreeReaderValue<int>   precuts(vertex_reader,    "PassPMTPrecut");
    TTreeReaderValue<int>   simplecuts(vertex_reader, "PassSimpleCuts");
    TTreeReaderValue<float> maxshrfrac(vertex_reader, "MaxShrFrac");
    TTreeReaderValue<float> openang(vertex_reader,    "OpenAng");
    TTreeReaderValue<float> qneartrunk(vertex_reader, "ChargeNearTrunk");
    TTreeReaderValue<int>   failboost(vertex_reader,  "FailedBoost");
    TTreeReaderValue<float> lep_edge(vertex_reader,   "Lepton_EdgeDist");
    TTreeReaderValue<float> pro_edge(vertex_reader,   "Proton_EdgeDist");
    TTreeReaderValue<float> cosmicbdt(vertex_reader,  "BDTscore_1mu1p_cosmic");
    TTreeReaderValue<float> nubdt(vertex_reader,      "BDTscore_1mu1p_nu");

    rse.clear();
    rsev.clear();

    int npass = 0;
    while(vertex_reader.Next()) {

      bool keepMe = false;
    
      int runNum = *run;
      int subNum = *subrun;
      int evtNum = *event;
      int vtxId  = *vtxid;

      if ( *precuts==1
           && *simplecuts==1
           && *maxshrfrac<0.2
           && *openang>0.5
           && *qneartrunk>0
           && *failboost!=1
           && *lep_edge>15.0
           && *pro_edge>15.0
           && *cosmicbdt>0.0
           && *nubdt>0.0 ) {
        keepMe = true;
        npass++;
      }

      rse[std::make_tuple(runNum,subNum,evtNum)] = keepMe;
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
      std::cout << "[Number passed: " << npass << "]" << std::endl;
    }

  }

}
}
