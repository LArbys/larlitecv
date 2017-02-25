
//PreCu stuff
#include "PreCuts.h"

// C++
//#include <algorithm>
#include <utility>

namespace larlitecv {

  bool PreCuts::CoordsContained(float x, float y, float z, float edge_x, float edge_y, float edge_z) {

    float xmax =  256.25;
    float xmin =  0.;
    float ymax =  116.5;
    float ymin = -116.5;
    float zmax =  1036.8;
    float zmin =  0;

    if (x<(xmax-edge_x) && x>(xmin+edge_x) && y<(ymax-edge_y) && y>(ymin+edge_y) && z<(zmax-edge_z) && z>(zmin+edge_z)){
      return true;
    } 
    else {
	return false;
    }

  }


  bool PreCuts::EventContained(std::vector<float> nuEndCoords,std::vector<float> lepStartCoords,std::vector<float> lepEndCoords) {

    float fidCut = 5.0;

    bool NuContained = CoordsContained(nuEndCoords[0],nuEndCoords[1],nuEndCoords[2],fidCut,fidCut,fidCut);
    bool LepStartContained = CoordsContained(lepStartCoords[0],lepStartCoords[1],lepStartCoords[2],fidCut,fidCut,fidCut);
    bool LepEndContained = CoordsContained(lepEndCoords[0],lepEndCoords[1],lepEndCoords[2],fidCut,fidCut,fidCut);

    if (NuContained && LepStartContained && LepEndContained) {
      return true;
    }
    else {
      return false;
    }

  }

  bool PreCuts::IsCCQE(int interactionCode) {
    if (interactionCode == 1001) {
      return true;
    }
    else {
      return false;
    }
  }


  std::vector<float> PreCuts::MakeTimeBin(const larlite::event_ophit& ophitList, int timeBinning,int beamWinStart,int beamWinEnd) {

    int numBins    = (beamWinEnd - beamWinStart)/timeBinning + 1;
    float tickSize = 0.015625;

    std::vector<float> PEbyTimeBin(numBins,0.);

    for (unsigned int i=0; i<ophitList.size(); i++) {
      
      if (ophitList[i].PeakTime()/tickSize > beamWinEnd || ophitList[i].PeakTime()/tickSize < beamWinStart) {
	continue;
      }
      
      PEbyTimeBin[(int)((ophitList[i].PeakTime()/tickSize - beamWinStart)/timeBinning)]+=ophitList[i].PE();

    }

    return PEbyTimeBin;

  }


  std::vector<float> PreCuts::GetTotalPE(float coincThresh, std::vector<float> flashes) {

    float totPE = 0;
    bool flashFound = false;
    std::vector<float> PEinfo;
    std::vector<float> flashBins;

    for (unsigned int i=0; i<flashes.size(); i++) {

      if (flashes[i] > coincThresh && flashFound == true) {
	flashBins.push_back(i);
	totPE+=flashes[i];
      }

      if (flashes[i] < coincThresh && flashFound == true) {
	break;
      }

      if (flashes[i] > coincThresh && flashFound == false) {
	totPE+=flashes[i];
	flashBins.push_back(i);
	flashFound = true;
      }

    }

    PEinfo.push_back(totPE);
    for (unsigned int i=0; i < flashBins.size(); i++) {
      PEinfo.push_back(flashBins[i]);
    }

    return PEinfo;

  }
    
  float PreCuts::PMTMaxFrac(const larlite::event_ophit& ophitList,std::vector<float> flashBins,float totPE,int timeBinning, float Winstart) {

    std::vector<float> PEFracbyPMT(32,0.);
    float maxFrac;

    for (unsigned int i=0; i<ophitList.size(); i++) {

      int ophitBin = (ophitList[i].PeakTime() - Winstart)/timeBinning;

      if (std::find(flashBins.begin(),flashBins.end(),ophitBin)!=flashBins.end()) {
	PEFracbyPMT[ophitList[i].OpChannel()]+=ophitList[i].PE()/totPE;
      }

    }
    
    maxFrac = *max_element(PEFracbyPMT.begin(),PEFracbyPMT.end());

    return maxFrac;

  }
    
    


}
