#ifndef __larlitecv_PreCuts_h__
#define __larlitecv_PreCuts_h__

//C++ stuff
#include <vector>

// config/storage from LArCV
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

//larlite stuff
#include "Base/DataFormatConstants.h"
#include "DataFormat/ophit.h"


namespace larlitecv {

  class PreCuts {
  public:
    PreCuts() {};
    virtual ~PreCuts() {};

    bool CoordsContained(float x, float y, float z, float edge_x, float edge_y, float edge_z);

    bool EventContained(std::vector<float> nuEndCoords,std::vector<float> leptStartCoords,std::vector<float> lepEndCoords);

    bool IsCCQE(int interactionCode);

    std::vector<float> MakeTimeBin( const larlite::event_ophit& ophitList,int timeBinning,int beamWinStart,int beamWinEnd);
    
    std::vector<float> GetTotalPE(float coincThresh, std::vector<float> flashes);

    float PMTMaxFrac( const larlite::event_ophit& ophitList, std::vector<float> flashBins,float totPE,int timeBinning,float Winstart);

  };

}
#endif

