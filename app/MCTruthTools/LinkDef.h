//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace larlitecv;

#pragma link C++ struct larlitecv::CrossingPointAnaData_t+;
#pragma link C++ function larlitecv::analyzeCrossingPoints+;
#pragma link C++ function larlitecv::analyzeCrossingMCTracks+;

#pragma link C++ struct larlitecv::TruthData_t+;
#pragma link C++ function larlitecv::extractTruth+;

//ADD_NEW_CLASS ... do not change this line
#endif
