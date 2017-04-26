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
#pragma link C++ class larlitecv::StopMuStart+;
#pragma link C++ class larlitecv::StopMuSkeleton+;
#pragma link C++ class larlitecv::StopMuAlgo+;
#pragma link C++ class larlitecv::Step3D+;
#pragma link C++ class larlitecv::StopMuTracker+;

#pragma link C++ class larlitecv::FoxStep+;
#pragma link C++ class larlitecv::FoxTrack+;
#pragma link C++ class larlitecv::FoxTrotTrackerAlgoConfig+;
#pragma link C++ class larlitecv::FoxTrotTrackerAlgo+;

#pragma link C++ class larlitecv::StopMuFoxTrotConfig+;
#pragma link C++ class larlitecv::StopMuFoxTrot+;

//ADD_NEW_CLASS ... do not change this line
#endif
