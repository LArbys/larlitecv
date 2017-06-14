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
//ADD_NEW_CLASS ... do not change this line
#pragma link C++ class larlitecv::ClusterGroupAlgoConfig+;
#pragma link C++ class larlitecv::ClusterGroupAlgo+;
#pragma link C++ class larlitecv::ClusterGroupMatchingAlgo+;
#endif
