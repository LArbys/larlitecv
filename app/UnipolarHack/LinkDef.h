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
#pragma link C++ class larlitecv::UnipolarROI_t+;
#pragma link C++ class std::vector< larlitecv::UnipolarROI_t >+;
#pragma link C++ class larlitecv::UnipolarHackAlgo+;

//ADD_NEW_CLASS ... do not change this line
#endif