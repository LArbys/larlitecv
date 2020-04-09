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
#pragma link C++ function larlitecv::make_evdisp;
#pragma link C++ function larlitecv::make_evdisp_single;
#pragma link C++ function larlitecv::make_evdisp_triple;
#pragma link C++ function larlitecv::is_close_enough;
#pragma link C++ function larlitecv::distance_between_pt;
#pragma link C++ function larlitecv::distance_between_pt2d;
#pragma link C++ function larlitecv::getProjectedPixel;
#pragma link C++ function larlitecv::print_signal;
#pragma link C++ class larlitecv::DQDXBuilder;

//ADD_NEW_CLASS ... do not change this line
#endif
