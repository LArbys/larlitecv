//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Classes
#pragma link C++ class std::map<string, std::vector<string> >+;
#pragma link C++ namespace larlitecv;
#pragma link C++ class larlitecv::LarliteFileManager+;
#pragma link C++ class larlitecv::LarcvFileManager+;
#pragma link C++ class larlitecv::DataCoordinator+;
//ADD_NEW_CLASS ... do not change this line

#endif







