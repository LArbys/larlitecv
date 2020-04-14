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
#pragma link C++ namespace larlitecv::ssnetshowerreco;
#pragma link C++ namespace std;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector< std::vector<float> >+;
#pragma link C++ class larlitecv::ssnetshowerreco::Utils+;
#pragma link C++ class larlitecv::ssnetshowerreco::SecondShower+;
#pragma link C++ class larlitecv::ssnetshowerreco::SSNetShowerReco+;
#pragma link C++ class larlitecv::ssnetshowerreco::ShowerRecoUtil+;

//ADD_NEW_CLASS ... do not change this line
#endif
