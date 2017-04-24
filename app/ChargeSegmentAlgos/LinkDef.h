//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ struct larlitecv::Segment3D_t+;
#pragma link C++ class std::vector<larlitecv::Segment3D_t>+;
#pragma link C++ struct larlitecv::Segment2D_t+;
#pragma link C++ class std::vector<larlitecv::Segment2D_t>+;
#pragma link C++ class larlitecv::Segment3DAlgo+;

#pragma link C++ namespace larlitecv;
#pragma link C++ class larlitecv::RadialSegmentSearch+;
#pragma link C++ class larlitecv::RadialHit_t+;
#pragma link C++ class std::vector<larlitecv::RadialHit_t>+;

//ADD_NEW_CLASS ... do not change this line
#endif
























