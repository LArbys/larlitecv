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

#pragma link C++ class larlitecv::FoxStep+;
#pragma link C++ class larlitecv::FoxTrack+;
#pragma link C++ class larlitecv::FoxTrotLeadStraight+;
#pragma link C++ class larlitecv::FoxTrotTrackerAlgoConfig+;
#pragma link C++ class larlitecv::FoxTrotTrackerAlgo+;

#pragma link C++ class larlitecv::PixelQPt+;
#pragma link C++ class std::vector<larlitecv::PixelQPt>+;
#pragma link C++ function larlitecV::getPixelQPts+;
#pragma link C++ function larlitecV::getTrackTotalPixelCharge+;

//ADD_NEW_CLASS ... do not change this line
#endif
























