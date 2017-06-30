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
#pragma link C++ enum larlitecv::BoundaryEnd_t+;
#pragma link C++ class larlitecv::BoundaryEndPt+;
#pragma link C++ class larlitecv::BoundarySpacePoint+;
#pragma link C++ class std::vector< larlitecv::BoundaryEndPt >+;
#pragma link C++ class std::vector< larlitecv::BoundarySpacePoint >+;
#pragma link C++ class larlitecv::BMTrackCluster2D+;
#pragma link C++ class larlitecv::BMTrackCluster3D+;

#pragma link C++ function larlitecv::getTrackPixelsFromImages+;
#pragma link C++ function larlitecv::getTrackPixelsFromImagesNoBadCh+;
#pragma link C++ function larlitecv::BoundaryEndNames+;

//ADD_NEW_CLASS ... do not change this line
#endif
