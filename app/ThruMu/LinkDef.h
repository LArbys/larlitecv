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
#pragma link C++ class larlitecv::BoundaryMatchArrays+;
#pragma link C++ class larlitecv::BoundaryEndPt+;
#pragma link C++ class larlitecv::BoundarySpacePoint+;
#pragma link C++ class std::vector< larlitecv::BoundaryEndPt >+;
#pragma link C++ class std::vector< larlitecv::BoundarySpacePoint >+;
#pragma link C++ class larlitecv::BMTrackCluster2D+;
#pragma link C++ class larlitecv::BMTrackCluster3D+;

#pragma link C++ class larlitecv::BoundaryMuonTaggerAlgoConfig+;
#pragma link C++ class larlitecv::BoundaryMuonTaggerAlgo+;
#pragma link C++ class larlitecv::BoundaryMuonTagger+;
#pragma link C++ class larlitecv::FlashMuonTaggerConfig+;
#pragma link C++ class larlitecv::FlashMuonTaggerConfigAlgo+;

#pragma link C++ class larlitecv::PointInfo+;
#pragma link C++ class larlitecv::PointInfoList+;
#pragma link C++ class larlitecv::Linear3DChargeTaggerConfig+;
#pragma link C++ class larlitecv::Linear3DChargeTagger+;

#pragma link C++ class larlitecv::AStar3DNode+;
#pragma link C++ class std::vector< larlitecv::AStar3DNode >+;
#pragma link C++ class larlitecv::AStar3DAlgoConfig+;
#pragma link C++ class larlitecv::AStar3DAlgo+;

//#pragma link C++ class larlitecv::BezierCurve+;


//ADD_NEW_CLASS ... do not change this line
#endif
