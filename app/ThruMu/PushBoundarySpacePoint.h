#ifndef __PushBoundaryPoints_h__

#include <vector>
#include "DataFormat/Image2D.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgoTypes.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgoConfig.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgo.h"

namespace larlitecv {
    class PushBoundarySpacePoint {
    public:
        PushBoundarySpacePoint();
        virtual ~PushBoundarySpacePoint() {};

    BoundarySpacePoint pushPoint( const larlitecv::BoundarySpacePoint& boundarypoint, const std::vector<larcv::Image2D>& img_v,
                                    const std::vector<larcv::Image2D>& badch_v );


    void clear();

    protected:
        // submethods
        larlitecv::FoxTrack runFoxTrot( const larlitecv::BoundarySpacePoint& sp, const std::vector<larcv::Image2D>& img_v,
                    const std::vector<larcv::Image2D>& badch_v ); //< runs foxtrottrackeralgo to try and go from spacepoint to end of track

        larlitecv::BoundarySpacePoint scanTrackForEndPoint( const BoundarySpacePoint& original, const FoxTrack& track,
                    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ); //< scans along foxtrot track to look for new end point

        larlitecv::BoundarySpacePoint evalEndPoint( const larlitecv::BoundarySpacePoint& sp ); //< evaluate end point quality -- and reclassify based on location

        // algo-variables
        std::vector<larlitecv::FoxTrack>           m_tracklist_v;
        std::vector<larlitecv::BoundarySpacePoint> m_endpoints_v;

        // used algos
        larlitecv::FoxTrotTrackerAlgoConfig m_foxalgo_cfg;
        larlitecv::FoxTrotTrackerAlgo       m_foxalgo;
        larlitecv::Segment3DAlgo            m_segment_algo;

    };
}

#endif