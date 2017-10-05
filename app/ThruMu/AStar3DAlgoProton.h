#ifndef __ASTAR_3D_ALGOPROTON__
#define __ASTAR_3D_ALGOPROTON__

/**

 AStar algorithm assuming 3D!! grid points.

 Uses Image2D to hold image.

 **/

#include <iostream>
#include <queue>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <array>

// larcv
#include "DataFormat/Image2D.h"
#include "Base/PSet.h"
#include "Reco3D/AStar3DAlgo.h"

namespace larlitecv {

    // ALGO
    class AStar3DAlgoProton {

        AStar3DAlgoProton() { verbose=2; doPxValEstimate = false; };
    public:

        AStar3DAlgoProton( larcv::AStar3DAlgoConfig config ) {
            _config = config;
            setVerbose( _config.verbosity );
        };
        virtual ~AStar3DAlgoProton() {};

        void setVerbose( int v ) { verbose = v; };
        void setPixelValueEsitmation( bool doIt ) { doPxValEstimate = doIt; };

        std::vector<larcv::AStar3DNode> findpath( const std::vector<larcv::Image2D>& img_v,
						  const std::vector<larcv::Image2D>& badch_v,
						  const std::vector<larcv::Image2D>& tagged_v,
						  const int start_row, const int goal_row,
						  const std::vector<int>& start_cols,
						  const std::vector<int>& goal_cols,
						  int& goal_reached );

        std::vector<larcv::AStar3DNode> makeRecoPath( larcv::AStar3DNode* start, larcv::AStar3DNode* goal, bool& path_completed );

        void evaluateNeighborNodes( larcv::AStar3DNode* current,
				    const larcv::AStar3DNode* start,
				    const larcv::AStar3DNode* goal,
				    const std::vector<larcv::Image2D>& img_v,
				    const std::vector<larcv::Image2D>& badch_v,
				    const std::vector<larcv::Image2D>& tagged_v,
				    larcv::AStar3DNodePtrList& openset,
				    larcv::AStar3DNodePtrList& closedset,
				    larcv::Lattice& lattice );

        bool evaluteLatticePoint( const larcv::A3DPixPos_t& latticept,
				  larcv::AStar3DNode* current,
				  const larcv::AStar3DNode* start,
				  const larcv::AStar3DNode* goal,
				  const std::vector<larcv::Image2D>& img_v,
				  const std::vector<larcv::Image2D>& badch_v,
				  const std::vector<larcv::Image2D>& tagged_v,
				  larcv::AStar3DNodePtrList& openset,
				  larcv::AStar3DNodePtrList& closedset,
				  larcv::Lattice& lattice );

        float distanceFromCentralLine( const std::vector<float>& start_tyz, const std::vector<float>& end_tyz, const std::vector<float>& testpt_tyz );
	
        const std::vector<larcv::Image2D>& getScoreImages() { return m_visualizedimgs; }
	
    protected:
	
	larcv::AStar3DAlgoConfig _config;
        int verbose;
        bool doPxValEstimate;
        
        std::vector< larcv::Image2D > m_visualizedimgs;
    };
    
    
}

#endif
