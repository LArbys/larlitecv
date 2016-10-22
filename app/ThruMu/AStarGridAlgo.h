#ifndef __ASTAR_GRID_ALGO__
#define __ASTAR_GRID_ALGO__

/** 

AStar algorithm assuming 2D grid points.

Uses Image2D to hold image.

 **/

#include <iostream>
#include <queue>
#include <set>
#include <map>

#include "DataFormat/Image2D.h"

namespace larlitecv {

  class AStarNode  {
  public:
    AStarNode() { col = 0; row=0; fscore=gscore=0.0; id=previd=-1; };
    AStarNode( int col_, int row_ ) { 
      col = col_; 
      row = row_; 
      fscore=0.0;
      gscore=0.0;
      id = -1;
      previd = -1;
    };
    AStarNode( const AStarNode& src ) {
      col = src.col;
      row = src.row;
      fscore = src.fscore;
      gscore = src.gscore;
      id = src.id;
      previd = src.previd;
    };
    virtual ~AStarNode() {};

    int col;
    int row;
    float fscore;
    float gscore;
    int id;
    int previd;

    bool operator <(const AStarNode& rhs ) const {
      //std::cout << "astar< " << fscore << " < " << rhs.fscore << std::endl;
      if ( fscore > rhs.fscore ) // lower scores get priority
	return true;
      return false;
    };
    bool operator==(const AStarNode& rhs ) const {
      if ( col==rhs.col && row==rhs.row ) return true;
      return false;
    };

    bool operator!=(const AStarNode& rhs ) const {
      if ( col!=rhs.col || row!=rhs.row ) return true;
      return false;
    };
    
  };

  struct node_pos_compare {
    bool operator() (const AStarNode& lhs, const AStarNode& rhs ) const {
      bool result = false;
      if ( lhs.row<rhs.row ) result=true;
      else if ( lhs.row==rhs.row && lhs.col<rhs.col ) result=true;
      //std::cout << "(" << lhs.row << "," << lhs.col << ") < (" << rhs.row << "," << rhs.col << ") result=" << result << std::endl;
      return result;
    };
  };

  class AStarSet : protected std::priority_queue< larlitecv::AStarNode > {
  public:
    AStarSet() {};
    virtual ~AStarSet() {};

    void add( AStarNode& anode ) { push(anode); nodes.insert(anode); };
    AStarNode gettop() { 
      AStarNode t=top(); 
      pop();
      auto it=nodes.find(t);
      if (it!=nodes.end() ) {
	nodes.erase( it ); 
      }
      return t; 
    };
    bool contains( const AStarNode& anode ) { if ( nodes.find(anode)!=nodes.end() ) return true; return false; };
    float get_gscore( const AStarNode& node ) { auto it=nodes.find(node); if ( it!=nodes.end() ) return (*it).gscore; return -1.0; };
    float get_fscore( const AStarNode& node ) { auto it=nodes.find(node); if ( it!=nodes.end() ) return (*it).fscore; return -1.0; };
    size_t nentries() { return size(); };

#ifndef __CINT__
#ifndef __CLING__
    void addemplace( AStarNode&& anode ) { nodes.insert(anode); emplace(anode); };
#endif
#endif

    std::set<AStarNode, node_pos_compare> nodes;
    

  };

  class AStarAlgoConfig {
  public:
    AStarAlgoConfig() {};
    virtual ~AStarAlgoConfig() {};

    std::vector<float> astar_threshold;
    std::vector<int>   astar_neighborhood;
  };

  class AStarGridAlgo {

    AStarGridAlgo() { verbose=2; };
  public:

    AStarGridAlgo( AStarAlgoConfig config ) { _config = config; m_badchimg = NULL; };
    virtual ~AStarGridAlgo() {};
    
    void setVerbose( int v ) { verbose = v; };
    void setBadChImage( const larcv::Image2D& badchimg ) { m_badchimg = &badchimg; };
    std::vector<AStarNode> findpath( const larcv::Image2D& img, int start_row, int start_col, int goal_row, int goal_col, float thresh, bool use_bad_chs=false );
    std::vector<AStarNode> makeRecoPath( const AStarNode& start, const AStarNode& goal, const std::map< AStarNode,AStarNode, node_pos_compare >& camefrom, 
					 int origin_row, int origin_col, bool& path_completed );

  protected:

    AStarAlgoConfig _config;
    int verbose;

    const larcv::Image2D* m_badchimg;
  };


}

#endif
