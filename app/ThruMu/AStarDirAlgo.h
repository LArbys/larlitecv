#ifndef __ASTAR_DIR_ALGO__
#define __ASTAR_DIR_ALGO__

/** 

AStar algorithm assuming 2D grid points. Also, exploration favors nodes that are in 
the current direction the path is on.

Uses Image2D to hold image.

 **/

#include <iostream>
#include <queue>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "DataFormat/Image2D.h"

namespace larlitecv {

  class AStarDirNode  {
  public:
    AStarDirNode() {
      col = 0;
      row=0;
      fscore=gscore=jump_cost=turn_cost=0.0;
      turn_cost_discrete = 0;
      id=previd=-1;
      closed = false;
      dir2d.resize(2,0.0);
      prev = NULL;
    };
    AStarDirNode( int col_, int row_, std::vector<float> dir2d_ ) { 
      col = col_; 
      row = row_; 
      fscore=0.0;
      gscore=0.0;
      jump_cost=0.0;
      turn_cost=0.0;
      turn_cost_discrete=0;
      id = -1;
      previd = -1;
      dir2d = dir2d_;
      closed = false;
      prev = NULL;
    };
    AStarDirNode( const AStarDirNode& src ) {
      col = src.col;
      row = src.row;
      fscore = src.fscore;
      gscore = src.gscore;
      jump_cost=src.jump_cost;
      turn_cost=src.turn_cost;
      turn_cost_discrete=src.turn_cost_discrete;
      id = src.id;
      previd = src.previd;
      dir2d = src.dir2d;
      closed = src.closed;
      prev = src.prev;
    };
    virtual ~AStarDirNode() {};

    int col; // position on image
    int row; // position on image
    float fscore; // score of path: past path + to goal heuristic
    float gscore; // past path score
    float jump_cost; // cost to jump to next node: this is control order with which we evaluate successor nodes
    float turn_cost; // costs of zig-zagging
    int   turn_cost_discrete; // number of turns (defined as node changes of direction above some threshold)

    std::vector<float> dir2d;
    int id;
    int previd;
    bool closed;
    AStarDirNode* prev;
    //AStarDirNode* next;

    bool operator <(const AStarDirNode& rhs ) const {
      //std::cout << "astar< " << fscore << " < " << rhs.fscore << std::endl;
      if ( fscore > rhs.fscore ) // lower scores get priority
	return true;
      return false;
    };
    bool operator==(const AStarDirNode& rhs ) const {
      if ( col==rhs.col && row==rhs.row ) return true;
      return false;
    };

    bool operator!=(const AStarDirNode& rhs ) const {
      if ( col!=rhs.col || row!=rhs.row ) return true;
      return false;
    };
    
  };

  struct dirnode_pos_compare {
    bool operator() (const AStarDirNode& lhs, const AStarDirNode& rhs ) const {
      bool result = false;
      if ( lhs.row<rhs.row ) result=true;
      else if ( lhs.row==rhs.row && lhs.col<rhs.col ) result=true;

      return result;
    };
  };

  class AStarDirNodePtrList : public std::vector< AStarDirNode* > {

    public:
      AStarDirNodePtrList() {};
      virtual ~AStarDirNodePtrList() {};

    protected:

    // we want to sort from largest to smallest, because we will pop the nodes off the back of the various sets
    struct dirnodeptr_compare_t {
      bool operator() (const AStarDirNode* lhs, const AStarDirNode* rhs ) const {
        if ( lhs->fscore > rhs->fscore )
          return true;
        return false;
      };
    } m_thecomparator;

    public:
    void sort() {
      std::sort(this->begin(),this->end(), m_thecomparator);
    }

    void addnode( AStarDirNode* elem ) {
      push_back( elem );
      sort();
    }
  };


  class AStarDirSet : protected std::priority_queue< larlitecv::AStarDirNode > {
  public:
    AStarDirSet() {};
    virtual ~AStarDirSet() {};

    void add( AStarDirNode& anode ) { push(anode); nodes.insert(anode); };
    AStarDirNode gettop() { 
      AStarDirNode t=top(); 
      pop();
      auto it=nodes.find(t);
      if (it!=nodes.end() ) {
	nodes.erase( it ); 
      }
      return t; 
    };
    bool contains( const AStarDirNode& anode ) { if ( nodes.find(anode)!=nodes.end() ) return true; return false; };
    float get_gscore( const AStarDirNode& node ) { auto it=nodes.find(node); if ( it!=nodes.end() ) return (*it).gscore; return -1.0; };
    float get_fscore( const AStarDirNode& node ) { auto it=nodes.find(node); if ( it!=nodes.end() ) return (*it).fscore; return -1.0; };
    size_t nentries() { return size(); };

#ifndef __CINT__
#ifndef __CLING__
    void addemplace( AStarDirNode&& anode ) { nodes.insert(anode); emplace(anode); };
#endif
#endif

    std::set<AStarDirNode, dirnode_pos_compare> nodes;
    

  };

  class AStarDirAlgoConfig {
  public:
    AStarDirAlgoConfig() {
      astar_start_padding = 0;
      astar_end_padding = 0;
    };
    virtual ~AStarDirAlgoConfig() {};

    std::vector<float> astar_threshold;
    std::vector<int>   astar_neighborhood;
    int astar_start_padding;
    int astar_end_padding;
    int image_padding;
  };

  class AStarDirAlgo {

    AStarDirAlgo() { verbose=2; };
  public:

    typedef std::pair<int,int> PixPos_t;
    typedef std::pair< PixPos_t, AStarDirNode* > pos_key_t;

    AStarDirAlgo( AStarDirAlgoConfig config ) { _config = config; m_badchimg = NULL; };
    virtual ~AStarDirAlgo() {};
    
    void setVerbose( int v ) { verbose = v; };
    void setBadChImage( const larcv::Image2D& badchimg ) { m_badchimg = &badchimg; };
    std::vector<AStarDirNode> findpath( const larcv::Image2D& img, int start_row, int start_col, int goal_row, int goal_col, float thresh, bool use_bad_chs=false );
    std::vector<AStarDirNode> makeRecoPath( AStarDirNode* start, AStarDirNode* goal, int origin_row, int origin_col, bool& path_completed );

    larcv::Image2D visualizeScores( std::string score_name, const larcv::Image2D& orig_img, 
      const int min_c, const int min_r, const int win_c, const int win_r, 
      const std::map<PixPos_t,AStarDirNode*>  position_lookup);

    const std::vector<larcv::Image2D>& getScoreImages() { return m_visualizedimgs; }

  protected:

    AStarDirAlgoConfig _config;
    int verbose;

    void evaluateNeighborNodes( AStarDirNode* current, const AStarDirNode* start, const AStarDirNode* goal,
      AStarDirNodePtrList& openset, AStarDirNodePtrList& closedset, 
      const int neighborhood_size, const int min_c, const int min_r, const int win_c, const int win_r, 
      const larcv::Image2D& img, const larcv::ImageMeta& meta, const bool use_bad_chs, std::map< PixPos_t, AStarDirNode* >& position_lookup );

   void evaluateBadChNeighbors( AStarDirNode* current, const AStarDirNode* start, const AStarDirNode* goal,
      AStarDirNodePtrList& openset, AStarDirNodePtrList& closedset, 
      const int neighborhood_size, const int min_c, const int min_r, const int win_c, const int win_r, 
      const larcv::Image2D& img, const larcv::ImageMeta& meta, const bool use_bad_chs, std::map< PixPos_t, AStarDirNode* >& position_lookup);


    const larcv::Image2D* m_badchimg;
    std::vector< larcv::Image2D > m_visualizedimgs;
  };


}

#endif
