#include "BoundaryMatchAlgo.h"

#include <ctime>
#include <iostream>
#include <assert.h>

#include "BoundaryMatchArrays.h"
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  BoundaryMatchList::~BoundaryMatchList() {
    delete [] uvec;
    delete [] vvec;
    delete [] yvec;
  }

  void BoundaryMatchList::load( BoundaryMatchArrays::Boundary_t t ) {
    const clock_t begin_time = clock();
    
    BoundaryMatchArrays match_arrays;
    
    m_combos.resize( match_arrays.nmatches(t) );

    const larcv::WireData& plane0data = larcv::UBWireTool::getWireData(0);
    const larcv::WireData& plane1data = larcv::UBWireTool::getWireData(1);
    const larcv::WireData& plane2data = larcv::UBWireTool::getWireData(2);

    for (int imatch=0; imatch<match_arrays.nmatches(t); imatch++) {
      int u, v, y;
      match_arrays.getMatch( t, imatch, u, v, y );

      // define and store the combo
      float pos[3] = { 0., 0., 0.};
      switch ( t ) {
      case BoundaryMatchArrays::kTop:
	pos[0] = top_positions[imatch][0];
	pos[1] = top_positions[imatch][1];
	pos[2] = top_positions[imatch][2];
	break;
      case BoundaryMatchArrays::kBottom:
	pos[0] = bottom_positions[imatch][0];
	pos[1] = bottom_positions[imatch][1];
	pos[2] = bottom_positions[imatch][2];
	break;
      case BoundaryMatchArrays::kUpstream:
	pos[0] = upstream_positions[imatch][0];
	pos[1] = upstream_positions[imatch][1];
	pos[2] = upstream_positions[imatch][2];
	break;
      case BoundaryMatchArrays::kDownstream:
	pos[0] = downstream_positions[imatch][0];
	pos[1] = downstream_positions[imatch][1];
	pos[2] = downstream_positions[imatch][2];
	break;
      default:
	assert(false);
	break;
      }
      BoundaryCombo combo( u, v, y, pos[0], pos[1], pos[2] );
      m_combos[imatch] = combo;

      // build index set trees
      auto it_u = u_indices.find((idx)u);
      auto it_v = v_indices.find((idx)v);
      auto it_y = y_indices.find((idx)y);
      if ( it_u==u_indices.end() ) {
	u_indices.insert( std::pair< idx, std::set<idx> >(u, std::set<idx>() ) );
      }
      if ( it_v==v_indices.end() ) {
	v_indices.insert( std::pair< idx, std::set<idx> >(v, std::set<idx>() ) );
      }
      if ( it_y==y_indices.end() ) {
	y_indices.insert( std::pair< idx, std::set<idx> >(y, std::set<idx>() ) );
      }
      u_indices[(idx)u].insert((idx)imatch);
      v_indices[(idx)v].insert((idx)imatch);
      y_indices[(idx)y].insert((idx)imatch);
      if ( u>uid_max ) uid_max = u;
      if ( v>vid_max ) vid_max = v;
      if ( y>yid_max ) yid_max = y;
    }
    uvec = new int[m_combos.size()];
    vvec = new int[m_combos.size()];
    yvec = new int[m_combos.size()];

    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    std::cout << "boundary match list loaded " << elapsed_secs << " secs" << std::endl;
    
  }
  
  std::vector<BoundaryCombo> BoundaryMatchList::findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires, 
							    int urange, int vrange, int yrange, const std::vector<larcv::Image2D>& badchs, bool use_badchs ) {
    const clock_t begin_time = clock();

    memset(uvec,0,sizeof(int)*m_combos.size());
    memset(vvec,0,sizeof(int)*m_combos.size());
    memset(yvec,0,sizeof(int)*m_combos.size());

    // get dsfactor
    int dsfactor = (int)badchs[0].meta().pixel_width();
    
    // filter u wires

    // loop over u-wires
    //std::cout << " u: ";
    for (int uwid=0; uwid<uwires.size(); uwid++) {
      if ( uwires[uwid]==0 ) continue;	
      idx uidx = (idx)( uwid );
      //std::cout << " " << uidx;
      // copy indices into passu
      for ( auto &iiu : u_indices[ uidx ] ) {
	//passu.insert( iiu );
	uvec[iiu] = 1;
      }
    }
    //std::cout << std::endl;
    
    // loop over v-wires, filter matched indices
    //std::cout << " v: ";
    for (int vwid=0; vwid<vwires.size(); vwid++) {
      if ( vwires[vwid]== 0 ) continue;
      idx vidx = (idx)( vwid );
      //std::cout << " " << vidx;
      // copy indices into passv if in passu
      for ( auto &iiv : v_indices[ vidx ] ) {
	vvec[iiv] = 1;
      }
    }
    //std::cout << std::endl;

    // loop over y-wires, filter matched indices
    //std::cout << " y: ";
    for (int ywid=0; ywid<ywires.size(); ywid++) {
      if ( ywires[ywid]== 0 ) continue;
      idx yidx = (idx)(ywid);
      //std::cout << " " << yidx;
      for ( auto &iiy : y_indices[ yidx ] ) {
	yvec[iiy] = 1;
      }
    }
    //std::cout << std::endl;

    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    //std::cout << "boundary match searched in " << elapsed_secs << " secs" << std::endl;    

    std::vector< BoundaryCombo > combo_v;
    combo_v.reserve( m_combos.size() );
    int ncombos=0;
    for (int idx=0; idx<m_combos.size(); idx++) {
      int votes = uvec[idx]+vvec[idx]+yvec[idx];
      if ( votes>=3 ) {
	combo_v.push_back( m_combos.at(idx) );
	ncombos++;
      }
      else if ( votes==2 && use_badchs ) {
	int badchp = -1;
	int col = -1;
	if ( uvec[idx]==0 ) { 
	  badchp = 0;
	  col = uvec[idx]/dsfactor;
	}
	else if ( vvec[idx]==0 ) {
	  badchp = 1;
	  col = vvec[idx]/dsfactor;
	}
	else if ( yvec[idx]==0 ) {
	  badchp = 2;
	  col = yvec[idx]/dsfactor;
	}
	if ( badchs[badchp].pixel( 0, col ) > 0 ) {
	  combo_v.push_back( m_combos.at(idx) );
	  ncombos++;
	}
      }
    }
    
    elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    //std::cout << "ncombos=" << ncombos << ". boundary match searched/stored in " << elapsed_secs << " secs" << std::endl;    
    
    return combo_v;
  }

  BoundaryMatchAlgo::BoundaryMatchAlgo() {
    matchlist[0] = new BoundaryMatchList( BoundaryMatchArrays::kTop );
    matchlist[1] = new BoundaryMatchList( BoundaryMatchArrays::kBottom );
    matchlist[2] = new BoundaryMatchList( BoundaryMatchArrays::kUpstream );
    matchlist[3] = new BoundaryMatchList( BoundaryMatchArrays::kDownstream );
  }

  BoundaryMatchAlgo::~BoundaryMatchAlgo() {
    for (int i=0; i<4; i++) {
      delete matchlist[i];
      matchlist[i] = NULL;
    }
  }

  void BoundaryMatchAlgo::findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires,
				      int urange, int vrange, int yrange,
				      const std::vector<larcv::Image2D>& badchs, bool use_badchs,
				      std::vector<  std::vector<BoundaryCombo>  >& boundary_combos ) {
    for (int i=0; i<4; i++) {
      std::vector<BoundaryCombo> combos = matchlist[i]->findCombos( uwires, vwires, ywires, urange, vrange, yrange, badchs, use_badchs );
      boundary_combos.emplace_back( std::move(combos) );
    }
  }
  
}
