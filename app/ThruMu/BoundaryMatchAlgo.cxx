#include "BoundaryMatchAlgo.h"

#include <ctime>
#include <iostream>

#include "BoundaryMatchArrays.h"

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

    for (int imatch=0; imatch<match_arrays.nmatches(t); imatch++) {
      int u, v, y;
      match_arrays.getMatch( t, imatch, u, v, y );
      BoundaryCombo combo( u, v, y );
      m_combos[imatch] = combo;
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
							     int urange, int vrange, int yrange ) {
    const clock_t begin_time = clock();

    memset(uvec,0,sizeof(int)*m_combos.size());
    memset(vvec,0,sizeof(int)*m_combos.size());
    memset(yvec,0,sizeof(int)*m_combos.size());
    
    // we keep filtering indices
    //std::set<idx> passu;
    //std::set<idx> passv;
    //std::set<idx> passy;
    
    // filter u wires

    // loop over u-wires
    //std::cout << "uhits: " << uwires.size() << std::endl;
    for (auto &uwid : uwires ) {

      // loop over neighborhood
      for (int iu=-urange; iu<=urange; iu++) {
	if ( uwid+iu<0 ) continue;
	if ( uwid+iu>=uid_max ) break;
	
	idx uidx = (idx)( uwid+iu );
	// copy indices into passu
	for ( auto &iiu : u_indices[ uidx ] ) {
	  //passu.insert( iiu );
	  uvec[iiu] = 1;
	}
      }
    }
    //std::cout << " passu=" << passu.size() << std::endl;

    // loop over v-wires, filter matched indices
    //std::cout << "vhits: " << vwires.size() << std::endl;
    for (auto &vwid : vwires ) {
      // loop over neighborhood of wires
      for (int iv=-vrange; iv<=vrange; iv++) {
        if ( vwid+iv<0 ) continue;
	if ( vwid+iv>=vid_max ) break;

	idx vidx = (idx)( vwid+iv );
	// copy indices into passv if in passu
	for ( auto &iiv : v_indices[ vidx ] ) {
	  //if ( passu.find( iiv )!=passu.end() ) {
	  //passv.insert( iiv );
	  //}
	  vvec[iiv] = 1;
	}
      }
    }
    //std::cout << " passv=" << passv.size() << std::endl;

    // loop over y-wires, filter matched indices
    //std::cout << "yhits: " << ywires.size() << std::endl;
    for (auto &ywid : ywires ) {
      // loop over neighborhood of wires
      for (int iy=-yrange; iy<=yrange; iy++) {
        if ( ywid+iy<0 ) continue;
	if ( ywid+iy>=yid_max ) break;

	idx yidx = (idx)( ywid+iy );
	//if ( y_indices.find( yidx )!=y_indices.end() ) {
	// copy indices into passy if in passu
	for ( auto &iiy : y_indices[ yidx ] ) {
	  //if ( passv.find( iiy )!=passv.end() ) {
	  //passy.insert( iiy );
	  //}
	  yvec[iiy] = 1;
	}
      }
    }
    //std::cout << " passy=" << passy.size() << std::endl;

    float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    //std::cout << "boundary match searched in " << elapsed_secs << " secs" << std::endl;    

    std::vector< BoundaryCombo > combo_v;
    combo_v.reserve( m_combos.size() );
    int ncombos=0;
    for (int idx=0; idx<m_combos.size(); idx++) {
      if ( uvec[idx] & vvec[idx] & yvec[idx]==1 ) {
	combo_v.push_back( m_combos.at(idx) );
	ncombos++;
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
				      std::vector<  std::vector<BoundaryCombo>  >& boundary_combos ) {
    for (int i=0; i<4; i++) {
      std::vector<BoundaryCombo> combos = matchlist[i]->findCombos( uwires, vwires, ywires, urange, vrange, yrange );
      boundary_combos.emplace_back( std::move(combos) );
    }
  }
  
}
