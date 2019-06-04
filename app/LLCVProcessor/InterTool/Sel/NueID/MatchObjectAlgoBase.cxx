#ifndef __MATCHOBJECTALGOBASE_CXX__
#define __MATCHOBJECTALGOBASE_CXX__

#include "MatchObjectAlgoBase.h"
#include "LLCVBase/llcv_err.h"

#include <cassert>

namespace llcv {

  void MatchObjectAlgoBase::ClearMatch() {
    LLCV_DEBUG() << "start" << std::endl;
    for(auto& v : _object_vv)
      v.clear();

    _MatchBookKeeper.Reset();
    _seed_v.clear();

    for(auto& v : _object_per_plane_vv) v.clear();
    for(auto& v : _valid_plane_v) v = false;

    _object_ptr_v.clear();
    _object_id_to_plane_v.clear();
    _match_v.clear();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void MatchObjectAlgoBase::ClearEvent() {
    LLCV_DEBUG() << "start" << std::endl;
    ClearMatch();
    LLCV_DEBUG() << "end" << std::endl;
  }

  void MatchObjectAlgoBase::Register(const Object2D& obj,size_t plane) {
    _object_vv.at(plane).push_back(&obj);
    _valid_plane_v[plane] = true;
    LLCV_DEBUG() << "@plane=" << plane << " have " << _object_vv.at(plane).size() << " clusters" << std::endl;
  }

  const Object2D* MatchObjectAlgoBase::Object(size_t plane,size_t id) {
    return _object_vv.at(plane).at(id);
  }
  
  std::vector<std::vector<std::pair<size_t,size_t> > >
  MatchObjectAlgoBase::MatchObjects(std::vector<float>& score_v) {

    std::vector<std::vector<std::pair<size_t,size_t> > > match_vv;

    assert(_seed_v.empty());
    _seed_v.resize(3);
    for(size_t plane = 0; plane < 3; ++plane)  {
      _seed_v[plane] = _object_vv.at(plane).size();
      LLCV_DEBUG() << "@plane=" << plane << " sz=" << _seed_v[plane] << std::endl;
    }

    size_t offset=0;
    for(size_t plane=0;plane<3;++plane) {
      auto& object_per_plane_v = _object_per_plane_vv[plane];
      const auto& object_v = _object_vv.at(plane);
      for(size_t parid=0; parid<object_v.size(); ++parid) {
	const auto par = object_v[parid];
	object_per_plane_v.push_back(parid+offset);
	_object_ptr_v.push_back(par);
	_object_id_to_plane_v.emplace_back(plane,parid);
      }
      offset += object_v.size();
    }

    std::vector<std::vector<std::pair<size_t,size_t> > > temp_comb_vv, comb_vv;
    
    temp_comb_vv = larocv::PlaneClusterCombinations(_seed_v);
    comb_vv.reserve(temp_comb_vv.size());

    for(auto& temp_comb_v : temp_comb_vv) {
      if (temp_comb_v.size()==2) {
	if (!_valid_plane_v[temp_comb_v[0].first]) continue;
	if (!_valid_plane_v[temp_comb_v[1].first]) continue;
      }

      if (temp_comb_v.size()==3) {
	if (!_valid_plane_v[temp_comb_v[0].first]) continue;
	if (!_valid_plane_v[temp_comb_v[1].first]) continue;
	if (!_valid_plane_v[temp_comb_v[2].first]) continue;
      }
      comb_vv.emplace_back(std::move(temp_comb_v));
    }
    

    LLCV_DEBUG() << "Calculated " << comb_vv.size() << " combinations" << std::endl;
    
    for(const auto& comb_v : comb_vv) {
      LLCV_DEBUG() << "@comb sz= "<< comb_v.size() << std::endl;
      if (comb_v.size()==2) {

	LLCV_DEBUG() << "@2" << std::endl;
	auto plane0 = comb_v[0].first;
	auto plane1 = comb_v[1].first;
	
	LLCV_DEBUG() << "planes=(" << plane0 << "," << plane1 << ")" << std::endl;

	auto cid0 = comb_v[0].second;
	auto cid1 = comb_v[1].second;

	auto id0 = _object_per_plane_vv[plane0][cid0];
	auto id1 = _object_per_plane_vv[plane1][cid1];
	
	auto score = this->Match(*(_object_ptr_v.at(id0)),
				 *(_object_ptr_v.at(id1)));
	
	_match_v.clear();
	_match_v.resize(2);
	_match_v[0] = id0;
	_match_v[1] = id1;
	
	if (plane0 == 2 or plane1 == 2)
	  score *= _plane_two_boost;

	if (_require_plane2) {
	  if (plane0 != 2 and plane1 != 2) {
	    score = 0.0;
	  }
	}

	LLCV_DEBUG() << "score=" << score << std::endl;
	if (score<_threshold) continue;
	LLCV_DEBUG() << "...pass" << std::endl;

	_MatchBookKeeper.Match(_match_v,score);
      }
      else if (comb_v.size()==3 && _match_three_planes) {
	LLCV_DEBUG() << "@3" << std::endl;

	auto plane0 = comb_v[0].first;
	auto plane1 = comb_v[1].first;
	auto plane2 = comb_v[2].first;

	LLCV_DEBUG() << "planes=(" << plane0 << "," << plane1 << "," << plane2 << ")" << std::endl;

	auto cid0 = comb_v[0].second;
	auto cid1 = comb_v[1].second;
	auto cid2 = comb_v[2].second;

	auto id0 = _object_per_plane_vv[plane0][cid0];
	auto id1 = _object_per_plane_vv[plane1][cid1];
	auto id2 = _object_per_plane_vv[plane2][cid2];

	LLCV_DEBUG() << "ids=(" << id0 << "," << id1 << "," << id2 << ")" << std::endl;

	auto score = this->Match(*(_object_ptr_v.at(id0)),
				 *(_object_ptr_v.at(id1)),
				 *(_object_ptr_v.at(id2)));
	
	_match_v.clear();
	_match_v.resize(3);
	_match_v[0] = id0;
	_match_v[1] = id1;
	_match_v[2] = id2;

	score *= _three_planes_boost;

	LLCV_DEBUG() << "score=" << score << std::endl;

	if (score<_threshold) continue;

	LLCV_DEBUG() << "...pass" << std::endl;

	_MatchBookKeeper.Match(_match_v,score);
      }
      else {
	LLCV_INFO() << "match_three_planes=" << _match_three_planes << std::endl;
      }

    } // end combos

    auto result_v = _MatchBookKeeper.GetResultAndScore(score_v);

    LLCV_DEBUG() << "Found " << result_v.size() << " particle matches" << std::endl;

    match_vv.reserve(result_v.size());

    std::vector<std::pair<size_t,size_t> > res_match_v;
    for(size_t mid=0; mid<result_v.size(); ++mid) {
      const auto& match  = result_v[mid];
      res_match_v.clear();
      res_match_v.reserve(match.size()); // the number of matches contours
      for(auto par_id : match) 
	res_match_v.emplace_back(_object_id_to_plane_v[par_id]);
      match_vv.emplace_back(std::move(res_match_v));
    }

    return match_vv;
  } 

}


#endif
