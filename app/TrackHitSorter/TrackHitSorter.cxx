#include "TrackHitSorter.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/DetectorProperties.h"

namespace larlitecv {

  void TrackHitSorter::buildSortedHitList( const larlite::vertex& vertex,
					   const larlite::track& track,
					   const std::vector<larlite::hit>& hit_v,
					   const float max_radius,
					   std::vector<int>& hitmask_v ) {
    
    // geo utility
    const larutil::Geometry* geo = larutil::Geometry::GetME();
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float cm_per_tick = driftv*0.5;
    
    // convert track into line segments
    for (int p=0; p<3; p++) {
      path3d[p].clear();
      dist3d[p].clear();
      seg_v[p].clear();
      segdist_v[p].clear();
    }
    
    // get plane position of vertex
    Double_t vtx_xyz[3];
    vertex.XYZ( vtx_xyz );
    float vtx_x = vtx_xyz[0];
    std::vector<float> vtx_w(3);
    std::vector<geo2d::Vector<float>> vtx_pt_v;
    for (int p=0; p<3; p++) {
      try {
	  vtx_w[p] = 0.3*geo->NearestWire( vtx_xyz, p );
	} 
      //catch (const larutil::InvalidWireError& err) 
      catch (...)  {
	  std::cout << "Could _not_ find nearest vtx_xyz wire @p=" << p 
		    << " for (" << vtx_xyz[0] << "," << vtx_xyz[1] << "," << vtx_xyz[2] << ")" << std::endl;
	  vtx_w[p] = -1;
	}
      vtx_pt_v.push_back( geo2d::Vector<float>( vtx_w[p], vtx_x ) );
    }

    int numpts = track.NumberTrajectoryPoints();
    int ipt = 0;
    float dist_s[3] = {0.0};
    while ( ipt+1<numpts ) {
      
      const TVector3& here = track.LocationAtPoint( ipt );
      const TVector3& next = track.LocationAtPoint( ipt+1 );

      //float tick1 = here.X()/cm_per_tick + 3200;
      //float tick2 = next.X()/cm_per_tick + 3200;

      std::cout << "[ipt" << ipt << "] ";
      for (int p=0; p<3; p++) {
	
	float wire1 = -1;
	float wire2 = -1;
	try { 
	  wire1 = geo->NearestWire( here, p );
	} 
	//catch (const larutil::InvalidWireError& err) {
	catch (...) {
	  std::cout << "Could _not_ find nearest wire1 @p=" << p 
		    << " for (" << here[0] << "," << here[1] << "," << here[2] << ")" << std::endl;
	}

	try {
	  wire2 = geo->NearestWire( next, p );
	}
	//catch (const larutil::InvalidWireError& err) {
	catch (...) {
	  std::cout << "Could _not_ find nearest wire2 @p=" << p 
		    << " for (" << next[0] << "," << next[1] << "," << next[2] << ")" << std::endl;
	}
	
	if (wire1 < 0) {
	  std::cout << "What are we supossed to do? wire1 < 0 continue" << std::endl;
	  continue;
	}
	if (wire2 < 0) {
	  std::cout << "What are we supossed to do? wire2 < 0 continue" << std::endl;
	  continue;
	}

	geo2d::LineSegment<float> ls( wire1*0.3, here.X(), wire2*0.3, next.X() ); // tick and wire in cm units. we want the distance to be in cm.
	std::cout << "p" << p << "=[ (" << wire1*0.3 << "," << here.X() << ") (" << wire2*0.3 << "," << next.X() << ") ] ";
	if ( geo2d::length2(ls)>0 ) {
	  seg_v[p].emplace_back( std::move(ls) );
	  segdist_v[p].push_back( dist_s[p] );
	  dist_s[p] += sqrt(geo2d::length2(ls));

	  // store 3d segment
	  std::vector<float> seg3d(3);
	  if ( path3d[p].size()==0 ) {
	    // set vertex
	    for (int i=0; i<3; i++)
	      seg3d[i] = here(i);
	    path3d[p].push_back( seg3d );
	  }
	  for (int i=0; i<3; i++)
	    seg3d[i]   = next(i);
	  path3d[p].push_back( seg3d );

	  // store distance between current and last point
	  float seglen3 = 0.;
	  size_t pos = path3d[p].size()-1;
	  for (int i=0; i<3; i++) {
	    float dd = path3d[p].at(pos)[i]-path3d[p].at(pos-1)[i];
	    seglen3 += dd*dd;
	  }
	  seglen3 = sqrt(seglen3);
	  dist3d[p].push_back(seglen3);
	} // end positive length

      } // end plane

      std::cout << std::endl;
      ipt++;
    }
    
    // convert hit positions to same coordinate system
    std::vector< geo2d::Vector<float> > hitco_v[3];
    std::vector< const larlite::hit* > hitp_v[3];
    std::vector< int > hitidx_v[3];
    int ihit=-1;
    for ( auto const& hit : hit_v ) {
      ihit++;
      if ( hitmask_v[ihit]==0 ) continue; // masked
      
      float peakx = (2400+hit.PeakTime() - 3200)*cm_per_tick; // x position
      int plane = (int)hit.WireID().Plane;
      int wire  = (int)hit.WireID().Wire;
      
      //std::cout << "[hit " << ihit << "] x=" << peakx << " p=" << plane << " w=" << wire*0.3 << std::endl;

      geo2d::Vector<float> pt( wire*0.3, peakx );
      hitco_v[plane].emplace_back( std::move(pt) );
      hitp_v[plane].push_back( &hit );
      hitidx_v[plane].push_back( ihit );
    }

    std::cout << "converted hits: p0=" << hitco_v[0].size() << " p1=" << hitco_v[1].size() << " p2=" << hitco_v[2].size() << std::endl;

    // ok, now we associate
    for (int p=0; p<3; p++) {
      // dumb N^2 loop...
      int ihit=-1;
      for ( auto& pt : hitco_v[p] ) {
	ihit++;

	int iseg=-1;
	for ( auto& seg : seg_v[p] ) {
	  iseg++;
	  
	  // check if part of seg
	  geo2d::Vector<float> ab = seg.pt2 - seg.pt1;
	  float s = ab.ddot( pt-seg.pt1 )/geo2d::length2(ab);
	  if ( s<0 || s > 1.0 || std::isnan(s) )
	    continue;

	  geo2d::Vector<float> pt1; // point on line
	  pt1 = seg.pt1 + s*ab;
	  s = geo2d::dist(pt1,seg.pt1);
	  float r = geo2d::dist(pt,pt1);

	  // std::cout << "(" << ihit << "," << iseg << "): "
	  // 	    << " hit[" << pt.x << "," << pt.y << "] "
	  // 	    << " seg[(" << seg.pt1.x << "," << seg.pt1.y << ")->(" << seg.pt2.x << "," << seg.pt2.y << ")] "
	  // 	    << " s= "<< s << " r=" << r << " len(ab)=" << geo2d::length2(ab) << std::endl;
	  
	  
	  if ( r > max_radius ) {
	    continue;
	  }
	  float d = geo2d::dist( pt1, vtx_pt_v[p] );
	  
	  HitOrder path_ho( hitp_v[p].at(ihit), s+segdist_v[p][iseg], r );
	  HitOrder dist_ho( hitp_v[p].at(ihit), d, r );	  
	  pathordered[p].emplace_back( std::move(path_ho) );
	  distordered[p].emplace_back( std::move(dist_ho) );
	  hitmask_v[ hitidx_v[p].at(ihit) ] = 0; // mask out
	  break;
	}//end of loop over track segment in plane
      }//end of loop over hits in plane
    }//end of loop over plane hits

    // sort hits
    for (int p=0; p<3; p++) {
      std::sort( pathordered[p].begin(), pathordered[p].end() );
      std::sort( distordered[p].begin(), distordered[p].end() ); 
    }
  }


  void TrackHitSorter::getPathBinneddEdx( const float binstep, const float binwidth, std::vector< std::vector<float> >& dedx_per_plane ) {
    // for each plane, we move through 3d track, taking steps of binwidth and collecting hits [-binwidth,binwidth] from center position
    // these are in 3D, so we need to decide what s-values to collect for he projected hits
    const larutil::Geometry* geo = larutil::Geometry::GetME();
    //const larutil::LArProperties* larp = larutil::LArProperties::GetME();

    for (int p=2; p<3; p++) {

      // get track segment information. both 3d and 2d projections
      const std::vector< std::vector<float> >& plpath3d = path3d[p];
      const std::vector< float >& pldist3d = dist3d[p];
      const std::vector< geo2d::LineSegment<float> >& plseg_v = seg_v[p]; // length of seg
      //const std::vector< float >& pldist2d = segdist_v[p]; // coordinate

      // std::cout << "pldist2d.size()=" << pldist2d.size() << std::endl;
      // std::cout << "pldist3d.size()=" << pldist3d.size() << std::endl;

      bincenters_xyz[p].clear();
      
      float pathdist = 0.;
      for (auto const& d : pldist3d )
	pathdist += d;

      float dcenter = binstep;      
      while ( dcenter<pathdist ) {
      
	float dstart  = dcenter-binwidth;
	float dend    = dcenter+binwidth;
	std::vector<float> bincenter_xyz(3);
	
	if ( dstart<0 )
	  dstart = 0;
	if ( dend > pathdist )
	  dend = pathdist;
	
	float d = 0; // 3d distance
	float s = 0; // 2d distance
	float sbin_start = 0;
	float sbin_end   = 0;
	for (int iseg=0; iseg<(int)pldist3d.size(); iseg++) {
	  // get d-values spanned by the segment
	  // 3d range
	  float segd1 = d;
	  float segd2 = d+pldist3d[iseg];
	  // 2d (projected) range
	  float segslen = sqrt( geo2d::length2( plseg_v[iseg] ) );
	  float segs1 = s;
	  float segs2 = s+segslen;

	  // check if 3d segment is in bin
	  if ( segd2<dstart ) {
	    s = segs2;
	    d = segd2;
	    continue; // below bin range: move on
	  }
	  //if ( segd1>dend )
	  //  break; // our bin has moved past the segments
	  
	  // otherwise, segment straddles either dstart or dend
	  std::vector<float> segdir(3);
	  for (int i=0; i<3; i++)
	    segdir[i] = (plpath3d[iseg+1][i]-plpath3d[iseg][i])/pldist3d[iseg];
	  
	  if ( segd1<dstart && dstart<segd2 ) {
	    // straddles dstart
	    Double_t dstart3d[3];
	    for (int i=0; i<3; i++)
	      dstart3d[i] = plpath3d[iseg][i] + segdir[i]*(dstart-segd1);
	    // project start into plane
	    float wire = -1;
	    try {
	      wire = geo->NearestWire( dstart3d, p );
	    }
	    catch (...) {
	      wire = geo->Nwires(p)-1;
	    }
	    geo2d::Vector<float> start2d( wire*0.3, dstart3d[0] );
	    float dels = geo2d::dist( start2d, plseg_v[iseg].pt1 );
	    sbin_start = segs1 + dels;
	    //std::cout << "sbin_start update: " << segs1 << "+" << dels << " of " << segslen << std::endl;
	  }
	  
	  if ( segd1<dend && dend<=segd2 ) {
	    // straddles dend
	    Double_t dend3d[3];
	    for (int i=0; i<3; i++)
	      dend3d[i] = plpath3d[iseg][i] + segdir[i]*(dend-segd1);
	    // project start into plane
	    float wire = -1;
	    try {
	      wire = geo->NearestWire( dend3d, p );
	    }
	    catch (...) {
	      wire = geo->Nwires(p)-1;
	    }
	    geo2d::Vector<float> end2d( wire*0.3, dend3d[0] );
	    float dels = geo2d::dist( end2d, plseg_v[iseg].pt1 );
	    sbin_end= s + dels;
	  }

	  if ( segd1<dcenter && dcenter<=segd2 ) {
	    for (int i=0; i<3; i++)
	      bincenter_xyz[i] = plpath3d[iseg][i] + segdir[i]*(dcenter-segd1);
	    bincenters_xyz[p].push_back( bincenter_xyz );
	  }

	  // update
	  s = segs2;
	  d = segd2;
	}//end of loop over segments

	
	// now we sum over hits in the srange we found
	float q = 0;
	int nhits = 0;
	for ( auto const& hitho : pathordered[p] ) {
	  if ( sbin_start < hitho.s && hitho.s < sbin_end ) {
	    q += hitho.phit->Integral();
	    nhits++;
	  }
	}

	float cm   = dend-dstart;
	float MeV  = q2MeV( q, bincenter_xyz );
	float dqdx = q/cm;
	float dEdx = MeV/cm;

	// std::cout << "bincenter:" << dcenter
	// 	  << " centerxyz=(" << bincenter_xyz[0] << "," << bincenter_xyz[1] << "," << bincenter_xyz[2] << ") "
	// 	  << "dbin=[" << dstart << "," << dend << "] "
	// 	  << "sbin=[" << sbin_start << "," << sbin_end << "] "
	// 	  << " nhits=" << nhits 
	// 	  << " dqdx=" << dqdx
	// 	  << " dEdx=" << dEdx
	// 	  << std::endl;      
	
	dedx_per_plane[p].push_back( dEdx );
	
	dcenter += binstep;
      }// end of bincenter loop
      
    }//end of loop over planes
  }

  float TrackHitSorter::q2MeV( const float q, const std::vector<float>& xyz ) {

    //const float _fC_to_e = 6250.; // e- / fC
    const float _e_to_eV = 23.6;  // eV / e-
    const float _eV_to_MeV = 1e-6;// eV / MeV
    //const float _ADC_to_mV = 0.5; // ADC -> mV conversion from gain measurements

    const float _clocktick  = larutil::DetectorProperties::GetME()->SamplingRate() * 1.e-3;    
    //const float _tau        = larutil::LArProperties::GetME()->ElectronLifetime();
    const float _tau        = 1.0e6; // what is the units?
    const float cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    
    float _tick = xyz[0]/cm_per_tick; // from trigger (not abs scale)
    
    float _lifetime_corr = exp( _tick * _clocktick / _tau );

    float _elec_gain = -1.0;
    
    if(_ismc)
      _elec_gain     = 200; // MC -- MCC 8.3 value
    else
      _elec_gain     = 240; // Data -- MCC 8.3 value
    
    float _recomb_factor = 0.62;
    
    //double qcorr = ChargeCorrection(h.charge,h.w,h.t,resultShower.fDCosStart,resultShower.fXYZStart);
    double qcorr = 1.0*q; // skipping correction for now;

    float _electrons = qcorr * _elec_gain;

    float _dQ = _electrons * _lifetime_corr * _e_to_eV * _eV_to_MeV;

    float _dE = _dQ / _recomb_factor;

    return _dE;

  }// loop over all hits

  void TrackHitSorter::clear() {
    for (int p=0; p<3; p++) {
      path3d[p].clear();
      dist3d[p].clear();
      seg_v[p].clear();
      segdist_v[p].clear();
      pathordered[p].clear();
      distordered[p].clear();
      bincenters_xyz[p].clear();
    }
  }
  
  void TrackHitSorter::dump() const {
    std::cout << "=========================================================" << std::endl;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    //const larutil::Geometry* geo = larutil::Geometry::GetME();
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float cm_per_tick = driftv*0.5;
    
    for (int p=0; p<3; p++) {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "Hits on Plane " << p << std::endl;
      for (auto const& ho : pathordered[p] ) {
	const larlite::hit* phit = ho.phit;
	std::cout << "  (" << ho.s << "," << ho.r << ") x=" << (2400+phit->PeakTime() - 3200)*cm_per_tick << " w=" << phit->WireID().Wire << std::endl;
      }
    }
    std::cout << "=========================================================" << std::endl;    
  }
  
}
