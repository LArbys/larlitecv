#ifndef __SELMICHELID_CXX__
#define __SELMICHELID_CXX__

#include "SelMichelID.h"

#include <iomanip>
#include <sstream>
#include <cstdlib>

#include "Geo2D/Core/LineSegment.h"
#include "Base/DataFormatConstants.h"

#include "TVector2.h"

#include "InterTool_Util/InterImageUtils.h"

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

namespace llcv {

  void SelMichelID::Configure (const larcv::PSet &pset) {
    set_verbosity((msg::Level_t)pset.get<int>("Verbosity",2));
    LLCV_DEBUG() << "start" << std::endl;
    _cropx = 400;
    _cropy = 400;

    _extension_cutoff = pset.get<float>("ExtensionFraction",0.8);
    _Match.Configure(pset.get<larcv::PSet>("MatchIOU"));
    
    _twatch.Stop();
  
    _CosmicTag_v.clear();
    _CosmicTag_v.resize(3);

    _LineExtension_v.clear();
    _LineExtension_v.resize(3);

    _brem_dist = pset.get<float>("BremDistance",60);
    _brem_size = pset.get<int>("BremSize",6);
    _tradius = pset.get<float>("TruncatedRadius",3);
    _tsigma  = pset.get<float>("TruncatedSigma",0.5);
    
    _ShowerTools._TruncMean.setRadius(_tradius);

    LLCV_DEBUG() << "end" << std::endl;
  }

  void SelMichelID::Initialize() {
    _outtree = new TTree("SelMichelID","");
    AttachRSEV(_outtree);
    
    _outtree->Branch("vertex_x", &_vertex_x, "vertex_x/F");
    _outtree->Branch("vertex_y", &_vertex_y, "vertex_y/F");
    _outtree->Branch("vertex_z", &_vertex_z, "vertex_z/F");

    //
    // vertex information
    //
    _outtree->Branch("n_par"       , &_n_par       , "n_par/I");
    _outtree->Branch("n_vtx_cosmic", &_n_vtx_cosmic, "n_vtx_cosmic/I");
    _outtree->Branch("par_score_v" , &_par_score_v);

    _outtree->Branch("vtx_xing_U", &_vtx_xing_U, "vtx_xing_U/I");
    _outtree->Branch("vtx_xing_V", &_vtx_xing_V, "vtx_xing_V/I");
    _outtree->Branch("vtx_xing_Y", &_vtx_xing_Y, "vtx_xing_Y/I");

    _outtree->Branch("vtx_charge_U", &_vtx_charge_U, "vtx_charge_U/I");
    _outtree->Branch("vtx_charge_V", &_vtx_charge_V, "vtx_charge_V/I");
    _outtree->Branch("vtx_charge_Y", &_vtx_charge_Y, "vtx_charge_Y/I");

    _outtree->Branch("vtx_dead_U", &_vtx_dead_U, "vtx_dead_U/I");
    _outtree->Branch("vtx_dead_V", &_vtx_dead_V, "vtx_dead_V/I");
    _outtree->Branch("vtx_dead_Y", &_vtx_dead_Y, "vtx_dead_Y/I");

    _outtree->Branch("vtx_linelen_U_v", &_vtx_linelen_U_v);
    _outtree->Branch("vtx_linelen_V_v", &_vtx_linelen_V_v);
    _outtree->Branch("vtx_linelen_Y_v", &_vtx_linelen_Y_v);

    _outtree->Branch("vtx_linefrac_U_v", &_vtx_linefrac_U_v);
    _outtree->Branch("vtx_linefrac_V_v", &_vtx_linefrac_V_v);
    _outtree->Branch("vtx_linefrac_Y_v", &_vtx_linefrac_Y_v);
    
    _outtree->Branch("edge_dist_v"                , &_edge_dist_v);
    _outtree->Branch("edge_n_cosmic_v"            , &_edge_n_cosmic_v);
    _outtree->Branch("edge_cosmic_vtx_dist_v"     , &_edge_cosmic_vtx_dist_v);
    _outtree->Branch("edge_cosmic_end_vtx_dist_v" , &_edge_cosmic_end_vtx_dist_v);
    _outtree->Branch("edge_cosmic_end_vtx_valid_v", &_edge_cosmic_end_vtx_valid_v);


    //
    // 3D information
    //
    _outtree->Branch("par1_theta"    , &_par1_theta  , "par1_theta/F");
    _outtree->Branch("par1_phi"      , &_par1_phi    , "par1_phi/F");
    _outtree->Branch("par1_length"   , &_par1_length , "par1_length/F");
    _outtree->Branch("par1_score"    , &_par1_score  , "par1_score/F");
    _outtree->Branch("par1_dx1"      , &_par1_dx1    , "par1_dx1/F");
    _outtree->Branch("par1_dy1"      , &_par1_dy1    , "par1_dy1/F");
    _outtree->Branch("par1_dz1"      , &_par1_dz1    , "par1_dz1/F");
    _outtree->Branch("par1_dx2"      , &_par1_dx2    , "par1_dx2/F");
    _outtree->Branch("par1_dy2"      , &_par1_dy2    , "par1_dy2/F");
    _outtree->Branch("par1_dz2"      , &_par1_dz2    , "par1_dz2/F");
    _outtree->Branch("par1_nplanes"  , &_par1_nplanes, "par1_nplanes/I");
    _outtree->Branch("par1_planes_v" , &_par1_planes_v);
    _outtree->Branch("par1_xdead_v"  , &_par1_xdead_v);
    _outtree->Branch("par1_end_pt_v" , &_par1_end_pt_v);

    _outtree->Branch("par1_cosmic_dist_v"        , &_par1_cosmic_dist_v);
    _outtree->Branch("par1_cosmic_dist_end_v"    , &_par1_cosmic_dist_end_v);
    _outtree->Branch("par1_cosmic_end_dist_v"    , &_par1_cosmic_end_dist_v);
    _outtree->Branch("par1_cosmic_end_dist_end_v", &_par1_cosmic_end_dist_end_v);


    _outtree->Branch("par2_theta"    , &_par2_theta  , "par2_theta/F");
    _outtree->Branch("par2_phi"      , &_par2_phi    , "par2_phi/F");
    _outtree->Branch("par2_length"   , &_par2_length , "par2_length/F");
    _outtree->Branch("par2_score"    , &_par2_score  , "par2_score/F");
    _outtree->Branch("par2_dx1"      , &_par2_dx1    , "par2_dx1/F");
    _outtree->Branch("par2_dy1"      , &_par2_dy1    , "par2_dy1/F");
    _outtree->Branch("par2_dz1"      , &_par2_dz1    , "par2_dz1/F");
    _outtree->Branch("par2_dx2"      , &_par2_dx2    , "par2_dx2/F");
    _outtree->Branch("par2_dy2"      , &_par2_dy2    , "par2_dy2/F");
    _outtree->Branch("par2_dz2"      , &_par2_dz2    , "par2_dz2/F");
    _outtree->Branch("par2_nplanes"  , &_par2_nplanes, "par2_nplanes/I");
    _outtree->Branch("par2_planes_v" , &_par2_planes_v);
    _outtree->Branch("par2_xdead_v"  , &_par2_xdead_v);
    _outtree->Branch("par2_end_pt_v" , &_par2_end_pt_v);
    _outtree->Branch("par2_cosmic_dist_v"        , &_par2_cosmic_dist_v);
    _outtree->Branch("par2_cosmic_dist_end_v"    , &_par2_cosmic_dist_end_v);
    _outtree->Branch("par2_cosmic_end_dist_v"    , &_par2_cosmic_end_dist_v);
    _outtree->Branch("par2_cosmic_end_dist_end_v", &_par2_cosmic_end_dist_end_v);


    //
    // 2D information
    //
    _outtree->Branch("par1_n_polygons_U", &_par1_n_polygons_U, "par1_n_polygons_U/I");
    _outtree->Branch("par1_n_polygons_V", &_par1_n_polygons_V, "par1_n_polygons_V/I");
    _outtree->Branch("par1_n_polygons_Y", &_par1_n_polygons_Y, "par1_n_polygons_Y/I");
    _outtree->Branch("par2_n_polygons_U", &_par2_n_polygons_U, "par2_n_polygons_U/I");
    _outtree->Branch("par2_n_polygons_V", &_par2_n_polygons_V, "par2_n_polygons_V/I");
    _outtree->Branch("par2_n_polygons_Y", &_par2_n_polygons_Y, "par2_n_polygons_Y/I");

    _outtree->Branch("par1_linelength_U", &_par1_linelength_U, "par1_linelength_U/F");
    _outtree->Branch("par1_linelength_V", &_par1_linelength_V, "par1_linelength_V/F");
    _outtree->Branch("par1_linelength_Y", &_par1_linelength_Y, "par1_linelength_Y/F");
    _outtree->Branch("par2_linelength_U", &_par2_linelength_U, "par2_linelength_U/F");
    _outtree->Branch("par2_linelength_V", &_par2_linelength_V, "par2_linelength_V/F");
    _outtree->Branch("par2_linelength_Y", &_par2_linelength_Y, "par2_linelength_Y/F");

    _outtree->Branch("par1_linefrac_U", &_par1_linefrac_U, "par1_linefrac_U/F");
    _outtree->Branch("par1_linefrac_V", &_par1_linefrac_V, "par1_linefrac_V/F");
    _outtree->Branch("par1_linefrac_Y", &_par1_linefrac_Y, "par1_linefrac_Y/F");
    _outtree->Branch("par2_linefrac_U", &_par2_linefrac_U, "par2_linefrac_U/F");
    _outtree->Branch("par2_linefrac_V", &_par2_linefrac_V, "par2_linefrac_V/F");
    _outtree->Branch("par2_linefrac_Y", &_par2_linefrac_Y, "par2_linefrac_Y/F");

    _outtree->Branch("par1_linefrac_empty_U", &_par1_linefrac_empty_U, "par1_linefrac_empty_U/F");
    _outtree->Branch("par1_linefrac_empty_V", &_par1_linefrac_empty_V, "par1_linefrac_empty_V/F");
    _outtree->Branch("par1_linefrac_empty_Y", &_par1_linefrac_empty_Y, "par1_linefrac_empty_Y/F");
    _outtree->Branch("par2_linefrac_empty_U", &_par2_linefrac_empty_U, "par2_linefrac_empty_U/F");
    _outtree->Branch("par2_linefrac_empty_V", &_par2_linefrac_empty_V, "par2_linefrac_empty_V/F");
    _outtree->Branch("par2_linefrac_empty_Y", &_par2_linefrac_empty_Y, "par2_linefrac_empty_Y/F");

    _outtree->Branch("par1_linedx_U", &_par1_linedx_U, "par1_linedx_U/F");
    _outtree->Branch("par1_linedx_V", &_par1_linedx_V, "par1_linedx_V/F");
    _outtree->Branch("par1_linedx_Y", &_par1_linedx_Y, "par1_linedx_Y/F");
    _outtree->Branch("par2_linedx_U", &_par2_linedx_U, "par2_linedx_U/F");
    _outtree->Branch("par2_linedx_V", &_par2_linedx_V, "par2_linedx_V/F");
    _outtree->Branch("par2_linedx_Y", &_par2_linedx_Y, "par2_linedx_Y/F");
    
    _outtree->Branch("par1_linedy_U", &_par1_linedy_U, "par1_linedy_U/F");
    _outtree->Branch("par1_linedy_V", &_par1_linedy_V, "par1_linedy_V/F");
    _outtree->Branch("par1_linedy_Y", &_par1_linedy_Y, "par1_linedy_Y/F");
    _outtree->Branch("par2_linedy_U", &_par2_linedy_U, "par2_linedy_U/F");
    _outtree->Branch("par2_linedy_V", &_par2_linedy_V, "par2_linedy_V/F");
    _outtree->Branch("par2_linedy_Y", &_par2_linedy_Y, "par2_linedy_Y/F");

    _outtree->Branch("par1_line_vtx_density_U", &_par1_line_vtx_density_U, "par1_line_vtx_density_U/F");
    _outtree->Branch("par1_line_vtx_density_V", &_par1_line_vtx_density_V, "par1_line_vtx_density_V/F");
    _outtree->Branch("par1_line_vtx_density_Y", &_par1_line_vtx_density_Y, "par1_line_vtx_density_Y/F");
    _outtree->Branch("par2_line_vtx_density_U", &_par2_line_vtx_density_U, "par2_line_vtx_density_U/F");
    _outtree->Branch("par2_line_vtx_density_V", &_par2_line_vtx_density_V, "par2_line_vtx_density_V/F");
    _outtree->Branch("par2_line_vtx_density_Y", &_par2_line_vtx_density_Y, "par2_line_vtx_density_Y/F");

    _outtree->Branch("par1_line_vtx_coverage_U", &_par1_line_vtx_coverage_U, "par1_line_vtx_coverage_U/F");
    _outtree->Branch("par1_line_vtx_coverage_V", &_par1_line_vtx_coverage_V, "par1_line_vtx_coverage_V/F");
    _outtree->Branch("par1_line_vtx_coverage_Y", &_par1_line_vtx_coverage_Y, "par1_line_vtx_coverage_Y/F");
    _outtree->Branch("par2_line_vtx_coverage_U", &_par2_line_vtx_coverage_U, "par2_line_vtx_coverage_U/F");
    _outtree->Branch("par2_line_vtx_coverage_V", &_par2_line_vtx_coverage_V, "par2_line_vtx_coverage_V/F");
    _outtree->Branch("par2_line_vtx_coverage_Y", &_par2_line_vtx_coverage_Y, "par2_line_vtx_coverage_Y/F");

    _outtree->Branch("par1_line_vtx_charge_U", &_par1_line_vtx_charge_U, "par1_line_vtx_charge_U/F");
    _outtree->Branch("par1_line_vtx_charge_V", &_par1_line_vtx_charge_V, "par1_line_vtx_charge_V/F");
    _outtree->Branch("par1_line_vtx_charge_Y", &_par1_line_vtx_charge_Y, "par1_line_vtx_charge_Y/F");
    _outtree->Branch("par2_line_vtx_charge_U", &_par2_line_vtx_charge_U, "par2_line_vtx_charge_U/F");
    _outtree->Branch("par2_line_vtx_charge_V", &_par2_line_vtx_charge_V, "par2_line_vtx_charge_V/F");
    _outtree->Branch("par2_line_vtx_charge_Y", &_par2_line_vtx_charge_Y, "par2_line_vtx_charge_Y/F");

    _outtree->Branch("par1_line_mean_dist_U", &_par1_line_mean_dist_U, "par1_line_mean_dist_U/F");
    _outtree->Branch("par1_line_mean_dist_V", &_par1_line_mean_dist_V, "par1_line_mean_dist_V/F");
    _outtree->Branch("par1_line_mean_dist_Y", &_par1_line_mean_dist_Y, "par1_line_mean_dist_Y/F");
    _outtree->Branch("par2_line_mean_dist_U", &_par2_line_mean_dist_U, "par2_line_mean_dist_U/F");
    _outtree->Branch("par2_line_mean_dist_V", &_par2_line_mean_dist_V, "par2_line_mean_dist_V/F");
    _outtree->Branch("par2_line_mean_dist_Y", &_par2_line_mean_dist_Y, "par2_line_mean_dist_Y/F");

    _outtree->Branch("par1_line_max_dist_U", &_par1_line_max_dist_U, "par1_line_max_dist_U/F");
    _outtree->Branch("par1_line_max_dist_V", &_par1_line_max_dist_V, "par1_line_max_dist_V/F");
    _outtree->Branch("par1_line_max_dist_Y", &_par1_line_max_dist_Y, "par1_line_max_dist_Y/F");
    _outtree->Branch("par2_line_max_dist_U", &_par2_line_max_dist_U, "par2_line_max_dist_U/F");
    _outtree->Branch("par2_line_max_dist_V", &_par2_line_max_dist_V, "par2_line_max_dist_V/F");
    _outtree->Branch("par2_line_max_dist_Y", &_par2_line_max_dist_Y, "par2_line_max_dist_Y/F");

    _outtree->Branch("par1_line_first_half_linefrac_U", &_par1_line_first_half_linefrac_U, "par1_line_first_half_linefrac_U");
    _outtree->Branch("par1_line_first_half_linefrac_V", &_par1_line_first_half_linefrac_V, "par1_line_first_half_linefrac_V");
    _outtree->Branch("par1_line_first_half_linefrac_Y", &_par1_line_first_half_linefrac_Y, "par1_line_first_half_linefrac_Y");
    _outtree->Branch("par2_line_first_half_linefrac_U", &_par2_line_first_half_linefrac_U, "par2_line_first_half_linefrac_U");
    _outtree->Branch("par2_line_first_half_linefrac_V", &_par2_line_first_half_linefrac_V, "par2_line_first_half_linefrac_V");
    _outtree->Branch("par2_line_first_half_linefrac_Y", &_par2_line_first_half_linefrac_Y, "par2_line_first_half_linefrac_Y");

    _outtree->Branch("par1_line_second_half_linefrac_U", &_par1_line_second_half_linefrac_U, "par1_line_second_half_linefrac_U");
    _outtree->Branch("par1_line_second_half_linefrac_V", &_par1_line_second_half_linefrac_V, "par1_line_second_half_linefrac_V");
    _outtree->Branch("par1_line_second_half_linefrac_Y", &_par1_line_second_half_linefrac_Y, "par1_line_second_half_linefrac_Y");
    _outtree->Branch("par2_line_second_half_linefrac_U", &_par2_line_second_half_linefrac_U, "par2_line_second_half_linefrac_U");
    _outtree->Branch("par2_line_second_half_linefrac_V", &_par2_line_second_half_linefrac_V, "par2_line_second_half_linefrac_V");
    _outtree->Branch("par2_line_second_half_linefrac_Y", &_par2_line_second_half_linefrac_Y, "par2_line_second_half_linefrac_Y");

    _outtree->Branch("par1_triangle_mid_to_edge_U", &_par1_triangle_mid_to_edge_U, "par1_triangle_mid_to_edge_U/F");
    _outtree->Branch("par1_triangle_mid_to_edge_V", &_par1_triangle_mid_to_edge_V, "par1_triangle_mid_to_edge_V/F");
    _outtree->Branch("par1_triangle_mid_to_edge_Y", &_par1_triangle_mid_to_edge_Y, "par1_triangle_mid_to_edge_Y/F");
    _outtree->Branch("par2_triangle_mid_to_edge_U", &_par2_triangle_mid_to_edge_U, "par2_triangle_mid_to_edge_U/F");
    _outtree->Branch("par2_triangle_mid_to_edge_V", &_par2_triangle_mid_to_edge_V, "par2_triangle_mid_to_edge_V/F");
    _outtree->Branch("par2_triangle_mid_to_edge_Y", &_par2_triangle_mid_to_edge_Y, "par2_triangle_mid_to_edge_Y/F");

    _outtree->Branch("par1_triangle_height_U", &_par1_triangle_height_U, "par1_triangle_height_U/F");
    _outtree->Branch("par1_triangle_height_V", &_par1_triangle_height_V, "par1_triangle_height_V/F");
    _outtree->Branch("par1_triangle_height_Y", &_par1_triangle_height_Y, "par1_triangle_height_Y/F");
    _outtree->Branch("par2_triangle_height_U", &_par2_triangle_height_U, "par2_triangle_height_U/F");
    _outtree->Branch("par2_triangle_height_V", &_par2_triangle_height_V, "par2_triangle_height_V/F");
    _outtree->Branch("par2_triangle_height_Y", &_par2_triangle_height_Y, "par2_triangle_height_Y/F");

    _outtree->Branch("par1_triangle_emptyarearatio_U", &_par1_triangle_emptyarearatio_U, "par1_triangle_emptyarearatio_U/F");
    _outtree->Branch("par1_triangle_emptyarearatio_V", &_par1_triangle_emptyarearatio_V, "par1_triangle_emptyarearatio_V/F");
    _outtree->Branch("par1_triangle_emptyarearatio_Y", &_par1_triangle_emptyarearatio_Y, "par1_triangle_emptyarearatio_Y/F");
    _outtree->Branch("par2_triangle_emptyarearatio_U", &_par2_triangle_emptyarearatio_U, "par2_triangle_emptyarearatio_U/F");
    _outtree->Branch("par2_triangle_emptyarearatio_V", &_par2_triangle_emptyarearatio_V, "par2_triangle_emptyarearatio_V/F");
    _outtree->Branch("par2_triangle_emptyarearatio_Y", &_par2_triangle_emptyarearatio_Y, "par2_triangle_emptyarearatio_Y/F");

    _outtree->Branch("par1_triangle_emptyarea_U", &_par1_triangle_emptyarea_U, "par1_triangle_emptyarea_U/F");
    _outtree->Branch("par1_triangle_emptyarea_V", &_par1_triangle_emptyarea_V, "par1_triangle_emptyarea_V/F");
    _outtree->Branch("par1_triangle_emptyarea_Y", &_par1_triangle_emptyarea_Y, "par1_triangle_emptyarea_Y/F");
    _outtree->Branch("par2_triangle_emptyarea_U", &_par2_triangle_emptyarea_U, "par2_triangle_emptyarea_U/F");
    _outtree->Branch("par2_triangle_emptyarea_V", &_par2_triangle_emptyarea_V, "par2_triangle_emptyarea_V/F");
    _outtree->Branch("par2_triangle_emptyarea_Y", &_par2_triangle_emptyarea_Y, "par2_triangle_emptyarea_Y/F");

    _outtree->Branch("par1_triangle_baselength_U", &_par1_triangle_baselength_U, "par1_triangle_baselength_U/F");
    _outtree->Branch("par1_triangle_baselength_V", &_par1_triangle_baselength_V, "par1_triangle_baselength_V/F");
    _outtree->Branch("par1_triangle_baselength_Y", &_par1_triangle_baselength_Y, "par1_triangle_baselength_Y/F");
    _outtree->Branch("par2_triangle_baselength_U", &_par2_triangle_baselength_U, "par2_triangle_baselength_U/F");
    _outtree->Branch("par2_triangle_baselength_V", &_par2_triangle_baselength_V, "par2_triangle_baselength_V/F");
    _outtree->Branch("par2_triangle_baselength_Y", &_par2_triangle_baselength_Y, "par2_triangle_baselength_Y/F");

    _outtree->Branch("par1_triangle_area_U", &_par1_triangle_area_U, "par1_triangle_area_U/F");
    _outtree->Branch("par1_triangle_area_V", &_par1_triangle_area_V, "par1_triangle_area_V/F");
    _outtree->Branch("par1_triangle_area_Y", &_par1_triangle_area_Y, "par1_triangle_area_Y/F");
    _outtree->Branch("par2_triangle_area_U", &_par2_triangle_area_U, "par2_triangle_area_U/F");
    _outtree->Branch("par2_triangle_area_V", &_par2_triangle_area_V, "par2_triangle_area_V/F");
    _outtree->Branch("par2_triangle_area_Y", &_par2_triangle_area_Y, "par2_triangle_area_Y/F");

    _outtree->Branch("par1_triangle_brem_U", &_par1_triangle_brem_U, "par1_triangle_brem_U/F");
    _outtree->Branch("par1_triangle_brem_V", &_par1_triangle_brem_V, "par1_triangle_brem_V/F");
    _outtree->Branch("par1_triangle_brem_Y", &_par1_triangle_brem_Y, "par1_triangle_brem_Y/F");
    _outtree->Branch("par2_triangle_brem_U", &_par2_triangle_brem_U, "par2_triangle_brem_U/F");
    _outtree->Branch("par2_triangle_brem_V", &_par2_triangle_brem_V, "par2_triangle_brem_V/F");
    _outtree->Branch("par2_triangle_brem_Y", &_par2_triangle_brem_Y, "par2_triangle_brem_Y/F");

    _outtree->Branch("par1_triangle_coverage_U", &_par1_triangle_coverage_U, "par1_triangle_coverage_U/F");
    _outtree->Branch("par1_triangle_coverage_V", &_par1_triangle_coverage_V, "par1_triangle_coverage_V/F");
    _outtree->Branch("par1_triangle_coverage_Y", &_par1_triangle_coverage_Y, "par1_triangle_coverage_Y/F");
    _outtree->Branch("par2_triangle_coverage_U", &_par2_triangle_coverage_U, "par2_triangle_coverage_U/F");
    _outtree->Branch("par2_triangle_coverage_V", &_par2_triangle_coverage_V, "par2_triangle_coverage_V/F");
    _outtree->Branch("par2_triangle_coverage_Y", &_par2_triangle_coverage_Y, "par2_triangle_coverage_Y/F");

    _outtree->Branch("par1_brem_triangle_coverage_U", &_par1_brem_triangle_coverage_U, "par1_brem_triangle_coverage_U/F");
    _outtree->Branch("par1_brem_triangle_coverage_V", &_par1_brem_triangle_coverage_V, "par1_brem_triangle_coverage_V/F");
    _outtree->Branch("par1_brem_triangle_coverage_Y", &_par1_brem_triangle_coverage_Y, "par1_brem_triangle_coverage_Y/F");
    _outtree->Branch("par2_brem_triangle_coverage_U", &_par2_brem_triangle_coverage_U, "par2_brem_triangle_coverage_U/F");
    _outtree->Branch("par2_brem_triangle_coverage_V", &_par2_brem_triangle_coverage_V, "par2_brem_triangle_coverage_V/F");
    _outtree->Branch("par2_brem_triangle_coverage_Y", &_par2_brem_triangle_coverage_Y, "par2_brem_triangle_coverage_Y/F");

    _outtree->Branch("par1_brem_idx_U", &_par1_brem_idx_U, "par1_brem_idx_U/I");
    _outtree->Branch("par1_brem_idx_V", &_par1_brem_idx_V, "par1_brem_idx_V/I");
    _outtree->Branch("par1_brem_idx_Y", &_par1_brem_idx_Y, "par1_brem_idx_Y/I");
    _outtree->Branch("par2_brem_idx_U", &_par2_brem_idx_U, "par2_brem_idx_U/I");
    _outtree->Branch("par2_brem_idx_V", &_par2_brem_idx_V, "par2_brem_idx_V/I");
    _outtree->Branch("par2_brem_idx_Y", &_par2_brem_idx_Y, "par2_brem_idx_Y/I");

    _outtree->Branch("par1_expand_charge_U", &_par1_expand_charge_U, "par1_expand_charge_U/F");
    _outtree->Branch("par1_expand_charge_V", &_par1_expand_charge_V, "par1_expand_charge_V/F");
    _outtree->Branch("par1_expand_charge_Y", &_par1_expand_charge_Y, "par1_expand_charge_Y/F");
    _outtree->Branch("par2_expand_charge_U", &_par2_expand_charge_U, "par2_expand_charge_U/F");
    _outtree->Branch("par2_expand_charge_V", &_par2_expand_charge_V, "par2_expand_charge_V/F");
    _outtree->Branch("par2_expand_charge_Y", &_par2_expand_charge_Y, "par2_expand_charge_Y/F");

    _outtree->Branch("par1_dqdx_U", &_par1_dqdx_U, "par1_dqdx_U/F");
    _outtree->Branch("par1_dqdx_V", &_par1_dqdx_V, "par1_dqdx_V/F");
    _outtree->Branch("par1_dqdx_Y", &_par1_dqdx_Y, "par1_dqdx_Y/F");
    _outtree->Branch("par2_dqdx_U", &_par2_dqdx_U, "par2_dqdx_U/F");
    _outtree->Branch("par2_dqdx_V", &_par2_dqdx_V, "par2_dqdx_V/F");
    _outtree->Branch("par2_dqdx_Y", &_par2_dqdx_Y, "par2_dqdx_Y/F");

    _outtree->Branch("par1_dqdx_step_U", &_par1_dqdx_step_U, "par1_dqdx_step_U/F");
    _outtree->Branch("par1_dqdx_step_V", &_par1_dqdx_step_V, "par1_dqdx_step_V/F");
    _outtree->Branch("par1_dqdx_step_Y", &_par1_dqdx_step_Y, "par1_dqdx_step_Y/F");
    _outtree->Branch("par2_dqdx_step_U", &_par2_dqdx_step_U, "par2_dqdx_step_U/F");
    _outtree->Branch("par2_dqdx_step_V", &_par2_dqdx_step_V, "par2_dqdx_step_V/F");
    _outtree->Branch("par2_dqdx_step_Y", &_par2_dqdx_step_Y, "par2_dqdx_step_Y/F");

    _outtree->Branch("par1_dqdx_pitch_U", &_par1_dqdx_pitch_U, "par1_dqdx_pitch_U/F");
    _outtree->Branch("par1_dqdx_pitch_V", &_par1_dqdx_pitch_V, "par1_dqdx_pitch_V/F");
    _outtree->Branch("par1_dqdx_pitch_Y", &_par1_dqdx_pitch_Y, "par1_dqdx_pitch_Y/F");
    _outtree->Branch("par2_dqdx_pitch_U", &_par2_dqdx_pitch_U, "par2_dqdx_pitch_U/F");
    _outtree->Branch("par2_dqdx_pitch_V", &_par2_dqdx_pitch_V, "par2_dqdx_pitch_V/F");
    _outtree->Branch("par2_dqdx_pitch_Y", &_par2_dqdx_pitch_Y, "par2_dqdx_pitch_Y/F");

    _outtree->Branch("par1_dqdx_U_v", &_par1_dqdx_U_v);
    _outtree->Branch("par1_dqdx_V_v", &_par1_dqdx_V_v);
    _outtree->Branch("par1_dqdx_Y_v", &_par1_dqdx_Y_v);
    _outtree->Branch("par2_dqdx_U_v", &_par2_dqdx_U_v);
    _outtree->Branch("par2_dqdx_V_v", &_par2_dqdx_V_v);
    _outtree->Branch("par2_dqdx_Y_v", &_par2_dqdx_Y_v);

    _outtree->Branch("par1_tdqdx_U_v", &_par1_tdqdx_U_v);
    _outtree->Branch("par1_tdqdx_V_v", &_par1_tdqdx_V_v);
    _outtree->Branch("par1_tdqdx_Y_v", &_par1_tdqdx_Y_v);
    _outtree->Branch("par2_tdqdx_U_v", &_par2_tdqdx_U_v);
    _outtree->Branch("par2_tdqdx_V_v", &_par2_tdqdx_V_v);
    _outtree->Branch("par2_tdqdx_Y_v", &_par2_tdqdx_Y_v);

    _outtree->Branch("par1_length3d_U", &_par1_length3d_U, "par1_length3d_U/F");
    _outtree->Branch("par1_length3d_V", &_par1_length3d_V, "par1_length3d_V/F");
    _outtree->Branch("par1_length3d_Y", &_par1_length3d_Y, "par1_length3d_Y/F");
    _outtree->Branch("par2_length3d_U", &_par2_length3d_U, "par2_length3d_U/F");
    _outtree->Branch("par2_length3d_V", &_par2_length3d_V, "par2_length3d_V/F");
    _outtree->Branch("par2_length3d_Y", &_par2_length3d_Y, "par2_length3d_Y/F");

    _outtree->Branch("par1_showerfrac_U", &_par1_showerfrac_U, "par1_showerfrac_U/F");
    _outtree->Branch("par1_showerfrac_V", &_par1_showerfrac_V, "par1_showerfrac_V/F");
    _outtree->Branch("par1_showerfrac_Y", &_par1_showerfrac_Y, "par1_showerfrac_Y/F");
    _outtree->Branch("par2_showerfrac_U", &_par2_showerfrac_U, "par2_showerfrac_U/F");
    _outtree->Branch("par2_showerfrac_V", &_par2_showerfrac_V, "par2_showerfrac_V/F");
    _outtree->Branch("par2_showerfrac_Y", &_par2_showerfrac_Y, "par2_showerfrac_Y/F");
    
    _outtree->Branch("par1_numberdefects_U_v", &_par1_numberdefects_U_v);
    _outtree->Branch("par1_numberdefects_V_v", &_par1_numberdefects_V_v);
    _outtree->Branch("par1_numberdefects_Y_v", &_par1_numberdefects_Y_v);
    _outtree->Branch("par2_numberdefects_U_v", &_par2_numberdefects_U_v);
    _outtree->Branch("par2_numberdefects_V_v", &_par2_numberdefects_V_v);
    _outtree->Branch("par2_numberdefects_Y_v", &_par2_numberdefects_Y_v);
    
    _outtree->Branch("par1_numberdefects_ns_U_v", &_par1_numberdefects_ns_U_v);
    _outtree->Branch("par1_numberdefects_ns_V_v", &_par1_numberdefects_ns_V_v);
    _outtree->Branch("par1_numberdefects_ns_Y_v", &_par1_numberdefects_ns_Y_v);
    _outtree->Branch("par2_numberdefects_ns_U_v", &_par2_numberdefects_ns_U_v);
    _outtree->Branch("par2_numberdefects_ns_V_v", &_par2_numberdefects_ns_V_v);
    _outtree->Branch("par2_numberdefects_ns_Y_v", &_par2_numberdefects_ns_Y_v);
    
    _outtree->Branch("par1_largestdefect_U_v", &_par1_largestdefect_U_v);
    _outtree->Branch("par1_largestdefect_V_v", &_par1_largestdefect_V_v);
    _outtree->Branch("par1_largestdefect_Y_v", &_par1_largestdefect_Y_v);
    _outtree->Branch("par2_largestdefect_U_v", &_par2_largestdefect_U_v);
    _outtree->Branch("par2_largestdefect_V_v", &_par2_largestdefect_V_v);
    _outtree->Branch("par2_largestdefect_Y_v", &_par2_largestdefect_Y_v);

    _outtree->Branch("par1_smallestdefect_U_v", &_par1_smallestdefect_U_v);
    _outtree->Branch("par1_smallestdefect_V_v", &_par1_smallestdefect_V_v);
    _outtree->Branch("par1_smallestdefect_Y_v", &_par1_smallestdefect_Y_v);
    _outtree->Branch("par2_smallestdefect_U_v", &_par2_smallestdefect_U_v);
    _outtree->Branch("par2_smallestdefect_V_v", &_par2_smallestdefect_V_v);
    _outtree->Branch("par2_smallestdefect_Y_v", &_par2_smallestdefect_Y_v);

    _outtree->Branch("par1_largestdefect_ns_U_v", &_par1_largestdefect_ns_U_v);
    _outtree->Branch("par1_largestdefect_ns_V_v", &_par1_largestdefect_ns_V_v);
    _outtree->Branch("par1_largestdefect_ns_Y_v", &_par1_largestdefect_ns_Y_v);
    _outtree->Branch("par2_largestdefect_ns_U_v", &_par2_largestdefect_ns_U_v);
    _outtree->Branch("par2_largestdefect_ns_V_v", &_par2_largestdefect_ns_V_v);
    _outtree->Branch("par2_largestdefect_ns_Y_v", &_par2_largestdefect_ns_Y_v);
    
    _outtree->Branch("par1_smallestdefect_ns_U_v", &_par1_smallestdefect_ns_U_v);
    _outtree->Branch("par1_smallestdefect_ns_V_v", &_par1_smallestdefect_ns_V_v);
    _outtree->Branch("par1_smallestdefect_ns_Y_v", &_par1_smallestdefect_ns_Y_v);
    _outtree->Branch("par2_smallestdefect_ns_U_v", &_par2_smallestdefect_ns_U_v);
    _outtree->Branch("par2_smallestdefect_ns_V_v", &_par2_smallestdefect_ns_V_v);
    _outtree->Branch("par2_smallestdefect_ns_Y_v", &_par2_smallestdefect_ns_Y_v); 

    _outtree->Branch("par1_emptyarearatio_U_v", &_par1_emptyarearatio_U_v);
    _outtree->Branch("par1_emptyarearatio_V_v", &_par1_emptyarearatio_V_v);
    _outtree->Branch("par1_emptyarearatio_Y_v", &_par1_emptyarearatio_Y_v);
    _outtree->Branch("par2_emptyarearatio_U_v", &_par2_emptyarearatio_U_v);
    _outtree->Branch("par2_emptyarearatio_V_v", &_par2_emptyarearatio_V_v);
    _outtree->Branch("par2_emptyarearatio_Y_v", &_par2_emptyarearatio_Y_v);

    _outtree->Branch("par1_emptyarea_U_v", &_par1_emptyarea_U_v);
    _outtree->Branch("par1_emptyarea_V_v", &_par1_emptyarea_V_v);
    _outtree->Branch("par1_emptyarea_Y_v", &_par1_emptyarea_Y_v);
    _outtree->Branch("par2_emptyarea_U_v", &_par2_emptyarea_U_v);
    _outtree->Branch("par2_emptyarea_V_v", &_par2_emptyarea_V_v);
    _outtree->Branch("par2_emptyarea_Y_v", &_par2_emptyarea_Y_v);

    _outtree->Branch("par1_pocketarea_U_v", &_par1_pocketarea_U_v);
    _outtree->Branch("par1_pocketarea_V_v", &_par1_pocketarea_V_v);
    _outtree->Branch("par1_pocketarea_Y_v", &_par1_pocketarea_Y_v);
    _outtree->Branch("par2_pocketarea_U_v", &_par2_pocketarea_U_v);
    _outtree->Branch("par2_pocketarea_V_v", &_par2_pocketarea_V_v);
    _outtree->Branch("par2_pocketarea_Y_v", &_par2_pocketarea_Y_v);

    _outtree->Branch("par1_pocketarea_ns_U_v", &_par1_pocketarea_ns_U_v);
    _outtree->Branch("par1_pocketarea_ns_V_v", &_par1_pocketarea_ns_V_v);
    _outtree->Branch("par1_pocketarea_ns_Y_v", &_par1_pocketarea_ns_Y_v);
    _outtree->Branch("par2_pocketarea_ns_U_v", &_par2_pocketarea_ns_U_v);
    _outtree->Branch("par2_pocketarea_ns_V_v", &_par2_pocketarea_ns_V_v);
    _outtree->Branch("par2_pocketarea_ns_Y_v", &_par2_pocketarea_ns_Y_v);

    _outtree->Branch("par1_poly_area_U_v", &_par1_polyarea_U_v);
    _outtree->Branch("par1_poly_area_V_v", &_par1_polyarea_V_v);
    _outtree->Branch("par1_poly_area_Y_v", &_par1_polyarea_Y_v);
    _outtree->Branch("par2_poly_area_U_v", &_par2_polyarea_U_v);
    _outtree->Branch("par2_poly_area_V_v", &_par2_polyarea_V_v);
    _outtree->Branch("par2_poly_area_Y_v", &_par2_polyarea_Y_v);

    _outtree->Branch("par1_polyperimeter_U_v", &_par1_polyperimeter_U_v);
    _outtree->Branch("par1_polyperimeter_V_v", &_par1_polyperimeter_V_v);
    _outtree->Branch("par1_polyperimeter_Y_v", &_par1_polyperimeter_Y_v);
    _outtree->Branch("par2_polyperimeter_U_v", &_par2_polyperimeter_U_v);
    _outtree->Branch("par2_polyperimeter_V_v", &_par2_polyperimeter_V_v);
    _outtree->Branch("par2_polyperimeter_Y_v", &_par2_polyperimeter_Y_v);

    _outtree->Branch("par1_polycharge_U_v", &_par1_polycharge_U_v);
    _outtree->Branch("par1_polycharge_V_v", &_par1_polycharge_V_v);
    _outtree->Branch("par1_polycharge_Y_v", &_par1_polycharge_Y_v);
    _outtree->Branch("par2_polycharge_U_v", &_par2_polycharge_U_v);
    _outtree->Branch("par2_polycharge_V_v", &_par2_polycharge_V_v);
    _outtree->Branch("par2_polycharge_Y_v", &_par2_polycharge_Y_v);

    _outtree->Branch("par1_polyedges_U_v", &_par1_polyedges_U_v);
    _outtree->Branch("par1_polyedges_V_v", &_par1_polyedges_V_v);
    _outtree->Branch("par1_polyedges_Y_v", &_par1_polyedges_Y_v);
    _outtree->Branch("par2_polyedges_U_v", &_par2_polyedges_U_v);
    _outtree->Branch("par2_polyedges_V_v", &_par2_polyedges_V_v);
    _outtree->Branch("par2_polyedges_Y_v", &_par2_polyedges_Y_v);

    _outtree->Branch("par1_polybranches_U_v", &_par1_polybranches_U_v);
    _outtree->Branch("par1_polybranches_V_v", &_par1_polybranches_V_v);
    _outtree->Branch("par1_polybranches_Y_v", &_par1_polybranches_Y_v);
    _outtree->Branch("par2_polybranches_U_v", &_par2_polybranches_U_v);
    _outtree->Branch("par2_polybranches_V_v", &_par2_polybranches_V_v);
    _outtree->Branch("par2_polybranches_Y_v", &_par2_polybranches_Y_v);

    _outtree->Branch("par1_showerfrac_U_v", &_par1_showerfrac_U_v);
    _outtree->Branch("par1_showerfrac_V_v", &_par1_showerfrac_V_v);
    _outtree->Branch("par1_showerfrac_Y_v", &_par1_showerfrac_Y_v);
    _outtree->Branch("par2_showerfrac_U_v", &_par2_showerfrac_U_v);
    _outtree->Branch("par2_showerfrac_V_v", &_par2_showerfrac_V_v);
    _outtree->Branch("par2_showerfrac_Y_v", &_par2_showerfrac_Y_v);
    
    
    //
    // segment information
    //
    _outtree->Branch("par1_electron_frac_U", &_par1_electron_frac_U, "par1_electron_frac_U/F");
    _outtree->Branch("par1_electron_frac_V", &_par1_electron_frac_V, "par1_electron_frac_V/F");
    _outtree->Branch("par1_electron_frac_Y", &_par1_electron_frac_Y, "par1_electron_frac_Y/F");

    _outtree->Branch("par2_electron_frac_U", &_par2_electron_frac_U, "par2_electron_frac_U/F");
    _outtree->Branch("par2_electron_frac_V", &_par2_electron_frac_V, "par2_electron_frac_V/F");
    _outtree->Branch("par2_electron_frac_Y", &_par2_electron_frac_Y, "par2_electron_frac_Y/F");

    _outtree->Branch("par1_muon_frac_U", &_par1_muon_frac_U, "par1_muon_frac_U/F");
    _outtree->Branch("par1_muon_frac_V", &_par1_muon_frac_V, "par1_muon_frac_V/F");
    _outtree->Branch("par1_muon_frac_Y", &_par1_muon_frac_Y, "par1_muon_frac_Y/F");

    _outtree->Branch("par2_muon_frac_U", &_par2_muon_frac_U, "par2_muon_frac_U/F");
    _outtree->Branch("par2_muon_frac_V", &_par2_muon_frac_V, "par2_muon_frac_V/F");
    _outtree->Branch("par2_muon_frac_Y", &_par2_muon_frac_Y, "par2_muon_frac_Y/F");

    _outtree->Branch("par1_proton_frac_U", &_par1_proton_frac_U, "par1_proton_frac_U/F");
    _outtree->Branch("par1_proton_frac_V", &_par1_proton_frac_V, "par1_proton_frac_V/F");
    _outtree->Branch("par1_proton_frac_Y", &_par1_proton_frac_Y, "par1_proton_frac_Y/F");

    _outtree->Branch("par2_proton_frac_U", &_par2_proton_frac_U, "par2_proton_frac_U/F");
    _outtree->Branch("par2_proton_frac_V", &_par2_proton_frac_V, "par2_proton_frac_V/F");
    _outtree->Branch("par2_proton_frac_Y", &_par2_proton_frac_Y, "par2_proton_frac_Y/F");


    return;
  }

  double SelMichelID::Select() {
    LLCV_DEBUG() << "start" << std::endl;
    
    ResetEvent();
      
    float Vertex_X = Data().PGraph()->ParticleArray().front().X();
    float Vertex_Y = Data().PGraph()->ParticleArray().front().Y();
    float Vertex_Z = Data().PGraph()->ParticleArray().front().Z();

    LLCV_DEBUG() << "(RSEV)=("<< Run() << "," << SubRun() << "," << Event() << "," << VertexID() << ")" << std::endl;    
    LLCV_DEBUG() << "VTX=(" 
		 << Vertex_X << "," 
		 << Vertex_Y << "," 
		 << Vertex_Z << ")" << std::endl;

    auto mat_v  = Image().Image<cv::Mat>(kImageADC,_cropx,_cropy);
    auto meta_v = Image().Image<larocv::ImageMeta>(kImageADC,_cropx,_cropy);
    auto img_v  = Image().Image<larcv::Image2D>(kImageADC,_cropx,_cropy);
    auto dead_v = Image().Image<cv::Mat>(kImageDead,_cropx,_cropy);
    auto shr_v  = Image().Image<cv::Mat>(kImageShower,_cropx,_cropy);
    
    for(const auto& meta : meta_v)
      _PixelScan3D.SetPlaneInfo(*meta);

    _white_img = larocv::BlankImage(*(mat_v.front()),0);
    _white_img.setTo(cv::Scalar(0));

    _ContourScan.Reset();
    
    for(const auto& meta : meta_v)
      _ContourScan.SetPlaneInfo(*meta);
    
    std::array<cv::Mat,3> aimg_v;   // adc image
    std::array<cv::Mat,3> timg_v;   // threshold image
    std::array<cv::Mat,3> cimg_v;   // cosmic rejected image
    std::array<cv::Mat,3> dimg_v;   // dead image
    std::array<cv::Mat,3> mat3d_v;  // threshold 8UC3 image

    for(size_t plane=0; plane<3; ++plane) {
      aimg_v[plane]   = *(mat_v[plane]);
      timg_v[plane]   = larocv::Threshold(*(mat_v[plane]),10,255);
      mat3d_v[plane]  = As8UC3(timg_v[plane]);
      dimg_v[plane]   = *(dead_v[plane]);
      
      // tag cosmics
      LLCV_DEBUG() << "Tag cosmic @plane=" << plane << std::endl;
      auto& _CosmicTag = _CosmicTag_v[plane];
      _CosmicTag.Reset();
      _CosmicTag.TagCosmic(timg_v[plane],dimg_v[plane]);

      cimg_v[plane] = timg_v[plane].clone();

      for(const auto& ctor : _CosmicTag.CosmicContours())
	cimg_v[plane] = larocv::MaskImage(cimg_v[plane],ctor,-1,true);  
      
      LLCV_DEBUG() << "Make extension @plane=" << plane << std::endl;
      auto& _LineExtension = _LineExtension_v[plane];
      _LineExtension.SetImageDimension(cimg_v[plane],dimg_v[plane]);
      _LineExtension.SetCosmicPixels(_CosmicTag.CosmicContours());
    }


    larocv::data::Vertex3D vtx3d;
    vtx3d.x = Vertex_X;
    vtx3d.y = Vertex_Y;
    vtx3d.z = Vertex_Z;

    _vertex_x = (float)vtx3d.x;
    _vertex_y = (float)vtx3d.y;
    _vertex_z = (float)vtx3d.z;

    //
    // Project the vertex onto the plane
    //
    larocv::GEO2D_Contour_t vertex_pt_v;
    vertex_pt_v.resize(3);

    for(size_t plane=0; plane<3; ++plane) {
      int px_x = kINVALID_INT;
      int px_y = kINVALID_INT;

      ProjectMat(img_v[plane]->meta(),
		 vtx3d.x,vtx3d.y,vtx3d.z,
    		 px_x, px_y);

      vertex_pt_v[plane] = cv::Point_<int>(px_y,px_x);
    }
    
    //
    // cosmic rejection
    //
    size_t vtx_cosmic_ctr = 0;
    std::array<size_t,3> vtx_cosmic_id_v;
    std::array<int,3> vtx_cosmic_valid_v;

    for(auto& v : vtx_cosmic_valid_v) v = 0;

    for(size_t plane=0; plane<3; ++plane) {
      const auto& vertex_pt = vertex_pt_v[plane];
      const auto& _CosmicTag = _CosmicTag_v[plane];

      auto& edge_n_cosmic = _edge_n_cosmic_v[plane];
      auto& edge_cosmic_vtx_dist = _edge_cosmic_vtx_dist_v[plane];
      auto& edge_cosmic_end_vtx_dist = _edge_cosmic_end_vtx_dist_v[plane];

      edge_n_cosmic = (int)_CosmicTag.CosmicContours().size();

      auto near_id = _CosmicTag.NearestCosmicToPoint(vertex_pt,edge_cosmic_vtx_dist);

      auto& vtx_cosmic_id = vtx_cosmic_id_v[plane];
      auto near_pt = _CosmicTag.NearestCosmicEndToPoint(vertex_pt,
							edge_cosmic_end_vtx_dist,
							vtx_cosmic_id);
      if (edge_cosmic_end_vtx_dist < 20) {
	_edge_cosmic_end_vtx_valid_v[plane] = 1;
	vtx_cosmic_ctr+=1;
      }
      else 
	_edge_cosmic_end_vtx_valid_v[plane] = 0;

    }

    // if there are no cosmics near the vertex on at least 2 planes
    // do not recontruct
    _n_vtx_cosmic = (int)vtx_cosmic_ctr;
    if (vtx_cosmic_ctr < 2) {
      _outtree->Fill();  
      return 0.0;
    }

    // add back in the cosmic pixels for the cosmic 
    // which comes close to the vertex @ cimg
    for(size_t plane=0; plane<3; ++plane) {
       const auto& vtx_cosmic_id = vtx_cosmic_id_v[plane];
       
       if (vtx_cosmic_id == larocv::kINVALID_SIZE) continue;
 
       const auto& _CosmicTag = _CosmicTag_v[plane];
       
       const auto& vtx_cosmic_ctor = _CosmicTag.CosmicContours().at(vtx_cosmic_id);
       
       const auto& timg = timg_v[plane];
       const auto& dimg = dimg_v[plane];
       auto& cimg = cimg_v[plane];

       // get the nonzero pixels inside this contour
       _white_img.setTo(cv::Scalar(255));
       auto vtx_cosmic_ctor_pix_v = larocv::FindNonZero(larocv::MaskImage(_white_img,vtx_cosmic_ctor,-1,false));
       
       for(const auto& vtx_cosmic_ctor_pix : vtx_cosmic_ctor_pix_v) {
	 // check if in dead region, if so, add a 10
	 auto row = vtx_cosmic_ctor_pix.y;
	 auto col = vtx_cosmic_ctor_pix.x;
	 
	 int dead_px = (int)(dimg.at<uchar>(row,col));

	 if (dead_px == 0) {
	   cimg.at<uchar>(row,col) = (uchar)255;
	   continue;
	 }
	 
	 cimg.at<uchar>(row,col) = (uchar)(timg.at<uchar>(row,col));
       } // end pixel point

    } // end inserting cosmic track @ plane

    
    //
    // Generate initial 2D clusters
    //
    std::array<larocv::GEO2D_ContourArray_t,3> plane_ctor_vv;
    std::array<larocv::GEO2D_Contour_t,3> vertex_ctor_v;
    std::array<size_t,3> close_id_v;
    
    std::array<larocv::GEO2D_ContourArray_t,3> par_ctor_v;
    std::array<larocv::GEO2D_ContourArray_t,3> line_contour_vv;
    std::array<std::vector<Triangle>,3> triangle_vv;
    std::array<geo2d::VectorArray<float>,3> edge_vv;
    
    float _distance_tol = 20;

    for(size_t plane=0; plane<3; ++plane) {
      
      // skip the plane which does not have a incoming muon
      if (_edge_cosmic_end_vtx_valid_v[plane] == 0)
	continue;
      
      auto& cimg         = cimg_v[plane];
      auto& close_id     = close_id_v[plane];
      auto& vertex_ctor  = vertex_ctor_v[plane];
      auto& plane_ctor_v = plane_ctor_vv[plane];
      
      plane_ctor_v = larocv::FindContours(cimg);
      
      LLCV_DEBUG() << "@plane=" << plane 
		   << " found " << plane_ctor_v.size() 
		   << " ctors" << std::endl;
      
      const auto vertex_pt = vertex_pt_v[plane];
      
      float distance = larocv::kINVALID_FLOAT;
      close_id = FindClosestContour(plane_ctor_v,vertex_pt,distance);

      if (close_id == larocv::kINVALID_SIZE) {
	LLCV_WARNING() << "No contour? Skip." << std::endl;
	continue;
      }

      if (distance > _distance_tol) {
	LLCV_WARNING() << "Distance=" << distance << ">" << _distance_tol << ". Skip." << std::endl;
	continue;
      }

      vertex_ctor = plane_ctor_v[close_id];
      
      // mask all but vertex_ctor
      auto cimg_vertex_mask = larocv::MaskImage(cimg,vertex_ctor,-1,false);

      // minimize over circle of radius 4
      auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(cimg_vertex_mask,geo2d::Circle<float>(vertex_pt,4),-1,false));

      float nratio_max = -1.0*larocv::kINVALID_FLOAT;

      for(const auto& nz_pt : nz_pt_v) {
	
	// mask out the vertex from image
	auto cimg_mask = larocv::MaskImage(cimg_vertex_mask,geo2d::Circle<float>(nz_pt,5),-1,true);
      
	// find contours
	auto& par_ctor = par_ctor_v[plane];
	par_ctor = larocv::FindContours(cimg_mask);

	auto nctor = par_ctor.size();

	// number particles == 1 veto
	if (nctor < 2) continue;
	
	larocv::GEO2D_ContourArray_t local_line_contour_v;
	std::vector<Triangle>        local_triangle_v;
	geo2d::VectorArray<float>    local_edge_v;
	std::vector<float>           local_line_frac_v;

	local_line_contour_v.reserve(nctor);
	local_triangle_v.reserve(nctor);
	local_edge_v.reserve(nctor);
	local_line_frac_v.reserve(nctor);

	for(size_t pid=0; pid<(size_t)nctor; ++pid) {

	  auto par = par_ctor[pid];

	  // estimate the end point
	  Triangle local_triangle(par,nz_pt);
	  
	  geo2d::Vector<float> local_edge;
	  float nline_pixels = 0;
	  float npar_pixels  = 0;
	  auto local_line_ctor = MaximizeTriangleLine(cimg,
						      local_triangle,
						      nline_pixels,
						      npar_pixels,
						      local_edge,
						      _white_img);
	  if (local_line_ctor.empty()) continue;
	  
	  float nratio = 0;
	  if (npar_pixels > 0)
	    nratio = nline_pixels / npar_pixels;
	  else
	    continue;
	  
	  local_line_frac_v.emplace_back(nratio);
	  local_line_contour_v.emplace_back(local_line_ctor);
	  local_triangle_v.emplace_back(local_triangle);
	  local_edge_v.emplace_back(local_edge);
	  
	} // end contour

	// number particles == 1 veto
	if (local_line_frac_v.size() < 2) 
	  continue;

	// check ratio product
	float ratio_prod = 1;
	for(auto local_line_frac : local_line_frac_v)
	  ratio_prod *= local_line_frac;
	
	// if it's higher, save
	if (ratio_prod > nratio_max) {
	  line_contour_vv[plane] = local_line_contour_v;
	  triangle_vv[plane]     = local_triangle_v;
	  edge_vv[plane]         = local_edge_v;
	  
	  nratio_max = ratio_prod;
	  LLCV_DEBUG() << "saved r=" << ratio_prod << std::endl;
	}

      } // end nz_pt

    } // end plane
    
    //
    // extend the lines across dead regions, keep going until we run it out
    // make object2d
    //
    LLCV_DEBUG() << "Object2D creation..." << std::endl;
    std::array<std::vector<Object2D>,3> object_vv;
    for(size_t plane=0; plane<3; ++plane) {
      LLCV_DEBUG() << "@plane=" << plane << std::endl;

      const auto& _LineExtension = _LineExtension_v[plane];
      const auto& plane_ctor_v   = plane_ctor_vv[plane];
      const auto& cimg           = cimg_v[plane];
      const auto& close_id       = close_id_v[plane];

      auto& triangle_v = triangle_vv[plane];
      auto& edge_v     = edge_vv[plane];
      auto& object_v   = object_vv[plane];

      object_v.resize(triangle_v.size());

      std::vector<bool> used_v;
      for(size_t lid=0; lid<triangle_v.size(); ++lid) {

	used_v.clear();
	used_v.resize(plane_ctor_v.size(),false);
	used_v[close_id] = true;

	auto& triangle = triangle_v[lid];
	auto& start    = triangle.Apex();
	auto& edge     = edge_v[lid];
	auto& object   = object_v[lid];

	object._polygon_v.emplace_back(triangle.Contour(),start);

	auto old_edge = edge;
	auto new_edge = _LineExtension.ExtendAcross(start,edge);
	
	LLCV_DEBUG() << "old_edge=(" << old_edge.x << "," << old_edge.y << ")" << std::endl;
	LLCV_DEBUG() << "new_edge=(" << new_edge.x << "," << new_edge.y << ")" << std::endl;

	geo2d::Vector<float> inv_pt(0,0);

	size_t ix=0;
	while ((new_edge.x != old_edge.x) and (new_edge.y != old_edge.y)) {

	  LLCV_DEBUG() << "@ix=" << ix << std::endl;
	  _white_img.setTo(cv::Scalar(0));
	  cv::line(_white_img,start,new_edge,cv::Scalar(255),3);
	  auto lc_v = larocv::FindContours(_white_img);
	  if (lc_v.empty()) continue;
	  const auto& lc = lc_v.front();
	  
	  for(size_t ict=0; ict<plane_ctor_v.size(); ++ict) {
	    if (used_v[ict]) continue;
	    if (!larocv::AreaOverlap(plane_ctor_v[ict],lc)) continue;
	    object._polygon_v.emplace_back(plane_ctor_v[ict],start);
	    used_v[ict] = true;
	  }

	  Triangle new_triangle(start,inv_pt,inv_pt);
	  float nline_pixels = 0;
	  float npar_pixels  = 0;
	  auto line_ctor = MaximizePolygonLine(cimg,
					       object._polygon_v,
					       new_triangle,
					       nline_pixels,
					       npar_pixels,
					       new_edge,
					       _white_img);

	  if ((new_edge.x == old_edge.x) and (new_edge.y == old_edge.y)) {
	    LLCV_DEBUG() << "circular logic, break" << std::endl;
	    break; 
	  }

	  old_edge = new_edge;
	  new_edge = _LineExtension.ExtendAcross(start,old_edge);

	  LLCV_DEBUG() << "old_edge=(" << old_edge.x << "," << old_edge.y << ")" << std::endl;
	  LLCV_DEBUG() << "new_edge=(" << new_edge.x << "," << new_edge.y << ")" << std::endl;
	  ix+=1;
	}

	Triangle new_triangle(start,inv_pt,inv_pt);
	float nline_pixels = 0;
	float npar_pixels  = 0;
	auto line_ctor = MaximizePolygonLine(cimg,
					     object._polygon_v,
					     new_triangle,
					     nline_pixels,
					     npar_pixels,
					     new_edge,
					     _white_img);

	// make sure the ege is _not_ outside the image
	if (new_edge.x >= (float)((_cropx)-0.5)) {
	  LLCV_WARNING() << "new_edge.x=" << new_edge.x << " outside cropx" << std::endl;
	  new_edge.x = _cropx - 1;
	}
	if (new_edge.y >= (float)((_cropy)-0.5))  {
	  LLCV_WARNING() << "new_edge.y=" << new_edge.y << " outside cropy" << std::endl;
	  new_edge.y = _cropy - 1;
	}

	float nratio = 0;
	if (npar_pixels > 0)
	  nratio = nline_pixels / npar_pixels;

	object._triangle  = new_triangle;
	object._line      = line_ctor;
	object._line_frac = nratio;
	object._edge      = new_edge;
	object._plane     = plane;
	
	LLCV_DEBUG() << "@lid=" << lid 
		     << " @start=[" << start.x << "," << start.y << "]"
		     << "& @edge=[" << new_edge.x << "," << new_edge.y << "]" << std::endl;
      } // end this triangle
    } // end plane

    LLCV_DEBUG() << "...done" << std::endl;

    _Match.ClearEvent();
    _Match.ClearMatch();

    LLCV_DEBUG() << "Match lines" << std::endl;
    for(size_t plane=0; plane<3; ++plane) {
      LLCV_DEBUG() << "@plane=" << plane << std::endl;
      for(const auto& object : object_vv[plane]) {
	_Match.Register(object,plane);
      }
      LLCV_DEBUG() << "num register=" << object_vv[plane].size() << std::endl;
    }
    
    std::vector<float> score_v;
    auto match_vv = _Match.MatchObjects(score_v);

    LLCV_DEBUG() << "& recieved " << match_vv.size() << " matched particles" << std::endl;
    
    std::vector<Object2DCollection> obj_col_v;
    obj_col_v.resize(match_vv.size());

    for(size_t mid=0; mid< match_vv.size(); ++mid) {
      auto match_v = match_vv[mid];
      auto score = score_v[mid];

      auto& obj_col = obj_col_v[mid];
      
      obj_col.SetStart(Vertex_X,
		       Vertex_Y,
		       Vertex_Z);
      
      obj_col.SetScore(score);
      
      // Fill the match
      for (auto match : match_v) {
	auto plane = match.first;
	auto id    = match.second;
	LLCV_DEBUG() << "@plane=" << plane << " id=" << id << std::endl;
	obj_col.emplace_back(object_vv[plane][id]);
      }

      LLCV_DEBUG() << "@mid=" << mid << std::endl;
     
      _ShowerTools.ReconstructAngle(img_v,aimg_v,obj_col);
      _ShowerTools.ReconstructLength(img_v,aimg_v,obj_col);
      _ShowerTools.ReconstructdQdx(img_v,aimg_v,obj_col,3+3); //offset by ~2cm for 5 pixel vertex mask out

      LLCV_DEBUG() << "theta=" << obj_col.Theta() << " phi=" << obj_col.Phi() << std::endl;
      LLCV_DEBUG() << "(" << obj_col.dX() << "," << obj_col.dY() << "," << obj_col.dZ() << ")" << std::endl;
      LLCV_DEBUG() << "length=" << obj_col.Length() << std::endl;
      LLCV_DEBUG() << "score=" << score << std::endl;

      auto pca_v = EstimateDirection(obj_col,_white_img,_ContourScan,_ShowerTools);

      LLCV_DEBUG() << "pca=(" << pca_v[0] << "," << pca_v[1] << "," << pca_v[2] << ")" << std::endl;
     
      obj_col.SetddX(pca_v[0]);
      obj_col.SetddY(pca_v[1]);
      obj_col.SetddZ(pca_v[2]);
      
      // determine near vertex parameters
      for(auto& obj2d : obj_col) {
	const auto plane = obj2d.Plane();
	_white_img.setTo(cv::Scalar(1));
	obj2d.LineVertex(*(img_v[plane]),aimg_v[plane],_white_img,5);
      }

      _ShowerTools.ReconstructdQdxProfile(img_v,aimg_v,obj_col);
      _ShowerTools.TruncatedQdxProfile(obj_col,_tsigma);

      // determine the end point
      for(const auto& obj2d : obj_col) 
	_ContourScan.RegisterEndPoint(obj2d.Edge(),obj2d.Plane());
      
      auto end_pt_v = _ContourScan.EndPoint();
      obj_col.SetEndX(end_pt_v[0]);
      obj_col.SetEndY(end_pt_v[1]);
      obj_col.SetEndZ(end_pt_v[2]);

      // determine line fraction parameters
      ComputeLineParameters(obj_col,cimg_v,_white_img);
      SplitLineParameters(obj_col,cimg_v,_white_img);
    }

    
    //
    // Detect brem function
    //
    LLCV_DEBUG() << "Detecting brem" << std::endl;
    // mask out all polygons from all particles, find contours in cimg
    std::array<larocv::GEO2D_ContourArray_t,3> cimg_ctor_vv;
    for(size_t plane=0; plane<3; ++plane) {
      auto cimg_mask = cimg_v[plane].clone();
      auto& cimg_ctor_v = cimg_ctor_vv[plane];
      
      for(size_t oid=0; oid<obj_col_v.size(); ++oid) {
	auto& obj_col = obj_col_v[oid];
	if (!obj_col.HasObject(plane)) continue;
	const auto& obj2d = obj_col.PlaneObject(plane);
	cimg_mask = larocv::MaskImage(cimg_mask,geo2d::Circle<float>(obj2d.triangle().Apex(),5),-1,true);
	for(size_t polyid=0; polyid<obj2d.Polygons().size(); ++polyid) {
	  const auto& polygon = obj2d.Polygons()[polyid];
	  cimg_mask = larocv::MaskImage(cimg_mask,polygon.Contour(),-1,true);
	} // end polygon
      } // end object collection

      cimg_ctor_v = larocv::FindContours(cimg_mask);
    } // end plane
    
    
    for(size_t oid=0; oid<obj_col_v.size(); ++oid) {
      auto& obj_col = obj_col_v[oid];

      for(size_t plane=0; plane<3; ++plane) {
	if (!obj_col.HasObject(plane)) continue;
	auto& obj2d = obj_col.PlaneObjectRW(plane);
	
	const auto& cimg_ctor_v = cimg_ctor_vv[plane];

	// orient a triangle in direction of line
	auto edge = obj2d.Edge();

	obj2d._brem_triangle = obj2d.triangle().RotateToPoint(edge,2.0);
	obj2d._brem_triangle.SetContour(obj2d.triangle().Contour());
	
	// look for brem
	std::vector<size_t> brem_v;
	std::vector<size_t> inside_v;
	auto nbrem = DetectBrem(obj2d._brem_triangle,
				cimg_ctor_v,
				brem_v,
				inside_v,
				_brem_dist,
				_brem_size);
	
	obj2d._n_brem = nbrem;

	// store all the contours inside the expanded triangle
	obj2d._expand_polygon_v = obj2d._polygon_v;

	for(auto iid : inside_v)
	  obj2d._expand_polygon_v.emplace_back(cimg_ctor_v[iid],obj2d.Start());
	
	if (!brem_v.empty())
	  obj2d._brem_index = obj2d._polygon_v.size();
	else
	  obj2d._brem_index = -1;

	// store the brem only into polygon
	for(auto bid : brem_v)
	  obj2d._polygon_v.emplace_back(cimg_ctor_v[bid],obj2d.Start());

      }
    }
    LLCV_DEBUG() << "...done" << std::endl;


    //
    // Detect branches for obj2d extended polygons
    //
    LLCV_DEBUG() << "Detecting branches..." << std::endl;
    for(auto& obj_col : obj_col_v) {
      for(size_t plane=0; plane<3; ++plane) {
	if (!obj_col.HasObject(plane)) continue;
	const auto& cimg = cimg_v[plane];
	auto& obj2d = obj_col.PlaneObjectRW(plane);
	for(auto& polygon : obj2d._expand_polygon_v) 
	  polygon.DetectBranching(cimg,4,2,5,5);
      }
    }
    LLCV_DEBUG() << "...done" << std::endl;

    //
    // distance of vertex contour to edge
    //
    for(size_t plane=0; plane<3; ++plane) {
      const auto& vertex_ctor = vertex_ctor_v[plane];
      auto& edge_dist = _edge_dist_v[plane];
      edge_dist = MinimizeToEdge(vertex_ctor,_cropx,_cropy);
    }

    //
    // vertex charge + dead region check
    //
    for(size_t plane=0; plane<3; ++plane) {
      // mask out the vertex from image
      const auto& vertex_pt = vertex_pt_v[plane];
      geo2d::Circle<float> vtx_circle(vertex_pt.x,vertex_pt.y,5);

      _white_img.setTo(cv::Scalar(255));
      auto cimg_mask  = larocv::MaskImage(cimg_v[plane],vtx_circle,-1,false);
      auto dead_mask  = larocv::MaskImage(dimg_v[plane],vtx_circle,-1,false);
      auto white_mask = larocv::MaskImage(_white_img,vtx_circle,-1,false);

      int cimg_count  = (int)larocv::CountNonZero(cimg_mask);
      int dead_count  = (int)larocv::CountNonZero(dead_mask);
      int white_count = (int)larocv::CountNonZero(white_mask);

      int is_charge = 0;
      int near_dead = 0;

      if (cimg_count > 0) is_charge = 1;
      if (white_count != dead_count) near_dead = 1;

      if (plane==0) {
	_vtx_charge_U = is_charge;
	_vtx_dead_U   = near_dead;
      }
      if (plane==1) {
	_vtx_charge_V = is_charge;
	_vtx_dead_V   = near_dead;
      }
      if (plane==2) {
	_vtx_charge_Y = is_charge;
	_vtx_dead_Y   = near_dead;
      }
    }


    // number of matched particles
    _n_par = (int)obj_col_v.size();
    _par_score_v.clear();
    _par_score_v.resize(_n_par,-1*larocv::kINVALID_FLOAT);
    for(size_t oid=0; oid<obj_col_v.size(); ++oid)
      _par_score_v[oid] = obj_col_v[oid].Score();

    for(size_t plane=0; plane<3; ++plane) {
      for(const auto& object : object_vv[plane]) {
	if (plane==0) {
	  _vtx_xing_U = (int)object_vv[plane].size();
	  _vtx_linelen_U_v.push_back(object.LineLength());
	  _vtx_linefrac_U_v.push_back(object.LineFrac());
	}
	if (plane==1) {
	  _vtx_xing_V = (int)object_vv[plane].size();
	  _vtx_linelen_V_v.push_back(object.LineLength());
	  _vtx_linefrac_V_v.push_back(object.LineFrac());
	}
	if (plane==2) {
	  _vtx_xing_Y = (int)object_vv[plane].size();
	  _vtx_linelen_Y_v.push_back(object.LineLength());
	  _vtx_linefrac_Y_v.push_back(object.LineFrac());
	}
      }
    }




    // object collection stuff
    auto& out_pgraph = Data().MakePGraph();

    for(size_t oid=0; oid<obj_col_v.size(); ++oid) {
      const auto& obj_col = obj_col_v[oid];
      
      LLCV_DEBUG() << "@oid=" << oid << std::endl;

      // note to future peoples: it's OK to skip particles > 2 because
      // the matching algorithm returns particles in descending score value

      if (oid > 1) continue;

      SetParticle(oid);
      
      *_par_theta    = obj_col.Theta();
      *_par_phi      = obj_col.Phi();
      *_par_length   = obj_col.Length();
      *_par_score    = obj_col.Score();
      *_par_dx1      = obj_col.dX();
      *_par_dy1      = obj_col.dY();
      *_par_dz1      = obj_col.dZ();
      *_par_dx2      = obj_col.ddX();
      *_par_dy2      = obj_col.ddY();
      *_par_dz2      = obj_col.ddZ();
      *_par_nplanes  = (int)obj_col.size();
      *_par_planes_v = obj_col.Planes();

      _white_img.setTo(cv::Scalar(1));
      *_par_xdead_v  = obj_col.XDead(dimg_v,_white_img);
      *_par_end_pt_v = obj_col.EndPoint();

      for(size_t plane=0; plane<3; ++plane) {
	LLCV_DEBUG() << "@plane=" << plane << std::endl;
	if (!obj_col.HasObject(plane)) continue;
	LLCV_DEBUG() << "...exists" << std::endl;

	const auto& obj2d = obj_col.PlaneObject(plane);

	size_t npolygons = obj2d.NPolygons();

	SetParticlePlane(oid,plane);
	
	*_par_n_polygons     = obj2d.NPolygons();
	*_par_linelength     = obj2d.LineLength();
	*_par_linefrac       = obj2d.LineFrac();   

	_white_img.setTo(cv::Scalar(255));
	*_par_linefrac_empty = obj2d.LineFracEmpty(_white_img);
	
	*_par_linedx         = obj2d.LinedX();
	*_par_linedy         = obj2d.LinedY();

	*_par_line_vtx_density  = obj2d.LineVertexDensity();
	*_par_line_vtx_coverage = obj2d.LineVertexCoverage();
	*_par_line_vtx_charge   = obj2d.LineVertexCharge();
	*_par_line_mean_dist    = obj2d.LineMeanDist();
	*_par_line_max_dist     = obj2d.LineMaxDist();

	*_par_line_first_half_linefrac  = obj2d.LineFirstHalfLineFrac();
	*_par_line_second_half_linefrac = obj2d.LineSecondHalfLineFrac();

	*_par_triangle_height         = obj2d.triangle().Height();
	*_par_triangle_mid_to_edge    = geo2d::dist(obj2d.triangle().MidPoint(),obj2d.Edge());
	*_par_triangle_emptyarearatio = obj2d.triangle().EmptyAreaRatio();
	*_par_triangle_emptyarea      = obj2d.triangle().EmptyArea();
	*_par_triangle_baselength     = obj2d.triangle().BaseLength();
	*_par_triangle_area           = obj2d.triangle().Area();
	*_par_triangle_brem           = obj2d.NBrem();
	*_par_triangle_coverage       = obj2d.triangle().Coverage(aimg_v[plane]);
	*_par_brem_triangle_coverage  = obj2d.brem_triangle().Coverage(aimg_v[plane]);
	*_par_expand_charge           = obj2d.Charge(*(img_v[plane]),aimg_v[plane]);
	*_par_dqdx                    = obj2d.dQdx();
	*_par_dqdx_v                  = obj2d.dQdxProfile();
	*_par_tdqdx_v                 = obj2d.TdQdxProfile();
	*_par_dqdx_step               = obj2d.dQdxStep(); 
	*_par_dqdx_pitch              = obj2d.dQdxPitch(); 
	*_par_length3d                = obj2d.Length();
	*_par_brem_idx                = obj2d.BremIndex();
	*_par_showerfrac              = obj2d.Fraction(*shr_v[plane],aimg_v[plane]);
	
	(*_par_cosmic_dist_v)[plane]         = NearestPolygonToCosmic(obj2d.Polygons(),_CosmicTag_v,plane);
	(*_par_cosmic_dist_end_v)[plane]     = PointCosmicDistance(obj2d.Edge(),_CosmicTag_v,plane);
	(*_par_cosmic_end_dist_v)[plane]     = NearestPolygonToCosmicEnd(obj2d.Polygons(),_CosmicTag_v,plane);
	(*_par_cosmic_end_dist_end_v)[plane] = PointCosmicEndDistance(obj2d.Edge(),_CosmicTag_v,plane);

	auto nexpanded_polygons = obj2d.ExpandedPolygons().size();
	ResizePlanePolygon(nexpanded_polygons);

	for(size_t polyid=0; polyid<nexpanded_polygons; ++polyid) {
	  const auto& polygon = obj2d.ExpandedPolygons()[polyid];
	  (*_par_numberdefects_v)[polyid]     = polygon.NumberDefects(5);
	  (*_par_numberdefects_ns_v)[polyid]  = polygon.NumberDefectsNoStart(5);
	  (*_par_largestdefect_v)[polyid]     = polygon.LargestDefect();
	  (*_par_smallestdefect_v)[polyid]    = polygon.SmallestDefect();
	  (*_par_largestdefect_ns_v)[polyid]  = polygon.LargestDefectNoStart();
	  (*_par_smallestdefect_ns_v)[polyid] = polygon.SmallestDefectNoStart();
	  (*_par_emptyarearatio_v)[polyid]    = polygon.EmptyAreaRatio();
	  (*_par_emptyarea_v)[polyid]         = polygon.EmptyArea();
	  (*_par_pocketarea_v)[polyid]        = polygon.PocketArea();
	  (*_par_pocketarea_ns_v)[polyid]     = polygon.PocketAreaNoStart();
	  (*_par_polyarea_v)[polyid]          = polygon.Area();
	  (*_par_polyperimeter_v)[polyid]     = polygon.Perimeter();
	  (*_par_polycharge_v)[polyid]        = polygon.Charge(*(img_v[plane]),aimg_v[plane]);
	  (*_par_polyedges_v)[polyid]         = (int)polygon.Edges().size();
	  (*_par_polybranches_v)[polyid]      = (int)polygon.Branches().size();
	  (*_par_showerfrac_v)[polyid]        = (float)polygon.Fraction(*shr_v[plane],aimg_v[plane]);

	} // end polygon on plane

      } // end plane

      //
      // fill a ROI per particle
      //
      larcv::ROI proi;
      proi.Position(Vertex_X,Vertex_Y,Vertex_Z,larocv::kINVALID_DOUBLE);
      
      for(size_t plane=0; plane<3; ++plane) 
	proi.AppendBB((*(img_v.at(plane))).meta());
      
      out_pgraph.Emplace(std::move(proi),oid);
      
      for(size_t plane=0; plane<3; ++plane) {
	std::vector<larcv::Pixel2D> pixel_v;
	
	if (obj_col.HasObject(plane)) {
	  const auto& obj2d = obj_col.PlaneObject(plane);
	  
	  for( const auto& polygon : obj2d.ExpandedPolygons()) {
	    _white_img.setTo(cv::Scalar(255));
	    auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(_white_img,polygon.Contour(),-1,false));
	    for(const auto& pt : nz_pt_v) {
	      auto pt_img2d = MatToImage2D(pt,*(mat_v.at(plane)));
	      
	      auto row = pt_img2d.x;
	      auto col = pt_img2d.y;
	      
	      pixel_v.emplace_back(row,col);
	      pixel_v.back().Intensity((*(img_v.at(plane))).pixel(row,col));
	      
	    } // cv::Mat pt
	  
	  } // end polygon
	} // no plane object
	
	LLCV_DEBUG() << "pixel2d sz=" << pixel_v.size() << std::endl;
	
	auto& out_pixel_cluster = Data().MakePixel2DCluster(plane,(*(img_v.at(plane))).meta());
	out_pixel_cluster = larcv::Pixel2DCluster(std::move(pixel_v));
	
      } // end plane
      
      
      //
      // fill a larlite track
      //
      auto& out_track = Data().MakeTrack();
      FillTrack(obj_col, out_track);
    } // end particle

    //
    // write out the cosmic removed image
    //
    for(size_t plane=0; plane<3; ++plane) {
      std::vector<larcv::Pixel2D> pixel_v;
      auto nz_pt_v = larocv::FindNonZero(cimg_v[plane]);
      for(const auto& nz_pt : nz_pt_v) {
	auto pt_img2d = MatToImage2D(nz_pt,*(mat_v.at(plane)));
	      
	auto row = pt_img2d.x;
	auto col = pt_img2d.y;
	      
	pixel_v.emplace_back(row,col);
	pixel_v.back().Intensity((*(img_v.at(plane))).pixel(row,col));
      }
      auto& out_img = Data().MakeImage(plane,(*(img_v[plane])).meta());
      out_img = larcv::Pixel2DCluster(std::move(pixel_v));
    }
    

    //
    // write out the interaction image
    //
    for(size_t plane=0; plane<3; ++plane) {

      std::vector<larcv::Pixel2D> pixel_v;
      geo2d::Vector<float> start_pt;
      
      // insert the particle clusters
      for(size_t oid=0; oid<obj_col_v.size(); ++oid) {
	const auto& obj_col = obj_col_v[oid];
	if (!obj_col.HasObject(plane)) continue;
	const auto& obj2d = obj_col.PlaneObject(plane);
	
	start_pt = obj2d.Start();

	for( const auto& polygon : obj2d.ExpandedPolygons()) {
	    _white_img.setTo(cv::Scalar(255));
	    auto nz_pt_v = larocv::FindNonZero(larocv::MaskImage(_white_img,polygon.Contour(),-1,false));
	    for(const auto& pt : nz_pt_v) {
	      auto pt_img2d = MatToImage2D(pt,*(mat_v.at(plane)));
	      
	      auto row = pt_img2d.x;
	      auto col = pt_img2d.y;
	      
	      pixel_v.emplace_back(row,col);
	      pixel_v.back().Intensity((*(img_v.at(plane))).pixel(row,col));
	    }
	  
	  }
	
      }

      // insert the vertex region
      _white_img.setTo(cv::Scalar(255));

      auto cimg_start = larocv::MaskImage(*(mat_v.at(plane)),geo2d::Circle<float>(start_pt,5),-1,false);
      auto nz_pt_v = larocv::FindNonZero(cimg_start);
      
      for(const auto& pt : nz_pt_v) {
	auto pt_img2d = MatToImage2D(pt,*(mat_v.at(plane)));
	
	auto row = pt_img2d.x;
	auto col = pt_img2d.y;
	      
	pixel_v.emplace_back(row,col);
	pixel_v.back().Intensity((*(img_v.at(plane))).pixel(row,col));
      }
      
      auto& out_inter = Data().MakeInteraction(plane,(*(img_v[plane])).meta());
      out_inter = larcv::Pixel2DCluster(std::move(pixel_v));
    }

    _outtree->Fill();
    
    //
    // Debug print out
    //
    if (_debug) {

      LLCV_DEBUG() << "DEBUG" << std::endl;

      for(size_t plane=0; plane<3; ++plane) {
	
	auto write_img = cimg_v[plane].clone();

	auto mat3d = As8UC3(write_img);

	const auto& _CosmicTag = _CosmicTag_v[plane];

	for(size_t oid=0; oid<obj_col_v.size(); ++oid) {
	  const auto& obj_col = obj_col_v[oid];
	  if (!obj_col.HasObject(plane)) continue;
	  const auto& obj2d = obj_col.PlaneObject(plane);
	  
	  for(const auto& polygon : obj2d.Polygons()) {
	    for(const auto& branch : polygon.Branches()) 
	      mat3d.at<cv::Vec3b>(branch.y,branch.x) = {255,51,255};
	    for(const auto& edge : polygon.Edges()) 
	      mat3d.at<cv::Vec3b>(edge.y,edge.x) = {0,128,255};
	  }

	  cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,obj2d.Line()),-1,cv::Scalar(255,255,0));
	  cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,obj2d.triangle().AsContour()),-1,cv::Scalar(0,255,255));
	  
	  const auto & tri_brem = obj2d.brem_triangle();
	  cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,tri_brem.AsContour()),-1,cv::Scalar(255,176,102));
	  
	  for(const auto& poly : obj2d.ExpandedPolygons()) 
	    cv::drawContours(mat3d,larocv::GEO2D_ContourArray_t(1,poly.Hull()),-1,cv::Scalar(0,0,255));

	  geo2d::Line<float> line(obj2d.Start(),obj2d.Dir());
	  
	  if (obj2d.Dir().x!=0) {
	    for(int i = -150; i < 150; ++i) {
	      int pt_i_x = obj2d.Start().x+i;
	      int pt_i_y = line.y(pt_i_x);
	      
	      if (pt_i_x>=400) continue;
	      if (pt_i_y>=400) continue;
	      
	      if (pt_i_x<0) continue;
	      if (pt_i_y<0) continue;
	      
	      mat3d.at<cv::Vec3b>(pt_i_y,pt_i_x) = {0,165,255};
	    }
	  }
	}
	
      
	// _CosmicTag.DrawLines(mat3d);
	_CosmicTag.DrawContours(mat3d);
	
	std::stringstream ss;
	ss.str("");
	ss << "cpng/plane_img_" << Run() << "_" << SubRun() << "_" << Event() << "_" << VertexID() << "_" << plane << ".png";
	cv::imwrite(ss.str(),mat3d);
      }
    }
    
    LLCV_DEBUG() << "end" << std::endl;

    return 0.0;
  }


  void SelMichelID::ResetEvent() {
    
    _n_par        = -1.0*larocv::kINVALID_INT;
    _n_vtx_cosmic = -1.0*larocv::kINVALID_INT;

    _edge_dist_v.clear();
    _edge_dist_v.resize(3,-1);

    _edge_n_cosmic_v.clear();
    _edge_n_cosmic_v.resize(3,-1);

    _edge_cosmic_vtx_dist_v.clear();
    _edge_cosmic_vtx_dist_v.resize(3,-1);

    _edge_cosmic_end_vtx_dist_v.clear();
    _edge_cosmic_end_vtx_dist_v.resize(3,-1);

    _edge_cosmic_end_vtx_valid_v.clear();
    _edge_cosmic_end_vtx_valid_v.resize(3,-1);

    _vtx_xing_U = -1.0*larocv::kINVALID_INT;
    _vtx_xing_V = -1.0*larocv::kINVALID_INT;
    _vtx_xing_Y = -1.0*larocv::kINVALID_INT;

    _vtx_charge_U = -1.0*larocv::kINVALID_INT;
    _vtx_charge_V = -1.0*larocv::kINVALID_INT;
    _vtx_charge_Y = -1.0*larocv::kINVALID_INT;

    _vtx_dead_U = -1.0*larocv::kINVALID_INT;
    _vtx_dead_V = -1.0*larocv::kINVALID_INT;
    _vtx_dead_Y = -1.0*larocv::kINVALID_INT;

    _vtx_linelen_U_v.clear();
    _vtx_linelen_V_v.clear();
    _vtx_linelen_Y_v.clear();

    _vtx_linefrac_U_v.clear();
    _vtx_linefrac_V_v.clear();
    _vtx_linefrac_Y_v.clear();

    //
    // 3D information
    //
    _par1_theta    = -1.0*larocv::kINVALID_FLOAT;
    _par1_phi      = -1.0*larocv::kINVALID_FLOAT;
    _par1_length   = -1.0*larocv::kINVALID_FLOAT;
    _par1_score    = -1.0*larocv::kINVALID_FLOAT;
    _par1_dx1      = -1.0*larocv::kINVALID_FLOAT;
    _par1_dy1      = -1.0*larocv::kINVALID_FLOAT;
    _par1_dz1      = -1.0*larocv::kINVALID_FLOAT;
    _par1_dx2      = -1.0*larocv::kINVALID_FLOAT;
    _par1_dy2      = -1.0*larocv::kINVALID_FLOAT;
    _par1_dz2      = -1.0*larocv::kINVALID_FLOAT;
    _par1_nplanes  = -1*larocv::kINVALID_INT;

    _par1_planes_v.clear();
    _par1_planes_v.resize(3,-1);

    _par1_xdead_v.clear();
    _par1_xdead_v.resize(3,-1);

    _par1_cosmic_dist_v.clear();
    _par1_cosmic_dist_v.resize(3,-1);    

    _par1_cosmic_dist_end_v.clear();
    _par1_cosmic_dist_end_v.resize(3,-1);    

    _par1_cosmic_end_dist_v.clear();
    _par1_cosmic_end_dist_v.resize(3,-1);    

    _par1_cosmic_end_dist_end_v.clear();
    _par1_cosmic_end_dist_end_v.resize(3,-1);    

    _par1_end_pt_v.clear();
    _par1_end_pt_v.resize(3,-1);

    _par2_theta   = -1.0*larocv::kINVALID_FLOAT;
    _par2_phi     = -1.0*larocv::kINVALID_FLOAT;
    _par2_length  = -1.0*larocv::kINVALID_FLOAT;
    _par2_score   = -1.0*larocv::kINVALID_FLOAT;
    _par2_dx1     = -1.0*larocv::kINVALID_FLOAT;
    _par2_dy1     = -1.0*larocv::kINVALID_FLOAT;
    _par2_dz1     = -1.0*larocv::kINVALID_FLOAT;
    _par2_dx2     = -1.0*larocv::kINVALID_FLOAT;
    _par2_dy2     = -1.0*larocv::kINVALID_FLOAT;
    _par2_dz2     = -1.0*larocv::kINVALID_FLOAT;
    _par2_nplanes = -1*larocv::kINVALID_INT;

    _par2_planes_v.clear();
    _par2_planes_v.resize(3,-1);

    _par2_xdead_v.clear();
    _par2_xdead_v.resize(3,-1);

    _par2_end_pt_v.clear();
    _par2_end_pt_v.resize(3,-1);

    _par2_cosmic_dist_v.clear();
    _par2_cosmic_dist_v.resize(3,-1);

    _par2_cosmic_dist_end_v.clear();
    _par2_cosmic_dist_end_v.resize(3,-1);

    _par2_cosmic_end_dist_v.clear();
    _par2_cosmic_end_dist_v.resize(3,-1);

    _par2_cosmic_end_dist_end_v.clear();
    _par2_cosmic_end_dist_end_v.resize(3,-1);

    _par_theta         = nullptr;
    _par_phi           = nullptr;
    _par_length        = nullptr;
    _par_score         = nullptr;
    _par_dx1           = nullptr;
    _par_dy1           = nullptr;
    _par_dz1           = nullptr;
    _par_dx2           = nullptr;
    _par_dy2           = nullptr;
    _par_dz2           = nullptr;
    _par_nplanes       = nullptr;
    _par_planes_v      = nullptr;
    _par_xdead_v       = nullptr;

    _par_cosmic_dist_v     = nullptr;
    _par_cosmic_dist_end_v = nullptr;
    _par_end_pt_v          = nullptr;

    //
    // 2D information
    //

    _par1_n_polygons_U = -1*larocv::kINVALID_INT;
    _par1_n_polygons_V = -1*larocv::kINVALID_INT;
    _par1_n_polygons_Y = -1*larocv::kINVALID_INT;
    _par2_n_polygons_U = -1*larocv::kINVALID_INT;
    _par2_n_polygons_V = -1*larocv::kINVALID_INT;
    _par2_n_polygons_Y = -1*larocv::kINVALID_INT;
    _par_n_polygons = nullptr;

    _par1_linelength_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_linelength_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_linelength_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_linelength_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_linelength_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_linelength_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_linelength = nullptr;

    _par1_linefrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_linefrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_linefrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_linefrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_linefrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_linefrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_linefrac = nullptr;

    _par1_linefrac_empty_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_linefrac_empty_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_linefrac_empty_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_linefrac_empty_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_linefrac_empty_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_linefrac_empty_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_linefrac_empty = nullptr;

    _par1_linedx_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_linedx_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_linedx_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_linedx_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_linedx_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_linedx_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_linedx = nullptr;

    _par1_linedy_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_linedy_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_linedy_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_linedy_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_linedy_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_linedy_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_linedy = nullptr;

    _par1_line_vtx_density_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_vtx_density_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_vtx_density_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_density_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_density_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_density_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_vtx_density = nullptr;

    _par1_line_vtx_coverage_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_vtx_coverage_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_vtx_coverage_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_coverage_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_coverage_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_coverage_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_vtx_coverage = nullptr;

    _par1_line_vtx_charge_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_vtx_charge_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_vtx_charge_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_charge_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_charge_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_vtx_charge_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_vtx_charge = nullptr;

    _par1_line_mean_dist_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_mean_dist_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_mean_dist_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_mean_dist_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_mean_dist_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_mean_dist_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_mean_dist = nullptr;
    
    _par1_line_max_dist_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_max_dist_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_max_dist_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_max_dist_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_max_dist_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_max_dist_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_max_dist = nullptr;
    
    _par1_line_first_half_linefrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_first_half_linefrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_first_half_linefrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_first_half_linefrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_first_half_linefrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_first_half_linefrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_first_half_linefrac = nullptr;
    
    _par1_line_second_half_linefrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_second_half_linefrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_line_second_half_linefrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_second_half_linefrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_second_half_linefrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_line_second_half_linefrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_line_second_half_linefrac = nullptr;


    _par1_triangle_mid_to_edge_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_mid_to_edge_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_mid_to_edge_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_mid_to_edge_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_mid_to_edge_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_mid_to_edge_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_mid_to_edge = nullptr;

    _par1_triangle_height_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_height_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_height_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_height_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_height_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_height_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_height = nullptr;

    _par1_triangle_emptyarearatio_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_emptyarearatio_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_emptyarearatio_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_emptyarearatio_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_emptyarearatio_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_emptyarearatio_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_emptyarearatio = nullptr;

    _par1_triangle_emptyarea_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_emptyarea_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_emptyarea_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_emptyarea_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_emptyarea_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_emptyarea_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_emptyarea = nullptr;
    
    _par1_triangle_baselength_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_baselength_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_baselength_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_baselength_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_baselength_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_baselength_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_baselength = nullptr;

    _par1_triangle_area_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_area_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_area_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_area_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_area_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_area_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_area = nullptr;

    _par1_triangle_brem_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_brem_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_brem_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_brem_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_brem_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_brem_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_brem = nullptr;

    _par1_triangle_coverage_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_coverage_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_triangle_coverage_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_coverage_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_coverage_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_triangle_coverage_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_triangle_coverage = nullptr;

    _par1_brem_triangle_coverage_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_brem_triangle_coverage_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_brem_triangle_coverage_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_brem_triangle_coverage_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_brem_triangle_coverage_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_brem_triangle_coverage_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_brem_triangle_coverage = nullptr;

    _par1_expand_charge_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_expand_charge_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_expand_charge_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_expand_charge_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_expand_charge_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_expand_charge_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_expand_charge = nullptr;

    _par1_dqdx_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_dqdx_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_dqdx_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_dqdx = nullptr;

    _par1_dqdx_U_v.clear();
    _par1_dqdx_V_v.clear();
    _par1_dqdx_Y_v.clear();
    _par2_dqdx_U_v.clear();
    _par2_dqdx_V_v.clear();
    _par2_dqdx_Y_v.clear();
    _par_dqdx_v = nullptr;

    _par1_tdqdx_U_v.clear();
    _par1_tdqdx_V_v.clear();
    _par1_tdqdx_Y_v.clear();
    _par2_tdqdx_U_v.clear();
    _par2_tdqdx_V_v.clear();
    _par2_tdqdx_Y_v.clear();
    _par_tdqdx_v = nullptr;

    _par1_dqdx_step_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_dqdx_step_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_dqdx_step_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_step_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_step_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_step_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_dqdx_step = nullptr;

    _par1_dqdx_pitch_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_dqdx_pitch_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_dqdx_pitch_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_pitch_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_pitch_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_dqdx_pitch_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_dqdx_pitch = nullptr;
    
    _par1_brem_idx_U = -1.0*larocv::kINVALID_INT;
    _par1_brem_idx_V = -1.0*larocv::kINVALID_INT;
    _par1_brem_idx_Y = -1.0*larocv::kINVALID_INT;
    _par2_brem_idx_U = -1.0*larocv::kINVALID_INT;
    _par2_brem_idx_V = -1.0*larocv::kINVALID_INT;
    _par2_brem_idx_Y = -1.0*larocv::kINVALID_INT;
    _par_brem_idx = nullptr;
    
    _par1_length3d_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_length3d_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_length3d_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_length3d_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_length3d_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_length3d_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_length3d = nullptr;

    _par1_showerfrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_showerfrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_showerfrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_showerfrac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_showerfrac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_showerfrac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_showerfrac = nullptr;
    
    _par1_numberdefects_U_v.clear();
    _par1_numberdefects_V_v.clear();
    _par1_numberdefects_Y_v.clear();
    _par2_numberdefects_U_v.clear();
    _par2_numberdefects_V_v.clear();
    _par2_numberdefects_Y_v.clear();
    _par_numberdefects_v = nullptr;

    _par1_numberdefects_ns_U_v.clear();
    _par1_numberdefects_ns_V_v.clear();
    _par1_numberdefects_ns_Y_v.clear();
    _par2_numberdefects_ns_U_v.clear();
    _par2_numberdefects_ns_V_v.clear();
    _par2_numberdefects_ns_Y_v.clear();
    _par_numberdefects_ns_v = nullptr;
    
    _par1_largestdefect_U_v.clear();
    _par1_largestdefect_V_v.clear();
    _par1_largestdefect_Y_v.clear();
    _par2_largestdefect_U_v.clear();
    _par2_largestdefect_V_v.clear();
    _par2_largestdefect_Y_v.clear();
    _par_largestdefect_v = nullptr;

    _par1_smallestdefect_U_v.clear();
    _par1_smallestdefect_V_v.clear();
    _par1_smallestdefect_Y_v.clear();
    _par2_smallestdefect_U_v.clear();
    _par2_smallestdefect_V_v.clear();
    _par2_smallestdefect_Y_v.clear();
    _par_smallestdefect_v = nullptr;

    _par1_largestdefect_ns_U_v.clear();
    _par1_largestdefect_ns_V_v.clear();
    _par1_largestdefect_ns_Y_v.clear();
    _par2_largestdefect_ns_U_v.clear();
    _par2_largestdefect_ns_V_v.clear();
    _par2_largestdefect_ns_Y_v.clear();
    _par_largestdefect_ns_v = nullptr;

    _par1_smallestdefect_ns_U_v.clear();
    _par1_smallestdefect_ns_V_v.clear();
    _par1_smallestdefect_ns_Y_v.clear();
    _par2_smallestdefect_ns_U_v.clear();
    _par2_smallestdefect_ns_V_v.clear();
    _par2_smallestdefect_ns_Y_v.clear();
    _par_smallestdefect_ns_v = nullptr;

    _par1_emptyarearatio_U_v.clear();
    _par1_emptyarearatio_V_v.clear();
    _par1_emptyarearatio_Y_v.clear();
    _par2_emptyarearatio_U_v.clear();
    _par2_emptyarearatio_V_v.clear();
    _par2_emptyarearatio_Y_v.clear();
    _par_emptyarearatio_v = nullptr;

    _par1_emptyarea_U_v.clear();
    _par1_emptyarea_V_v.clear();
    _par1_emptyarea_Y_v.clear();
    _par2_emptyarea_U_v.clear();
    _par2_emptyarea_V_v.clear();
    _par2_emptyarea_Y_v.clear();
    _par_emptyarea_v = nullptr;

    _par1_pocketarea_U_v.clear();
    _par1_pocketarea_V_v.clear();
    _par1_pocketarea_Y_v.clear();
    _par2_pocketarea_U_v.clear();
    _par2_pocketarea_V_v.clear();
    _par2_pocketarea_Y_v.clear();
    _par_pocketarea_v = nullptr;

    _par1_pocketarea_ns_U_v.clear();
    _par1_pocketarea_ns_V_v.clear();
    _par1_pocketarea_ns_Y_v.clear();
    _par2_pocketarea_ns_U_v.clear();
    _par2_pocketarea_ns_V_v.clear();
    _par2_pocketarea_ns_Y_v.clear();
    _par_pocketarea_ns_v = nullptr;

    _par1_polyarea_U_v.clear();
    _par1_polyarea_V_v.clear();
    _par1_polyarea_Y_v.clear();
    _par2_polyarea_U_v.clear();
    _par2_polyarea_V_v.clear();
    _par2_polyarea_Y_v.clear();
    _par_polyarea_v = nullptr;

    _par1_polyperimeter_U_v.clear();
    _par1_polyperimeter_V_v.clear();
    _par1_polyperimeter_Y_v.clear();
    _par2_polyperimeter_U_v.clear();
    _par2_polyperimeter_V_v.clear();
    _par2_polyperimeter_Y_v.clear();
    _par_polyperimeter_v = nullptr;

    _par1_polycharge_U_v.clear();
    _par1_polycharge_V_v.clear();
    _par1_polycharge_Y_v.clear();
    _par2_polycharge_U_v.clear();
    _par2_polycharge_V_v.clear();
    _par2_polycharge_Y_v.clear();
    _par_polycharge_v = nullptr;

    _par1_polyedges_U_v.clear();
    _par1_polyedges_V_v.clear();
    _par1_polyedges_Y_v.clear();
    _par2_polyedges_U_v.clear();
    _par2_polyedges_V_v.clear();
    _par2_polyedges_Y_v.clear();
    _par_polyedges_v = nullptr;

    _par1_polybranches_U_v.clear();
    _par1_polybranches_V_v.clear();
    _par1_polybranches_Y_v.clear();
    _par2_polybranches_U_v.clear();
    _par2_polybranches_V_v.clear();
    _par2_polybranches_Y_v.clear();
    _par_polybranches_v = nullptr;

    _par1_showerfrac_U_v.clear();
    _par1_showerfrac_V_v.clear();
    _par1_showerfrac_Y_v.clear();
    _par2_showerfrac_U_v.clear();
    _par2_showerfrac_V_v.clear();
    _par2_showerfrac_Y_v.clear();
    _par_showerfrac_v = nullptr;

    _par1_electron_frac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_electron_frac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_electron_frac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_electron_frac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_electron_frac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_electron_frac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_electron_frac = nullptr;

    _par1_muon_frac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_muon_frac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_muon_frac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_muon_frac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_muon_frac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_muon_frac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_muon_frac = nullptr;

    _par1_proton_frac_U = -1.0*larocv::kINVALID_FLOAT;
    _par1_proton_frac_V = -1.0*larocv::kINVALID_FLOAT;
    _par1_proton_frac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par2_proton_frac_U = -1.0*larocv::kINVALID_FLOAT;
    _par2_proton_frac_V = -1.0*larocv::kINVALID_FLOAT;
    _par2_proton_frac_Y = -1.0*larocv::kINVALID_FLOAT;
    _par_proton_frac = nullptr;

    return;
  }

  void SelMichelID::ResizePlanePolygon(size_t sz) {

    if(!_par_numberdefects_v) {
      LLCV_CRITICAL() << "ptr to _par_numberdefects_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_numberdefects_ns_v) {
      LLCV_CRITICAL() << "ptr to _par_numberdefects_ns_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_largestdefect_v) {
      LLCV_CRITICAL() << "ptr to _par_largestdefect_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_smallestdefect_v) {
      LLCV_CRITICAL() << "ptr to _par_smallestdefect_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_largestdefect_ns_v) {
      LLCV_CRITICAL() << "ptr to _par_largestdefect_ns_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_smallestdefect_ns_v) {
      LLCV_CRITICAL() << "ptr to _par_smallestdefect_ns_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_emptyarearatio_v) {
      LLCV_CRITICAL() << "ptr to _par_emptyarearatio_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_emptyarea_v) {
      LLCV_CRITICAL() << "ptr to _par_emptyarea_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_pocketarea_v) {
      LLCV_CRITICAL() << "ptr to _par_pocketarea_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_pocketarea_ns_v) {
      LLCV_CRITICAL() << "ptr to _par_pocketarea_ns_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_polyarea_v) {
      LLCV_CRITICAL() << "ptr to _par_polyarea_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_polyperimeter_v) {
      LLCV_CRITICAL() << "ptr to _par_polyperimeter_v is null" << std::endl;
      throw llcv_err("die");
    }
    if(!_par_polycharge_v) {
      LLCV_CRITICAL() << "ptr to _par_polycharge_v is null" << std::endl;
      throw llcv_err("die");
    }

    if(!_par_polyedges_v) {
      LLCV_CRITICAL() << "ptr to _par_polyedges_v is null" << std::endl;
      throw llcv_err("die");
    }

    if(!_par_polybranches_v) {
      LLCV_CRITICAL() << "ptr to _par_branches_v is null" << std::endl;
      throw llcv_err("die");
    }

    if(!_par_showerfrac_v) {
      LLCV_CRITICAL() << "ptr to _par_showerfrac_v is null" << std::endl;
      throw llcv_err("die");
    }

    
    _par_numberdefects_v->resize(sz,-1.0*larocv::kINVALID_INT);

    _par_numberdefects_ns_v->resize(sz,-1.0*larocv::kINVALID_INT);

    _par_largestdefect_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_smallestdefect_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_largestdefect_ns_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_smallestdefect_ns_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_emptyarearatio_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);
    
    _par_emptyarea_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_pocketarea_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_pocketarea_ns_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_polyarea_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_polyperimeter_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);

    _par_polycharge_v->resize(sz,-1.0*larocv::kINVALID_FLOAT);
    
    _par_polyedges_v->resize(sz,-1.0*larocv::kINVALID_INT);

    _par_polybranches_v->resize(sz,-1.0*larocv::kINVALID_INT);

    _par_showerfrac_v->resize(sz,-1.0*larocv::kINVALID_INT);

    return;
  }
  
  void SelMichelID::SetParticle(size_t pid) {

    switch(pid) {
    case 0 : {
      _par_theta    = &_par1_theta;
      _par_phi      = &_par1_phi;
      _par_length   = &_par1_length;
      _par_score    = &_par1_score;
      _par_dx1      = &_par1_dx1;
      _par_dy1      = &_par1_dy1;
      _par_dz1      = &_par1_dz1;
      _par_dx2      = &_par1_dx2;
      _par_dy2      = &_par1_dy2;
      _par_dz2      = &_par1_dz2;
      _par_nplanes  = &_par1_nplanes;
      _par_planes_v = &_par1_planes_v;
      _par_xdead_v  = &_par1_xdead_v;
      _par_end_pt_v = &_par1_end_pt_v;

      _par_cosmic_dist_v         = &_par1_cosmic_dist_v;
      _par_cosmic_dist_end_v     = &_par1_cosmic_dist_end_v;
      _par_cosmic_end_dist_v     = &_par1_cosmic_end_dist_v;
      _par_cosmic_end_dist_end_v = &_par1_cosmic_end_dist_end_v;

      break;
    }
    case 1 : {
      _par_theta    = &_par2_theta;
      _par_phi      = &_par2_phi;
      _par_length   = &_par2_length;
      _par_score    = &_par2_score;
      _par_dx1      = &_par2_dx1;
      _par_dy1      = &_par2_dy1;
      _par_dz1      = &_par2_dz1;
      _par_dx2      = &_par2_dx2;
      _par_dy2      = &_par2_dy2;
      _par_dz2      = &_par2_dz2;
      _par_nplanes  = &_par2_nplanes;
      _par_planes_v = &_par2_planes_v;
      _par_xdead_v  = &_par2_xdead_v;
      _par_end_pt_v = &_par2_end_pt_v;

      _par_cosmic_dist_v         = &_par2_cosmic_dist_v;
      _par_cosmic_dist_end_v     = &_par2_cosmic_dist_end_v;
      _par_cosmic_end_dist_v     = &_par2_cosmic_end_dist_v;
      _par_cosmic_end_dist_end_v = &_par2_cosmic_end_dist_end_v;

      
      break;
    }
    default : { break; }
    }

    return;
  }

  void SelMichelID::SetParticlePlane(size_t pid, size_t plane) {

    LLCV_DEBUG() << "@pid=" << pid << " @plane=" << plane << std::endl;

    switch(pid) {
    case 0 : {
      switch(plane) {
      case 0 : {
	_par_n_polygons                = &_par1_n_polygons_U;
	_par_linelength                = &_par1_linelength_U;
	_par_linefrac                  = &_par1_linefrac_U;
	_par_linefrac_empty            = &_par1_linefrac_empty_U;
	_par_linedx                    = &_par1_linedx_U;
	_par_linedy                    = &_par1_linedy_U;
	_par_line_vtx_density          = &_par1_line_vtx_density_U;
	_par_line_vtx_coverage         = &_par1_line_vtx_coverage_U;
	_par_line_vtx_charge           = &_par1_line_vtx_charge_U;
	_par_line_mean_dist            = &_par1_line_mean_dist_U;
	_par_line_max_dist             = &_par1_line_max_dist_U;
	_par_line_first_half_linefrac  = &_par1_line_first_half_linefrac_U;
	_par_line_second_half_linefrac = &_par1_line_second_half_linefrac_U;
	_par_triangle_mid_to_edge      = &_par1_triangle_mid_to_edge_U;
	_par_triangle_height           = &_par1_triangle_height_U;
	_par_triangle_emptyarearatio   = &_par1_triangle_emptyarearatio_U;
	_par_triangle_emptyarea        = &_par1_triangle_emptyarea_U;
	_par_triangle_baselength       = &_par1_triangle_baselength_U;
	_par_triangle_area             = &_par1_triangle_area_U;
	_par_triangle_brem             = &_par1_triangle_brem_U;	
	_par_triangle_coverage         = &_par1_triangle_coverage_U;
	_par_brem_triangle_coverage    = &_par1_brem_triangle_coverage_U;
	_par_brem_idx                  = &_par1_brem_idx_U;
	_par_expand_charge             = &_par1_expand_charge_U;
	_par_dqdx                      = &_par1_dqdx_U;
	_par_dqdx_step                 = &_par1_dqdx_step_U;
	_par_dqdx_pitch                = &_par1_dqdx_pitch_U;
	_par_dqdx_v                    = &_par1_dqdx_U_v;
	_par_tdqdx_v                   = &_par1_tdqdx_U_v;
	_par_length3d                  = &_par1_length3d_U;
	_par_showerfrac                = &_par1_showerfrac_U;

	_par_numberdefects_v     = &_par1_numberdefects_U_v;
	_par_numberdefects_ns_v  = &_par1_numberdefects_ns_U_v;
	_par_largestdefect_v     = &_par1_largestdefect_U_v;
	_par_smallestdefect_v    = &_par1_smallestdefect_U_v;
	_par_largestdefect_ns_v  = &_par1_largestdefect_ns_U_v;
	_par_smallestdefect_ns_v = &_par1_smallestdefect_ns_U_v;
	_par_emptyarearatio_v    = &_par1_emptyarearatio_U_v;
	_par_emptyarea_v         = &_par1_emptyarea_U_v;
	_par_pocketarea_v        = &_par1_pocketarea_U_v;
	_par_pocketarea_ns_v     = &_par1_pocketarea_ns_U_v;
	_par_polyarea_v          = &_par1_polyarea_U_v;
	_par_polyperimeter_v     = &_par1_polyperimeter_U_v;
	_par_polycharge_v        = &_par1_polycharge_U_v;
	_par_polyedges_v         = &_par1_polyedges_U_v;
	_par_polybranches_v      = &_par1_polybranches_U_v;
	_par_showerfrac_v        = &_par1_showerfrac_U_v;

	break;
      }
      case 1: {
	_par_n_polygons                = &_par1_n_polygons_V;
	_par_linelength                = &_par1_linelength_V;
	_par_linefrac                  = &_par1_linefrac_V;
	_par_linefrac_empty            = &_par1_linefrac_empty_V;
	_par_linedx                    = &_par1_linedx_V;
	_par_linedy                    = &_par1_linedy_V;
	_par_line_vtx_density          = &_par1_line_vtx_density_V;
	_par_line_vtx_coverage         = &_par1_line_vtx_coverage_V;
	_par_line_vtx_charge           = &_par1_line_vtx_charge_V;
	_par_line_mean_dist            = &_par1_line_mean_dist_V;
	_par_line_max_dist             = &_par1_line_max_dist_V;
	_par_line_first_half_linefrac  = &_par1_line_first_half_linefrac_V;
	_par_line_second_half_linefrac = &_par1_line_second_half_linefrac_V;
	_par_triangle_mid_to_edge      = &_par1_triangle_mid_to_edge_V;
	_par_triangle_height           = &_par1_triangle_height_V;
	_par_triangle_emptyarearatio   = &_par1_triangle_emptyarearatio_V;
	_par_triangle_emptyarea        = &_par1_triangle_emptyarea_V;
	_par_triangle_baselength       = &_par1_triangle_baselength_V;
	_par_triangle_area             = &_par1_triangle_area_V;
	_par_triangle_brem             = &_par1_triangle_brem_V;
	_par_triangle_coverage         = &_par1_triangle_coverage_V;
	_par_brem_triangle_coverage    = &_par1_brem_triangle_coverage_V;
	_par_brem_idx                  = &_par1_brem_idx_V;
	_par_expand_charge             = &_par1_expand_charge_V;
	_par_dqdx                      = &_par1_dqdx_V;
	_par_dqdx_step                 = &_par1_dqdx_step_V;
	_par_dqdx_pitch                = &_par1_dqdx_pitch_V;
	_par_dqdx_v                    = &_par1_dqdx_V_v;
	_par_tdqdx_v                   = &_par1_tdqdx_V_v;
	_par_length3d                  = &_par1_length3d_V;
	_par_showerfrac                = &_par1_showerfrac_V;

	_par_numberdefects_v     = &_par1_numberdefects_V_v;
	_par_numberdefects_ns_v  = &_par1_numberdefects_ns_V_v;
	_par_largestdefect_v     = &_par1_largestdefect_V_v;
	_par_smallestdefect_v    = &_par1_smallestdefect_V_v;
	_par_largestdefect_ns_v  = &_par1_largestdefect_ns_V_v;
	_par_smallestdefect_ns_v = &_par1_smallestdefect_ns_V_v;
	_par_emptyarearatio_v    = &_par1_emptyarearatio_V_v;
	_par_emptyarea_v         = &_par1_emptyarea_V_v;
	_par_pocketarea_v        = &_par1_pocketarea_V_v;
	_par_pocketarea_ns_v     = &_par1_pocketarea_ns_V_v;
	_par_polyarea_v          = &_par1_polyarea_V_v;
	_par_polyperimeter_v     = &_par1_polyperimeter_V_v;
	_par_polycharge_v        = &_par1_polycharge_V_v;
	_par_polyedges_v         = &_par1_polyedges_V_v;
	_par_polybranches_v      = &_par1_polybranches_V_v;
	_par_showerfrac_v        = &_par1_showerfrac_V_v;

	break;
      }
      case 2: {
	_par_n_polygons                = &_par1_n_polygons_Y;
	_par_linelength                = &_par1_linelength_Y;
	_par_linefrac                  = &_par1_linefrac_Y;
	_par_linefrac_empty            = &_par1_linefrac_empty_Y;
	_par_linedx                    = &_par1_linedx_Y;
	_par_linedy                    = &_par1_linedy_Y;
	_par_line_vtx_density          = &_par1_line_vtx_density_Y;
	_par_line_vtx_coverage         = &_par1_line_vtx_coverage_Y;
	_par_line_vtx_charge           = &_par1_line_vtx_charge_Y;
	_par_line_mean_dist            = &_par1_line_mean_dist_Y;
	_par_line_max_dist             = &_par1_line_max_dist_Y;
	_par_line_first_half_linefrac  = &_par1_line_first_half_linefrac_Y;
	_par_line_second_half_linefrac = &_par1_line_second_half_linefrac_Y;
	_par_triangle_mid_to_edge      = &_par1_triangle_mid_to_edge_Y;
	_par_triangle_height           = &_par1_triangle_height_Y;
	_par_triangle_emptyarearatio   = &_par1_triangle_emptyarearatio_Y;
	_par_triangle_emptyarea        = &_par1_triangle_emptyarea_Y;
	_par_triangle_baselength       = &_par1_triangle_baselength_Y;
	_par_triangle_area             = &_par1_triangle_area_Y;
	_par_triangle_brem             = &_par1_triangle_brem_Y;
	_par_triangle_coverage         = &_par1_triangle_coverage_Y;
	_par_brem_triangle_coverage    = &_par1_brem_triangle_coverage_Y;
	_par_brem_idx                  = &_par1_brem_idx_Y;
	_par_expand_charge             = &_par1_expand_charge_Y;
	_par_dqdx                      = &_par1_dqdx_Y;
	_par_dqdx_step                 = &_par1_dqdx_step_Y;
	_par_dqdx_pitch                = &_par1_dqdx_pitch_Y;
	_par_dqdx_v                    = &_par1_dqdx_Y_v;
	_par_tdqdx_v                   = &_par1_tdqdx_Y_v;
	_par_length3d                  = &_par1_length3d_Y;
	_par_showerfrac                = &_par1_showerfrac_Y;

	_par_numberdefects_v     = &_par1_numberdefects_Y_v;
	_par_numberdefects_ns_v  = &_par1_numberdefects_ns_Y_v;
	_par_largestdefect_v     = &_par1_largestdefect_Y_v;
	_par_smallestdefect_v    = &_par1_smallestdefect_Y_v;
	_par_largestdefect_ns_v  = &_par1_largestdefect_ns_Y_v;
	_par_smallestdefect_ns_v = &_par1_smallestdefect_ns_Y_v;
	_par_emptyarearatio_v    = &_par1_emptyarearatio_Y_v;
	_par_emptyarea_v         = &_par1_emptyarea_Y_v;
	_par_pocketarea_v        = &_par1_pocketarea_Y_v;
	_par_pocketarea_ns_v     = &_par1_pocketarea_ns_Y_v;
	_par_polyarea_v          = &_par1_polyarea_Y_v;
	_par_polyperimeter_v     = &_par1_polyperimeter_Y_v;
	_par_polycharge_v        = &_par1_polycharge_Y_v;
	_par_polyedges_v         = &_par1_polyedges_Y_v;
	_par_polybranches_v      = &_par1_polybranches_Y_v;
	_par_showerfrac_v        = &_par1_showerfrac_Y_v;

	break;
      }
      default : { break; }
      } // end plane
      break;
    } // end particle 1

    case 1 : {
      switch(plane) {
      case 0 : {
	_par_n_polygons                = &_par2_n_polygons_U;
	_par_linelength                = &_par2_linelength_U;
	_par_linefrac                  = &_par2_linefrac_U;
	_par_linefrac_empty            = &_par2_linefrac_empty_U;
	_par_linedx                    = &_par2_linedx_U;
	_par_linedy                    = &_par2_linedy_U;
	_par_line_vtx_density          = &_par2_line_vtx_density_U;
	_par_line_vtx_coverage         = &_par2_line_vtx_coverage_U;
	_par_line_vtx_charge           = &_par2_line_vtx_charge_U;
	_par_line_mean_dist            = &_par2_line_mean_dist_U;
	_par_line_max_dist             = &_par2_line_max_dist_U;
	_par_line_first_half_linefrac  = &_par2_line_first_half_linefrac_U;
	_par_line_second_half_linefrac = &_par2_line_second_half_linefrac_U;
	_par_triangle_mid_to_edge      = &_par2_triangle_mid_to_edge_U;
	_par_triangle_height           = &_par2_triangle_height_U;
	_par_triangle_emptyarearatio   = &_par2_triangle_emptyarearatio_U;
	_par_triangle_emptyarea        = &_par2_triangle_emptyarea_U;
	_par_triangle_baselength       = &_par2_triangle_baselength_U;
	_par_triangle_area             = &_par2_triangle_area_U;
	_par_triangle_brem             = &_par2_triangle_brem_U;
	_par_triangle_coverage         = &_par2_triangle_coverage_U;
	_par_brem_triangle_coverage    = &_par2_brem_triangle_coverage_U;
	_par_brem_idx                  = &_par2_brem_idx_U;
	_par_expand_charge             = &_par2_expand_charge_U;
	_par_dqdx                      = &_par2_dqdx_U;
	_par_dqdx_step                 = &_par2_dqdx_step_U;
	_par_dqdx_pitch                = &_par2_dqdx_pitch_U;
	_par_dqdx_v                    = &_par2_dqdx_U_v;
	_par_tdqdx_v                   = &_par2_tdqdx_U_v;
	_par_length3d                  = &_par2_length3d_U;
	_par_showerfrac                = &_par2_showerfrac_U;

	_par_numberdefects_v     = &_par2_numberdefects_U_v;
	_par_numberdefects_ns_v  = &_par2_numberdefects_ns_U_v;
	_par_largestdefect_v     = &_par2_largestdefect_U_v;
	_par_smallestdefect_v    = &_par2_smallestdefect_U_v;
	_par_largestdefect_ns_v  = &_par2_largestdefect_ns_U_v;
	_par_smallestdefect_ns_v = &_par2_smallestdefect_ns_U_v;
	_par_emptyarearatio_v    = &_par2_emptyarearatio_U_v;
	_par_emptyarea_v         = &_par2_emptyarea_U_v;
	_par_pocketarea_v        = &_par2_pocketarea_U_v;
	_par_pocketarea_ns_v     = &_par2_pocketarea_ns_U_v;
	_par_polyarea_v          = &_par2_polyarea_U_v;
	_par_polyperimeter_v     = &_par2_polyperimeter_U_v;
	_par_polycharge_v        = &_par2_polycharge_U_v;
	_par_polyedges_v         = &_par2_polyedges_U_v;
	_par_polybranches_v      = &_par2_polybranches_U_v;
	_par_showerfrac_v        = &_par2_showerfrac_U_v;

	break;
      }
      case 1: {
	_par_n_polygons                = &_par2_n_polygons_V;
	_par_linelength                = &_par2_linelength_V;
	_par_linefrac                  = &_par2_linefrac_V;
	_par_linefrac_empty            = &_par2_linefrac_empty_V;
	_par_linedx                    = &_par2_linedx_V;
	_par_linedy                    = &_par2_linedy_V;
	_par_line_vtx_density          = &_par2_line_vtx_density_V;
	_par_line_vtx_coverage         = &_par2_line_vtx_coverage_V;
	_par_line_vtx_charge           = &_par2_line_vtx_charge_V;
	_par_line_mean_dist            = &_par2_line_mean_dist_V;
	_par_line_max_dist             = &_par2_line_max_dist_V;
	_par_line_first_half_linefrac  = &_par2_line_first_half_linefrac_V;
	_par_line_second_half_linefrac = &_par2_line_second_half_linefrac_V;
	_par_triangle_mid_to_edge      = &_par2_triangle_mid_to_edge_V;
	_par_triangle_height           = &_par2_triangle_height_V;
	_par_triangle_emptyarearatio   = &_par2_triangle_emptyarearatio_V;
	_par_triangle_emptyarea        = &_par2_triangle_emptyarea_V;
	_par_triangle_baselength       = &_par2_triangle_baselength_V;
	_par_triangle_area             = &_par2_triangle_area_V;
	_par_triangle_brem             = &_par2_triangle_brem_V;
	_par_triangle_coverage         = &_par2_triangle_coverage_V;
	_par_brem_triangle_coverage    = &_par2_brem_triangle_coverage_V;
	_par_brem_idx                  = &_par2_brem_idx_V;
	_par_expand_charge             = &_par2_expand_charge_V;
	_par_dqdx                      = &_par2_dqdx_V;
	_par_dqdx_step                 = &_par2_dqdx_step_V;
	_par_dqdx_pitch                = &_par2_dqdx_pitch_V;
	_par_dqdx_v                    = &_par2_dqdx_V_v;
	_par_tdqdx_v                   = &_par2_tdqdx_V_v;
	_par_length3d                  = &_par2_length3d_V;
	_par_showerfrac                = &_par2_showerfrac_V;

	_par_numberdefects_v     = &_par2_numberdefects_V_v;
	_par_numberdefects_ns_v  = &_par2_numberdefects_ns_V_v;
	_par_largestdefect_v     = &_par2_largestdefect_V_v;
	_par_smallestdefect_v    = &_par2_smallestdefect_V_v;
	_par_largestdefect_ns_v  = &_par2_largestdefect_ns_V_v;
	_par_smallestdefect_ns_v = &_par2_smallestdefect_ns_V_v;
	_par_emptyarearatio_v    = &_par2_emptyarearatio_V_v;
	_par_emptyarea_v         = &_par2_emptyarea_V_v;
	_par_pocketarea_v        = &_par2_pocketarea_V_v;
	_par_pocketarea_ns_v     = &_par2_pocketarea_ns_V_v;
	_par_polyarea_v          = &_par2_polyarea_V_v;
	_par_polyperimeter_v     = &_par2_polyperimeter_V_v;
	_par_polycharge_v        = &_par2_polycharge_V_v;
	_par_polyedges_v         = &_par2_polyedges_V_v;
	_par_polybranches_v      = &_par2_polybranches_V_v;
	_par_showerfrac_v        = &_par2_showerfrac_V_v;

	break;
      }
      case 2: {
	_par_n_polygons                = &_par2_n_polygons_Y;
	_par_linelength                = &_par2_linelength_Y;
	_par_linefrac                  = &_par2_linefrac_Y;
	_par_linefrac_empty            = &_par2_linefrac_empty_Y;
	_par_linedx                    = &_par2_linedx_Y;
	_par_linedy                    = &_par2_linedy_Y;
	_par_line_vtx_density          = &_par2_line_vtx_density_Y;
	_par_line_vtx_coverage         = &_par2_line_vtx_coverage_Y;
	_par_line_vtx_charge           = &_par2_line_vtx_charge_Y;
	_par_line_mean_dist            = &_par2_line_mean_dist_Y;
	_par_line_max_dist             = &_par2_line_max_dist_Y;
	_par_line_first_half_linefrac  = &_par2_line_first_half_linefrac_Y;
	_par_line_second_half_linefrac = &_par2_line_second_half_linefrac_Y;
	_par_triangle_mid_to_edge      = &_par2_triangle_mid_to_edge_Y;
	_par_triangle_height           = &_par2_triangle_height_Y;
	_par_triangle_emptyarearatio   = &_par2_triangle_emptyarearatio_Y;
	_par_triangle_emptyarea        = &_par2_triangle_emptyarea_Y;
	_par_triangle_baselength       = &_par2_triangle_baselength_Y;
	_par_triangle_area             = &_par2_triangle_area_Y;
	_par_triangle_brem             = &_par2_triangle_brem_Y;
	_par_triangle_coverage         = &_par2_triangle_coverage_Y;
	_par_brem_triangle_coverage    = &_par2_brem_triangle_coverage_Y;
	_par_brem_idx                  = &_par2_brem_idx_Y;
	_par_expand_charge             = &_par2_expand_charge_Y;
	_par_dqdx                      = &_par2_dqdx_Y;
	_par_dqdx_step                 = &_par2_dqdx_step_Y;
	_par_dqdx_pitch                = &_par2_dqdx_pitch_Y;
	_par_dqdx_v                    = &_par2_dqdx_Y_v;
	_par_tdqdx_v                   = &_par2_tdqdx_Y_v;
	_par_length3d                  = &_par2_length3d_Y;
	_par_showerfrac                = &_par2_showerfrac_Y;

	_par_numberdefects_v     = &_par2_numberdefects_Y_v;
	_par_numberdefects_ns_v  = &_par2_numberdefects_ns_Y_v;
	_par_largestdefect_v     = &_par2_largestdefect_Y_v;
	_par_smallestdefect_v    = &_par2_smallestdefect_Y_v;
	_par_largestdefect_ns_v  = &_par2_largestdefect_ns_Y_v;
	_par_smallestdefect_ns_v = &_par2_smallestdefect_ns_Y_v;
	_par_emptyarearatio_v    = &_par2_emptyarearatio_Y_v;
	_par_emptyarea_v         = &_par2_emptyarea_Y_v;
	_par_pocketarea_v        = &_par2_pocketarea_Y_v;
	_par_pocketarea_ns_v     = &_par2_pocketarea_ns_Y_v;
	_par_polyarea_v          = &_par2_polyarea_Y_v;
	_par_polyperimeter_v     = &_par2_polyperimeter_Y_v;
	_par_polycharge_v        = &_par2_polycharge_Y_v;
	_par_polyedges_v         = &_par2_polyedges_Y_v;
	_par_polybranches_v      = &_par2_polybranches_Y_v;
	_par_showerfrac_v        = &_par2_showerfrac_Y_v;

	break;
      }
      default : { break; }
      } // end plane
      break;
    } // end particle 2
    default : { break; }
    } // end particle
    return;
  }

  void SelMichelID::SetSegmentPlane(size_t pid, size_t plane) {

    LLCV_DEBUG() << "@pid=" << pid << " @plane=" << plane << std::endl;

    switch(pid) {
    case 0 : {
      switch(plane) {
      case 0 : {
	_par_electron_frac = &_par1_electron_frac_U;
	_par_muon_frac     = &_par1_muon_frac_U;
	_par_proton_frac   = &_par1_proton_frac_U;
	break;
      }
      case 1: {
	_par_electron_frac = &_par1_electron_frac_V;
	_par_muon_frac     = &_par1_muon_frac_V;
	_par_proton_frac   = &_par1_proton_frac_V;
	break;
      }
      case 2: {
	_par_electron_frac = &_par1_electron_frac_Y;
	_par_muon_frac     = &_par1_muon_frac_Y;
	_par_proton_frac   = &_par1_proton_frac_Y;
	break;
      }
      default : { break; }
      } // end plane
      break;
    } // end particle 1

    case 1 : {
      switch(plane) {
      case 0 : {
	_par_electron_frac = &_par2_electron_frac_U;
	_par_muon_frac     = &_par2_muon_frac_U;
	_par_proton_frac   = &_par2_proton_frac_U;
	break;
      }
      case 1: {
	_par_electron_frac = &_par2_electron_frac_V;
	_par_muon_frac     = &_par2_muon_frac_V;
	_par_proton_frac   = &_par2_proton_frac_V;
	break;
      }
      case 2: {
	_par_electron_frac = &_par2_electron_frac_Y;
	_par_muon_frac     = &_par2_muon_frac_Y;
	_par_proton_frac   = &_par2_proton_frac_Y;
	break;
      }
      default : { break; }
      } // end plane
      break;
    } // end particle 2
    default : { break; }
    } // end particle

    return;
  }

  
  void SelMichelID::Finalize() {
    LLCV_DEBUG() << "start" << std::endl;
    _outtree->Write();
    LLCV_DEBUG() << "end" << std::endl;
  }

}


#endif
  
