#ifndef __SELTRACKSCATTER_H__
#define __SELTRACKSCATTER_H__

#include "InterTool_Core/InterSelBase.h"
#include "LArOpenCV/ImageCluster/AlgoClass/DefectBreaker.h"
#include "LArOpenCV/ImageCluster/AlgoClass/PixelScan3D.h"

#include "DBSCAN.h"
#include "Skeletonize.h"

#include "Object3D.h"
#include <array>

#include "TStopwatch.h"

namespace llcv {
  

  class SelTrackScatter : public InterSelBase { 

  public:

  SelTrackScatter(std::string name="SelTrackScatter") : InterSelBase(name), _outtree(nullptr) {}
    ~SelTrackScatter(){}
    
    void Configure (const larcv::PSet &pset);
    void Initialize();
    double Select();
    void Finalize();
    
  private:
    TTree* _outtree;

    //
    // Members
    //
  private:

    DBSCAN _DBSCAN;
    Skeletonize _Skeletonize;

    size_t _cropx;
    size_t _cropy;

    larocv::DefectBreaker _LR_DefectBreaker;
    larocv::DefectBreaker _SR_DefectBreaker;

    larocv::PixelScan3D _PixelScan3D;

    TStopwatch _twatch;

    bool _debug;
    bool _fill2d;
    bool _skeletonize;
    bool _sub_skeleton;
    bool _allow_dead_image;

    larocv::GEO2D_ContourArray_t FindAndMaskVertex(const cv::Mat& mat,const cv::Point_<int> vertex);
    larocv::GEO2D_ContourArray_t FindAndBreakVertex(const cv::Mat& mat,const cv::Point_<int> vertex);
    
    std::vector<std::vector<size_t> > AssociateToTracks(const larocv::GEO2D_ContourArray_t& ctor_v,
							const std::vector<larocv::GEO2D_ContourArray_t>& track_ctor_vv,
							const size_t plane);

    bool ContainsTrack(const std::vector<size_t>& tin_v);

    float TrackFraction(const std::vector<size_t>& tin_v, size_t tid);
    float TrackLength(const larlite::track& track);
    std::pair<float,float> TrackAngle(const larlite::track& track); 

    std::vector<int> Cluster(const Object3D& obj);

    std::array<float,4> ComputeMeans(const std::vector<float>& data_v);

    float Average(const std::vector<float>& data_v, size_t start, size_t end);

    int CountClusters(const std::vector<int>& cid_v);

    void ResizeOutput(size_t sz);


    //
    // TTree
    //
  private:

    std::vector<std::vector<float> > _track_x_vv;
    std::vector<std::vector<float> > _track_y_vv;
    std::vector<std::vector<float> > _track_z_vv;

    std::vector<std::vector<float> > _shower_x_vv;
    std::vector<std::vector<float> > _shower_y_vv;
    std::vector<std::vector<float> > _shower_z_vv;

    std::vector<std::vector<float> > _shower_p0_x_vv;
    std::vector<std::vector<float> > _shower_p0_y_vv;
    std::vector<std::vector<float> > _shower_p0_z_vv;

    std::vector<std::vector<float> > _shower_p1_x_vv;
    std::vector<std::vector<float> > _shower_p1_y_vv;
    std::vector<std::vector<float> > _shower_p1_z_vv;

    std::vector<std::vector<float> > _shower_p2_x_vv;
    std::vector<std::vector<float> > _shower_p2_y_vv;
    std::vector<std::vector<float> > _shower_p2_z_vv;

    std::vector<std::vector<float> > _shower_skel_x_vv;
    std::vector<std::vector<float> > _shower_skel_y_vv;
    std::vector<std::vector<float> > _shower_skel_z_vv;

    std::vector<std::vector<float> > _shower_start_x_vv;
    std::vector<std::vector<float> > _shower_start_y_vv;
    std::vector<std::vector<float> > _shower_start_z_vv;

    std::vector<std::vector<float> > _shower_end_x_vv;
    std::vector<std::vector<float> > _shower_end_y_vv;
    std::vector<std::vector<float> > _shower_end_z_vv;

    std::vector<std::vector<float> > _shower_center_x_vv;
    std::vector<std::vector<float> > _shower_center_y_vv;
    std::vector<std::vector<float> > _shower_center_z_vv;

    std::vector<std::vector<int> > _shower_cid_vv;
    
    std::vector<std::vector<float> > _shower_pca_dev_vv;
    std::vector<std::vector<float> > _shower_trk_dev_vv;

    std::vector<std::vector<float> > _shower_edge1_x_vv;
    std::vector<std::vector<float> > _shower_edge1_y_vv;
    std::vector<std::vector<float> > _shower_edge1_z_vv;

    std::vector<std::vector<float> > _shower_edge2_x_vv;
    std::vector<std::vector<float> > _shower_edge2_y_vv;
    std::vector<std::vector<float> > _shower_edge2_z_vv;

    //
    // per track -- shower
    //
    std::vector<int>  _shower3D_n_points_v;

    std::vector<float> _shower3D_length_v;
    std::vector<float> _shower3D_width_v;
    std::vector<float> _shower3D_width1_v;
    std::vector<float> _shower3D_width2_v;

    std::vector<float> _shower3D_theta_v;
    std::vector<float> _shower3D_phi_v;

    std::vector<float> _shower3D_opening_v;
    std::vector<float> _shower3D_opening1_v;
    std::vector<float> _shower3D_opening2_v;

    std::vector<float> _shower3D_pca_mean_dev_v;
    std::vector<float> _shower3D_start_pca_mean_dev_v;
    std::vector<float> _shower3D_middle_pca_mean_dev_v;
    std::vector<float> _shower3D_end_pca_mean_dev_v;

    std::vector<float> _shower3D_track_mean_dev_v;
    std::vector<float> _shower3D_start_track_mean_dev_v;
    std::vector<float> _shower3D_middle_track_mean_dev_v;
    std::vector<float> _shower3D_end_track_mean_dev_v;

    std::vector<int> _shower3D_n_clusters_v;
    
    //
    // per track per 3D shower -- cluster
    //
    std::vector<std::vector<int> > _shower3D_cluster_n_points_vv;

    std::vector<std::vector<float> > _shower3D_cluster_length_vv;
    std::vector<std::vector<float> > _shower3D_cluster_width_vv;
    std::vector<std::vector<float> > _shower3D_cluster_width1_vv;
    std::vector<std::vector<float> > _shower3D_cluster_width2_vv;
    
    std::vector<std::vector<float> > _shower3D_cluster_theta_vv;
    std::vector<std::vector<float> > _shower3D_cluster_phi_vv;
    
    std::vector<std::vector<float> > _shower3D_cluster_opening_vv;
    std::vector<std::vector<float> > _shower3D_cluster_opening1_vv;
    std::vector<std::vector<float> > _shower3D_cluster_opening2_vv;
    
    std::vector<std::vector<float> > _shower3D_cluster_distance_vv;
    
    //
    // per track per cluster -- 2D stuffs
    //
    std::vector<int> _shower2D_n_clusters_U_v;
    std::vector<int> _shower2D_n_clusters_V_v;
    std::vector<int> _shower2D_n_clusters_Y_v;

    std::vector<std::vector<float> > _shower2D_area_U_vv;
    std::vector<std::vector<float> > _shower2D_length_U_vv;
    std::vector<std::vector<float> > _shower2D_width_U_vv;
    std::vector<std::vector<float> > _shower2D_npixel_U_vv;
    std::vector<std::vector<float> > _shower2D_qsum_U_vv;

    std::vector<std::vector<float> > _shower2D_area_V_vv;
    std::vector<std::vector<float> > _shower2D_length_V_vv;
    std::vector<std::vector<float> > _shower2D_width_V_vv;
    std::vector<std::vector<float> > _shower2D_npixel_V_vv;
    std::vector<std::vector<float> > _shower2D_qsum_V_vv;

    std::vector<std::vector<float> > _shower2D_area_Y_vv;
    std::vector<std::vector<float> > _shower2D_length_Y_vv;
    std::vector<std::vector<float> > _shower2D_width_Y_vv;
    std::vector<std::vector<float> > _shower2D_npixel_Y_vv;
    std::vector<std::vector<float> > _shower2D_qsum_Y_vv;

    std::vector<std::vector<int> > _shower2D_n_defects_U_vv;
    std::vector<std::vector<int> > _shower2D_n_defects_V_vv;
    std::vector<std::vector<int> > _shower2D_n_defects_Y_vv;

  };

}


#endif
