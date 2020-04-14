#include "ShowerRecoUtil.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#include "TMath.h"

#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/EventPGraph.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/track.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {
namespace ssnetshowerreco {

  bool ShowerRecoUtil::_setup_numpy = false;

  PyObject* ShowerRecoUtil::process_event_get_ndarray( larcv::IOManager& iolcv )
  {
    
    if ( !ShowerRecoUtil::_setup_numpy ) {
      ShowerRecoUtil::_setup_numpy = true;            
      import_array1(0);
    }


    // get image2d
    std::vector< std::vector<larcv::Image2D> > adccrop_vv;
    std::vector< std::vector<larcv::Image2D> > statuscrop_vv;

    process_event( iolcv, adccrop_vv, statuscrop_vv );
    

    PyObject* vtx_list = PyList_New(0);

    for ( int ivtx=0; ivtx<(int)adccrop_vv.size(); ivtx++ ) {
      auto& adc_v    = adccrop_vv[ivtx];
      auto& status_v = statuscrop_vv[ivtx];

      // make dict
      PyObject *d = PyDict_New();
      PyObject *str_adc_key    = Py_BuildValue("s", "adc" );
      PyObject *str_status_key = Py_BuildValue("s", "status" );

      PyObject *adc_list = PyList_New(0);
      PyObject *status_list = PyList_New(0);

      if ( adc_v.size()>0 ) {

        for (int p=0; p<3; p++ ) {

          // make numpy array
          npy_intp dims1[] = { (long)adc_v[0].meta().rows(), (long)adc_v[0].meta().cols() };
          npy_intp dims2[] = { (long)status_v[0].meta().rows(), (long)status_v[0].meta().cols() };
          
          PyArrayObject* adcimg    = (PyArrayObject*)PyArray_SimpleNew( 2, dims1, NPY_FLOAT );
          PyArrayObject* statusimg = (PyArrayObject*)PyArray_SimpleNew( 2, dims2, NPY_FLOAT );
          
          for (size_t r=0; r<adc_v[p].meta().rows(); r++ ) {
            for (size_t c=0; c<adc_v[p].meta().cols(); c++ ) {
              *((float*)PyArray_GETPTR2( adcimg, r, c ))    = adc_v[p].pixel(r,c);
              *((float*)PyArray_GETPTR2( statusimg, r, c )) = status_v[p].pixel(r,c);
            }
          }

          PyList_Append( adc_list,    (PyObject*)adcimg );
          PyList_Append( status_list, (PyObject*)statusimg );

          Py_DECREF( adcimg );
          Py_DECREF( statusimg );
          
        }//end of loop over planes
      }
      
      PyDict_SetItem(d, str_adc_key,    (PyObject*)adc_list);        
      PyDict_SetItem(d, str_status_key, (PyObject*)status_list );
    
      Py_DECREF( str_adc_key );
      Py_DECREF( str_status_key );
      Py_DECREF( adc_list );
      Py_DECREF( status_list );

      PyList_Append( vtx_list, d );
      Py_DECREF( d );
      
    }//end of vertex loop


    return vtx_list;
  }
  
  bool ShowerRecoUtil::process_event( larcv::IOManager& iolcv,
                                      std::vector< std::vector<larcv::Image2D> >& adccrop_vv,
                                      std::vector< std::vector<larcv::Image2D> >& statuscrop_vv )
  {


    adccrop_vv.clear();
    statuscrop_vv.clear();

    std::string _adc_tree_name = "wire";
    std::string _chstatus_tree_name = "wire";
    std::string _ssnet_shower_image_stem = "ubspurn";
    std::string _vertex_tree_name = "test";

    larcv::EventImage2D* ev_adc
      = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _adc_tree_name );
    const std::vector<larcv::Image2D>& adc_v = ev_adc->Image2DArray();

    larcv::EventImage2D* ev_shower_score[3] = { nullptr };
    std::vector<larcv::Image2D> shower_img_v;
    for ( size_t p=0; p<3; p++ ) {
      char treename[50];
      sprintf( treename, "%s_plane%d", _ssnet_shower_image_stem.c_str(), (int)p );
      ev_shower_score[p] =
        (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, treename );
      shower_img_v.push_back( ev_shower_score[p]->Image2DArray()[0] );
    }
    
    const auto ev_wirestatus =
      (larcv::EventChStatus*)iolcv.get_data( larcv::kProductChStatus, _chstatus_tree_name );
    
    larcv::EventPGraph* ev_vtx
      = (larcv::EventPGraph*)iolcv.get_data( larcv::kProductPGraph, _vertex_tree_name );

    // get candidate vertices, make crops around said vertex

    for ( auto const& pgraph : ev_vtx->PGraphArray() ) {

      std::vector<larcv::Image2D> adccrop_v;
      std::vector<larcv::Image2D> statuscrop_v;

      if ( pgraph.ParticleArray().size()>0 ) {
      
        auto const& roi = pgraph.ParticleArray().front();
        std::vector<double> vtx3d = { roi.X(), roi.Y(), roi.Z() };
        std::cout << "Vertex Pos (" << vtx3d[0] << "," << vtx3d[1] << "," << vtx3d[2] << ")" << std::endl;
        
        if ((vtx3d[0]>    0.001) && (vtx3d[0] <  255.999)
            && (vtx3d[1]> -116.499) && (vtx3d[1] < 116.499)
            && (vtx3d[2]>    0.001) && (vtx3d[2] < 1036.999) ) {
          
          makeVertexImageCrop( adc_v,
                               shower_img_v,
                               *ev_wirestatus,
                               pgraph,
                               10.0,
                               0.5,
                               adccrop_v,
                               statuscrop_v );
          
        }
        else {
          std::cout<<"Vertex Outside fiducial volume!"<<std::endl;
        }
        
      }
      
      adccrop_vv.emplace_back( std::move(adccrop_v));
      statuscrop_vv.emplace_back( std::move(statuscrop_v) );
      
    }//end of loop over vertices/pgraph
    
    
    return true;
  }
  
  void ShowerRecoUtil::makeVertexImageCrop( const std::vector<larcv::Image2D>& adc_v,
                                            const std::vector<larcv::Image2D>& shower_img_v,
                                            const larcv::EventChStatus& ev_wirestatus,
                                            const larcv::PGraph& vtxinfo,
                                            const float adc_threshold,
                                            const float showerscore_threshold,
                                            std::vector<larcv::Image2D>& cropimg_v,
                                            std::vector<larcv::Image2D>& statusimg_v )
  {

    cropimg_v.clear();
    statusimg_v.clear();

    auto const& roi = vtxinfo.ParticleArray().front();    
    std::vector<double> vtx3d = { roi.X(), roi.Y(), roi.Z() };
    std::cout << "Cropping image for Vertex Pos (" << vtx3d[0] << "," << vtx3d[1] << "," << vtx3d[2] << ")" << std::endl;

    /// make image coordinates
    std::vector<int> imgcoord_v(4);  // (U,V,Y,tick)
    imgcoord_v[3] = 3200 + vtx3d[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
    std::cout << "Image coordinate vertex positions (u,v,y,tick):\n";
    for ( size_t p=0; p<3; p++ ) {
      imgcoord_v[p] = larutil::Geometry::GetME()->NearestWire( vtx3d, (int)p );
      std::cout << imgcoord_v[p] << ", ";
    }
    std::cout << imgcoord_v[3] << std::endl;
    
    // define crops
    
    for ( size_t p=0; p<3; p++ ) {
      // define the bounds for this plane
      auto const& meta = adc_v[p].meta();
      int tbounds[2] = { (int)(imgcoord_v[3]-256*meta.pixel_height()),
                         (int)(imgcoord_v[3]+256*meta.pixel_height()) };
      int wbounds[2] = { (int)(imgcoord_v[p]-256*meta.pixel_width()),
                         (int)(imgcoord_v[p]+256*meta.pixel_width()) };
      
      // correct bounds
      int wireMax = 2400;
      if(p == 2) wireMax = meta.max_x(); 
      
      int vtx_row = imgcoord_v[3];
      int vtx_col = imgcoord_v[p];
      
      if ( tbounds[0]<=meta.min_y() ) {
        //vtx_row -= meta.min_y() + meta.pixel_height() - tbounds[0];
        tbounds[0] = meta.min_y()+meta.pixel_height();
        tbounds[1] = tbounds[0] + 512*meta.pixel_height();
      }
      if ( tbounds[1]>=meta.max_y() ) {
        //vtx_row += tbounds[1] - meta.max_y() + meta.pixel_height(); 
        tbounds[1] = meta.max_y() - meta.pixel_height();
        tbounds[0] = tbounds[1] - 512*meta.pixel_height();
          }
      if ( wbounds[0]<=meta.min_x() ) {
        //vtx_col -= meta.min_x() + meta.pixel_width() - wbounds[0];
        wbounds[0] = meta.min_x()+meta.pixel_width();
        wbounds[1] = wbounds[0] + 512*meta.pixel_width();
      }
      if ( wbounds[1]>=wireMax ) {
        //vtx_col += wbounds[1] - wireMax + meta.pixel_width();
        wbounds[1] = wireMax-meta.pixel_width();
        wbounds[0] = wbounds[1] - 512*meta.pixel_width();
      }
      
      // get origin coordinates
      vtx_row = 512 + meta.row(vtx_row) - meta.row(tbounds[0]);
      vtx_col = meta.col(vtx_col) - meta.col(wbounds[0]);

      std::cout << "Plane " << p << " crop VTX R/C: " << vtx_row << ", " << vtx_col << std::endl;
      
      //Channel status saving
      larcv::Image2D chstatusimg = GetChStatusImage(512,ev_wirestatus.Status(p),meta.col(wbounds[0]));

      // make image marking where vertex is
      //larcv::Image2D chstatusimg = GetChStatusVtxImage(512,ev_wirestatus.Status(p),meta.col(wbounds[0]),vtx_row,vtx_col);
      
      // define meta for crop
      larcv::ImageMeta cropmeta( 512*meta.pixel_width(), 512*meta.pixel_height(),
                                 512, 512, wbounds[0], tbounds[1],
                                 meta.plane() );
      
      larcv::Image2D crop   = adc_v[p].crop(cropmeta);
      larcv::Image2D sscrop = shower_img_v[p].crop(cropmeta);
      
      // mask the ADC image using ssnet
      for ( size_t r=0; r<crop.meta().rows(); r++ ) {
        for ( size_t c=0; c<crop.meta().cols(); c++ ) {
          if ( crop.pixel(r,c)<adc_threshold ) {
            crop.set_pixel(r,c,0.0);
          }
          if ( sscrop.pixel(r,c)<=showerscore_threshold ) {
            crop.set_pixel(r,c,0.0);
          }
        }
      }
      
      cropimg_v.emplace_back( std::move(crop) );
      statusimg_v.emplace_back( std::move(chstatusimg) );
      
    }//end of loop over planes
  }
  

  larcv::Image2D ShowerRecoUtil::GetChStatusImage(int L, const larcv::ChStatus& ch, int leftCol){
    larcv::Image2D chstatusimg(L,L);
    chstatusimg.paint(0);
    for(int col = 0; col < L; ++col){
      if(ch.Status(col+leftCol) != 4) chstatusimg.paint_col(col,1);
    }
    return chstatusimg;
  }
  
  larcv::Image2D ShowerRecoUtil::GetChStatusVtxImage(int L, const larcv::ChStatus& ch, int leftCol, int vtxRow, int vtxCol){
    larcv::Image2D chstatusimg(L,L);
    chstatusimg.paint(0);
    for(int col = 0; col < L; ++col){
      if(ch.Status(col+leftCol) != 4) chstatusimg.paint_col(col,1);
    }
    chstatusimg.set_pixel(vtxRow,vtxCol,10);
    return chstatusimg;
  }
  
}
}
