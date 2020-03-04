#include "SecondShower.h"

#include "TMath.h"

#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPGraph.h"
#include "DataFormat/track.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {
namespace ssnetshowerreco {

  /**
   * constructor.
   *
   * sets some defaults for parameters
   */
//------------------------------------------------------------------------------
//copy of  sparse enclosed charge function for use in second shower search
std::vector< std::vector<int> >  SecondShower::_enclosedCharge( std::vector<std::vector<float>> chargeMap,
                                                                   float theta,
                                                                   float& sumIn,
                                                                   std::vector< std::vector<float> >& triangle,
                                                                   bool calcE,
                                                                   int vtx_col,
                                                                   int vtx_row,
                                                                   float shLen,
                                                                   float shOpen) {

  std::vector< std::vector<int> > pix_v;

  sumIn = 0.;
  triangle.resize(3);
  for (int i=0; i<3; i++ ) {
    triangle[i].resize(2,0);
  }

  std::vector<int> vtx = { vtx_col, vtx_row };
  float t1X = vtx[0];
  float t1Y = vtx[1];

  float t2X = vtx[0] + shLen*cos(theta+shOpen);
  float t2Y = vtx[1] + shLen*sin(theta+shOpen);

  float t3X = vtx[0] + shLen*cos(theta-shOpen);
  float t3Y = vtx[1] + shLen*sin(theta-shOpen);

  triangle[0] = { t1X, t1Y };
  triangle[1] = { t2X, t2Y };
  triangle[2] = { t3X, t3Y };

  if ( shLen==0.0 ) {
    // weird things when we have zero area triangle
    bool hasCharge = false;
    // std::cout<<"------"<<vtx_row<<","<<vtx_col<<"------------"<<std::endl;
    for (int ii =0;ii<chargeMap.size();ii++){
      // std::cout<<chargeMap[ii][0]<<","<<chargeMap[ii][1]<<std::endl;
      if(chargeMap[ii][0]==vtx_row&&chargeMap[ii][1]==vtx_col){
        sumIn = chargeMap[ii][2];
        hasCharge =true;
      }
    }
    //default to fall back on
    if (!hasCharge) sumIn = 10.0;
    return pix_v;
  }

  sumIn = 0;
  for (int ii =0;ii<chargeMap.size();ii++){
    if ( _isInside2(t1X,t1Y,t2X,t2Y,t3X,t3Y,chargeMap[ii][1],chargeMap[ii][0])){
      if (calcE) sumIn+=chargeMap[ii][2];
      else sumIn+=1;
      std::vector<int> pix = { (int)chargeMap[ii][0], (int)chargeMap[ii][1] };
      pix_v.push_back( pix );
    }
  }
  return pix_v;
}//end of function

//------------------------------------------------------------------------------

//copy of function to check if inside triangle
bool SecondShower::_isInside2( float x1, float y1,
                                  float x2, float y2,
                                  float x3, float y3,
                                  float x, float y ) {

  float d1, d2, d3;
  bool has_neg, has_pos;

  d1 = _sign( x, y, x1, y1, x2, y2 );
  d2 = _sign( x, y, x2, y2, x3, y3 );
  d3 = _sign( x, y, x3, y3, x1, y1 );

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  return !(has_neg && has_pos);

}//end of function

//------------------------------------------------------------------------------

//copy of sign function for use here
float SecondShower::_sign( float x1, float y1,
                              float x2, float y2,
                              float x3, float y3 ) {
  // (x1,y1) is the test point
  return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
}//end of function

//------------------------------------------------------------------------------

//function to choose optimal gap size
int SecondShower::ChooseGapSize(int vtx_col,int vtx_row,float shangle,
                std::vector<std::vector<float>> crop_v, int showernum){

  int nsteps = 1;
  //need to optimize this
  if (showernum == 1) nsteps = 50;
  else if (showernum == 2) nsteps = 300;
  std::vector<float> enclosedCharge_v;

  for (int step = 0; step < nsteps; step++){
    float xpos = vtx_col + (step * cos(shangle));
    float ypos = vtx_row + (step * sin(shangle));
    float tmpQ;
    triangle_t tmptri;
    _enclosedCharge( crop_v, shangle, tmpQ, tmptri,false, xpos, ypos );
    enclosedCharge_v.push_back(tmpQ);
  }
  //get max -rough for now
  float tmp_maxcharge = -1;
  int max_step = 0;
  for (int ii = 0; ii<enclosedCharge_v.size(); ii++){
    if (enclosedCharge_v[ii]>tmp_maxcharge){
      tmp_maxcharge = enclosedCharge_v[ii];
      max_step = ii;
    }
  }

  return max_step;

}//end of functions

//------------------------------------------------------------------------------

std::vector<double> SecondShower::SaveTrueEnergies(larlite::event_mcshower* ev_mcshower,
            std::vector<double> _scex){

  //retrieve and save true shower energies in order of closeness to vertex
  std::vector<double> energy_v(2,0);

  std::vector< std::vector<double>> showerstart_v;
  for (int shower = 0; shower<ev_mcshower->size();shower++){
    larlite::mcstep start = ev_mcshower->at(shower).DetProfile();
    std::vector<double> showerstart ={start.X(),start.Y(),start.Z()};
    showerstart_v.push_back(showerstart);
  }

  std::vector<double> vtx3d_true = { _scex[0], _scex[1], _scex[2]};
  // std::cout << " True Vertex Pos (" << vtx3d_true[0] << "," << vtx3d_true[1] << "," << vtx3d_true[2] << ")" << std::endl;
  std::vector<double> showerdist_v(2,0);
  if (showerstart_v[0][0]<10000&&showerstart_v[0][1]<10000&&showerstart_v[0][2]<10000){
    showerdist_v[0] = std::sqrt(std::pow(showerstart_v[0][0]-vtx3d_true[0],2)+std::pow(showerstart_v[0][1]-vtx3d_true[1],2)+std::pow(showerstart_v[0][2]-vtx3d_true[2],2));
  }
  else showerdist_v[0]= -1.0;
  if (showerstart_v[1][0]<10000&&showerstart_v[1][1]<10000&&showerstart_v[1][2]<10000){
    showerdist_v[1] = std::sqrt(std::pow(showerstart_v[1][0]-vtx3d_true[0],2)+std::pow(showerstart_v[1][1]-vtx3d_true[1],2)+std::pow(showerstart_v[1][2]-vtx3d_true[2],2));
  }
  else showerdist_v[1]= -1.0;
  // std::cout<<"Shower distances: "<<showerdist_v[0]<< " , "<<showerdist_v[1]<<std::endl;
  //shower energy in order of closeness
  double energy0 = ev_mcshower->at(0).DetProfile().E();
  double energy1 = ev_mcshower->at(1).DetProfile().E();

  int close_shower_id;
  if (showerdist_v[1]>=showerdist_v[0]&&showerdist_v[0]>=0){
    energy_v[0] = energy0;
    energy_v[1] = energy1;
    close_shower_id = 0;
  }
  else if (showerdist_v[0]>=showerdist_v[1]&&showerdist_v[1]>=0){
    energy_v[0] = energy1;
    energy_v[1] = energy0;
    close_shower_id = 1;
  }
  else if(showerdist_v[0]==-1&&showerdist_v[1]==-1){
    energy_v[0] = energy0;
    energy_v[1] = energy1;
    close_shower_id = 0;
  }
  else if (showerdist_v[0]==-1){
    energy_v[0] = energy1;
    energy_v[1] = energy0;
    close_shower_id = 1;
  }
  else if (showerdist_v[1]==-1){
    energy_v[1] = energy1;
    energy_v[0] = energy0;
    close_shower_id = 1;
  }

  return energy_v;
}//end of function

//------------------------------------------------------------------------------

std::vector<std::vector<double>> SecondShower::SaveTrueStarts(larlite::event_mcshower* ev_mcshower){

  //retrieve and save true shower start points in order of closeness to vertex
  std::vector<std::vector<double>> starts_v(ev_mcshower->size(),std::vector<double> (3,-1));

  for (int shower = 0; shower<ev_mcshower->size();shower++){
    larlite::mcstep start = ev_mcshower->at(shower).Start();
    std::vector<double> showerstart ={start.X(),start.Y(),start.Z()};
    starts_v[shower] = showerstart;
  }

  return starts_v;
}//end of function

//------------------------------------------------------------------------------

std::vector<std::vector<float>> SecondShower::Match_3D(std::vector<std::vector<float>> triangle_vv,
    std::vector<std::vector<std::vector<float>>> img_vvv,int shower_num){

  std::vector<std::vector<float>> match_vv(2,std::vector<float>(2,-1));
  std::vector<float> y_points_v = triangle_vv[3+shower_num];

  std::vector<float> ytimes_v;
  //loop through adc sparse image

  for (int ypix = 0; ypix<int(img_vvv[2].size());ypix++){
    //is it inside the shower?
    bool yinside = _isInside2( y_points_v[0], y_points_v[1],y_points_v[2], y_points_v[3],
                                      y_points_v[4], y_points_v[5],img_vvv[2][ypix][1], img_vvv[2][ypix][0]);
    if (yinside) ytimes_v.push_back(img_vvv[2][ypix][0]);
  }//end of loop through y sparse

  //loop through u and v planes
  for(int p = 0; p<2; p++){
    int sum_match_1 = 0;
    int sum_match_2 = 0;
    int tot_in_shower1 = 0;
    int tot_in_shower2 = 0;
    //loop through pixels in plane
    for(int pix = 0; pix<int(img_vvv[p].size()); pix++){
      //is it inside shower one?
      bool inside_shower1 = _isInside2(triangle_vv[2*p][0],triangle_vv[2*p][1],triangle_vv[2*p][2],
                                        triangle_vv[2*p][3],triangle_vv[2*p][4],triangle_vv[2*p][5],
                                        img_vvv[p][pix][1],img_vvv[p][pix][0]);
      if (inside_shower1){
        tot_in_shower1++;
        bool addone = false;
        float time = img_vvv[p][pix][0];
        //loop through times in yshower
        for (int ytime = 0; ytime<(int)ytimes_v.size();ytime++){
          if (ytimes_v[ytime] == time) addone =true;
        }
        if (addone) sum_match_1++;
      }
      //is is inside shower two?
      bool inside_shower2 = _isInside2(triangle_vv[2*p+1][0],triangle_vv[2*p+1][1],triangle_vv[2*p+1][2],
                                        triangle_vv[2*p+1][3],triangle_vv[2*p+1][4],triangle_vv[2*p+1][5],
                                        img_vvv[p][pix][1],img_vvv[p][pix][0]);
      if (inside_shower2){
        tot_in_shower2++;
        bool addone = false;
        float time = img_vvv[p][pix][0];
        //loop through times in yshower
        for (int ytime = 0; ytime<(int)ytimes_v.size();ytime++){
          if (ytimes_v[ytime] == time) addone =true;
        }
        if (addone) sum_match_2++;
      }//end of if inside
    }//end of pixel loop
    //calculate fraction match in time
    float frac1 = 0;
    if(tot_in_shower1>0) frac1 = (float)sum_match_1/(float)tot_in_shower1;
    float frac2 = 0;
    if(tot_in_shower2>0) frac2 = (float)sum_match_2/(float)tot_in_shower2;
    match_vv[p][0] = frac1;
    match_vv[p][1] = frac2;
  }//end of plane loop

  return match_vv;
}//end of function

//------------------------------------------------------------------------------

std::vector<std::vector<int>> SecondShower::ChooseBestMatch(std::vector<std::vector<float>> match_vv){
  //function to choose the highest fraction. Returns 1 in the index of any good,
  //returns 2 in index of best.
  std::vector<std::vector<int>> best_vv(2,std::vector<int>(2,0));
  float cutoff = .5;

  for(int plane = 0;plane<match_vv.size(); plane++){
    for (int shower = 0; shower<match_vv[plane].size();shower++){
        if (match_vv[plane][shower] > cutoff) best_vv[plane][shower] = 1;
    }
  }

  if (match_vv[0][0]>match_vv[0][1]&& match_vv[0][0]>match_vv[1][0]&& match_vv[0][0] >match_vv[1][1] &&match_vv[0][0] > cutoff)
    best_vv[0][0] = 2;
  else if (match_vv[0][1]>match_vv[1][0]&& match_vv[0][1] >match_vv[1][1]&&match_vv[0][1] > cutoff)
    best_vv[0][1] = 2;
  else if (match_vv[1][0] >match_vv[1][1]&&match_vv[1][0] > cutoff) best_vv[1][0] =2;
  else if (match_vv[1][1] > 0 &&match_vv[1][1] > cutoff) best_vv[1][1] =2;

  return best_vv;
}
//------------------------------------------------------------------------------

std::vector<std::vector<float>> SecondShower::SaveTrueProfile(int plane,
      larcv::EventImage2D* ev_segment, larcv::EventImage2D* ev_instance,
      larlite::event_mcshower* ev_mcshower, std::vector<std::vector<std::vector<float>>> img_vvv){

  //function to return the profile of a true shower
  std::vector<std::vector<float>> profile_vv;
  Utils Utils;
  //get shower start and 2d direction

  larlite::mcstep start = ev_mcshower->at(0).DetProfile();
  std::vector<double> showerstart3D ={start.X(),start.Y(),start.Z()};
  std::vector<int> showerstart2D = Utils.getProjectedPixel(showerstart3D, ev_segment->Image2DArray()[2].meta(), 3);

  std::vector<double> direction3D = {ev_mcshower->at(0).StartDir().X(),ev_mcshower->at(0).StartDir().Y(),ev_mcshower->at(0).StartDir().Z()};
  std::vector<double> secondpoint3D = {direction3D[0]*10+showerstart3D[0],direction3D[1]*10+showerstart3D[1],direction3D[2]*10+showerstart3D[2]};
  std::vector<int> secondpoint2D = Utils.getProjectedPixel(secondpoint3D, ev_segment->Image2DArray()[2].meta(), 3);

  float shangle =0;
  float x2 = float(secondpoint2D[plane+1])*.3;
  float x1 = float(showerstart2D[plane+1])*.3;
  float y2 = float(secondpoint2D[0])*6*.5*.1098;
  float y1 = float(showerstart2D[0])*6*.5*.1098;
  float a = (y2-y1)/(x2-x1);
  float b = y1- a*x1;

  std::vector<std::vector<float>> img_sparse = img_vvv[plane];
  for (int pixid =0;pixid<(int)img_sparse.size();pixid++){
    float ynew = float(img_sparse[pixid][0])*6*.5*.1098;
    float xnew = float(img_sparse[pixid][1])*.3;
    float segpix = ev_segment->Image2DArray()[plane].pixel(ynew,xnew);
    float instpix = ev_instance->Image2DArray()[plane].pixel(ynew,xnew);

    if (instpix>0&&segpix==3){
      float r = abs((y2-y1)*xnew - (x2-x1)*ynew +x2*y1+y2*x1)/sqrt(std::pow(y2-y1,2)+std::pow(x2-x1,2));

      float xint = (b*(b*xnew-a*ynew))/(a*a + b*b);
      float yint = (a*(-b*xnew+a*ynew))/(a*a +b*b);
      float d = sqrt(xint*xint + yint*yint);

      profile_vv.push_back({d,r});

    }
  }//end of loop through sparse image

  return profile_vv;
}//end of function
//------------------------------------------------------------------------------

std::vector<float> SecondShower::TruthMatchNCPi0(std::vector<float> triangle_v,
        std::vector<std::vector<float>> img_vv, larcv::EventImage2D* ev_segment,
        larcv::EventImage2D* ev_instance, int plane){


  //function to see if a ncpi0 shower is "good"
  //(purity,efficiency)
  std::vector<float> metrics(2,0);
  float total_gamma = 0.0;
  float efficiency =0.0;
  float fraction_good = 0.0;
  float total_in = 0.0;
  float total_good = 0.0;

  //loop through sparse image, see if in shower
  for (int pix = 0; pix<(int)img_vv.size();pix++){
    float segpix = ev_segment->Image2DArray()[plane].pixel(img_vv[pix][0],img_vv[pix][1]);
    float instpix = ev_instance->Image2DArray()[plane].pixel(img_vv[pix][0],img_vv[pix][1]);

    if (segpix ==4 && instpix > 0 ) total_gamma++;
    bool inside_shower = _isInside2(triangle_v[0],triangle_v[1],triangle_v[2],
                                      triangle_v[3],triangle_v[4],triangle_v[5],
                                      img_vv[pix][1],img_vv[pix][0]);
    //if good, get segment and instance value
    if (inside_shower){
      total_in++;
      //check if true gamma pixel
      if (segpix ==4 && instpix > 0) total_good++;
    }
  }

  if (total_in > 0.0) fraction_good = (float)total_good/(float)total_in;
  if (total_gamma >0.0) efficiency = (float)total_good/(float)total_gamma;

  metrics[0]=fraction_good;
  metrics[1]=efficiency;

  return metrics;

}//end of function
//------------------------------------------------------------------------------


std::vector<float> SecondShower::TruthMatchNueint(std::vector<float> triangle_v,
        std::vector<std::vector<float>> img_vv, larcv::EventImage2D* ev_segment,
        larcv::EventImage2D* ev_instance, int plane){


  //function to see if a ncpi0 shower is "good"
  //(purity,efficiency)
  std::vector<float> metrics(2,0);
  float total_electron = 0.0;
  float efficiency =0.0;
  float fraction_good = 0.0;
  float total_in = 0.0;
  float total_good = 0.0;


  //loop through sparse image, see if in shower
  for (int pix = 0; pix<(int)img_vv.size();pix++){
    float segpix = ev_segment->Image2DArray()[plane].pixel(img_vv[pix][0],img_vv[pix][1]);
    float instpix = ev_instance->Image2DArray()[plane].pixel(img_vv[pix][0],img_vv[pix][1]);

    if (segpix ==3 && instpix > 0) total_electron++;
    bool inside_shower = _isInside2(triangle_v[0],triangle_v[1],triangle_v[2],
                                      triangle_v[3],triangle_v[4],triangle_v[5],
                                      img_vv[pix][1],img_vv[pix][0]);
    //if good, get segment and instance value
    if (inside_shower){
      total_in++;
      //check if true gamma pixel
      if (segpix ==3 && instpix > 0) total_good++;
    }
  }
  // std::cout<<"total good: "<<total_good<<std::endl;
  // std::cout<<"total_in: "<<total_in<<std::endl;
  // std::cout<<"total_in: "<<total_electron<<std::endl;

  if (total_in > 0.0) fraction_good = (float)total_good/(float)total_in;
  if (total_electron >0.0) efficiency = (float)total_good/(float)total_electron;

  metrics[0]=fraction_good;
  metrics[1]=efficiency;

  return metrics;

}//end of function
//------------------------------------------------------------------------------

bool SecondShower::RunMassCalc(std::vector<std::vector<float>> _match_y1_vv,
      std::vector<std::vector<float>>_match_y2_vv,
      std::vector< std::vector<float> > _shower_energy_vv,
      std::vector< std::vector<float> >_secondshower_energy_vv){
  //function that runs checks to determined if this event should be used for
  // determining a pi0 mass peak. 1) are both reco showers >35MeV? 2) for each shower
  // is there an overlap >.2 3) do they match to different showers
  bool useevent = false;
  //check shower energies
  bool goodenergy = false;
  for( int vtx = 0;vtx<_shower_energy_vv.size(); vtx++){
    // std::cout<<"energies: "<< _shower_energy_vv.at(vtx)[2]<<" "<<_secondshower_energy_vv.at(vtx)[2]<<std::endl;
    if (_shower_energy_vv.at(vtx)[2] > 35 &&_secondshower_energy_vv.at(vtx)[2] > 35){
      goodenergy = true;
    }
  }
  // if (goodenergy) std::cout<<"Passes energy cut"<<std::endl;
  //check overlap
  bool hasoverlap1= false;
  bool hasoverlap2= false;

  for(int plane = 0;plane<_match_y1_vv.size(); plane++){
    for (int shower = 0; shower<_match_y1_vv[plane].size();shower++){
      if (_match_y1_vv[plane][shower]> .5){
        hasoverlap1 = true;
      }
      if (_match_y2_vv[plane][shower]> .5){
        hasoverlap2 = true;
      }
    }//end of showerloop
  }//end of plane loop


  //use this?
  if(goodenergy&& hasoverlap1 && hasoverlap2) useevent=true;
  return useevent;
}//end of function

//------------------------------------------------------------------------------

std::vector<std::vector<float>> SecondShower::Get3DPoints(std::vector<float> shower_points1,
        std::vector<float> shower_points2, std::vector<std::vector<std::vector<float>>> masked_adc_vvv,
        int planeofmatch, larcv::ImageMeta wire_meta){

  //function that returns the 3d positions for all points that overlap between the two showers
  std::vector<std::vector<float>> position_3d_v;
  std::vector<std::vector<float>> position_2d_v;
  std::vector<std::vector<float>> yinside_v;
  //loop through adc sparse image

  for (int ypix = 0; ypix<int(masked_adc_vvv[2].size());ypix++){
    //is it inside the shower?
    bool yinside = _isInside2( shower_points2[0], shower_points2[1], shower_points2[2],
            shower_points2[3], shower_points2[4], shower_points2[5], masked_adc_vvv[2][ypix][1],
            masked_adc_vvv[2][ypix][0]);
    // std::cout<<masked_adc_vvv[2][ypix][1]<<std::endl;
    if (yinside) yinside_v.push_back(masked_adc_vvv[2][ypix]);
  }//end of loop through y sparse

  //loop through other plane
  for (int pix = 0; pix<int(masked_adc_vvv[planeofmatch].size());pix++){
    bool isinside = _isInside2( shower_points1[0], shower_points1[1], shower_points1[2],
            shower_points1[3], shower_points1[4], shower_points1[5], masked_adc_vvv[planeofmatch][pix][1],
            masked_adc_vvv[planeofmatch][pix][0]);

    int ypixmatch = -1;
    for (int ypix = 0; ypix<int(yinside_v.size());ypix++){
      if (yinside_v[ypix][0] == masked_adc_vvv[planeofmatch][pix][0]) ypixmatch = ypix;
    }//end of loop through ypoints
    if (isinside && ypixmatch >= 0){
      //get 2d position - things I'll need for intersection point
      std::vector<float> position2d =  {yinside_v[ypixmatch][1], masked_adc_vvv[planeofmatch][pix][1], yinside_v[ypixmatch][0]};
      position_2d_v.push_back(position2d);
    }
  }//end of loop through other points

  //loop through 2d points.
  for(int pt = 0;pt<position_2d_v.size();pt++){
    double x,y,z;
    double tick = wire_meta.pos_y(position_2d_v[pt][2]);
    x = (tick -3200)*larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    larutil::Geometry::GetME()->IntersectionPoint( position_2d_v[pt][0], position_2d_v[pt][1], (UChar_t)2, (UChar_t)planeofmatch, y, z );
    // std::cout<<"x: "<<x<<", y: "<<y<<", z:"<<z<<std::endl;
    std::vector<float> position_3d = {(float)x,(float)y,(float)z};
    position_3d_v.push_back(position_3d);
  }

  return position_3d_v;
}//end of function

//------------------------------------------------------------------------------
float SecondShower::GetOpeningAngle(std::vector<std::vector<float>> firstshower,
        std::vector<std::vector<float>> secondshower){
  //function to get the opening angle of 2 clusters of 3d points
  float alpha = -1;
  //first calculate direction of showers
  std::vector<float> firstdirection = GetDirection(firstshower);
  std::vector<float> seconddirection = GetDirection(secondshower);

  //cos(alpha) = (u dot v)/(mag(u)*mag(v))
  float dotprod = (firstdirection[0]*seconddirection[0])+(firstdirection[1]*
                    seconddirection[1])+(firstdirection[2]*seconddirection[2]);

  //mag should be 1 because eigen vectors

  alpha = acos(dotprod);

  return alpha;
}//end of function
//------------------------------------------------------------------------------
std::vector<float> SecondShower::GetDirection(std::vector<std::vector<float>> shower){
  //function to calculate direction of 3d points using PCA
  std::vector<float> direction_v (3,0); //x,y,z
  TPrincipal* principal = new TPrincipal(3,"ND");
  //add data points to root object
  for (int ii = 0; ii<shower.size(); ii++){
    Double_t* data = new Double_t[3];
    for (int iii = 0; iii<shower[ii].size();iii++){
      data[iii] = (double) shower[ii][iii];
    }
    principal->AddRow(data);
  }
  //run PCA
  principal->MakePrincipals();
  // principal->Print("MSEV");
  const TMatrixD* eigenvectors = principal->GetEigenVectors();

  direction_v[0] = (float) (*eigenvectors)(0,0);
  direction_v[1] = (float) (*eigenvectors)(1,0);
  direction_v[2] = (float) (*eigenvectors)(2,0);

  return direction_v;
}//end of function
//------------------------------------------------------------------------------



}//end of ssnetshowerreco namespace
}//end of larlitecv namespace
