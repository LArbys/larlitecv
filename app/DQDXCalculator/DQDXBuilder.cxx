#ifndef DQDXBUILDER_CXX
#define DQDXBUILDER_CXX

#include "DQDXBuilder.h"


namespace larlitecv {
  DQDXBuilder::DQDXBuilder()
  : _img_v(3,larcv::Image2D(1008,3456))
  {
    Initialize_DQDXBuilder();
  }


  DQDXBuilder::DQDXBuilder(std::vector<larcv::Image2D> img_v)
  :  _img_v(img_v)
  {
    Initialize_DQDXBuilder();
  }

  void DQDXBuilder::Initialize_DQDXBuilder(){
    views[0] = larlite::geo::kU;
    views[1] = larlite::geo::kV;
    views[2] = larlite::geo::kZ;
    geo = larutil::Geometry::GetME();
    Initialize_Wirept_Vector();
    residual_dqdx = TH2D("residual_dqdx","residual_dqdx ",20,0,100,30,0,1500);
    residual_dqdx.SetXTitle("DQ/DX");
    residual_dqdx.SetYTitle("Distance From End (cm)");
  }


  void DQDXBuilder::Initialize_Wirept_Vector(){
    // Create plane Equation variables for each wire for each wire Plane
    // Equation ax+by+cz+d = 0
    wire_planar_equations.resize(3);
    for (int p =0; p < 3; p++){
      // std::cout << "Plane " << p << "\n";
      wire_planar_equations[p].reserve(4000);
      for (int w=0;w<(int)_img_v[p].meta().cols(); w++){
        double xyzstart[3];
        double xyzend[3];
        if ((p == 2) || (w < 2400)) {
          geo->WireEndPoints( (UChar_t)p, (UInt_t)w, xyzstart, xyzend );
          TVector3 tvec_start(xyzstart[0],xyzstart[1],xyzstart[2]);
          TVector3 tvec_end(xyzend[0],xyzend[1],xyzend[2]);
          double shifted_x = xyzend[0]+10;
          TVector3 pt_3(shifted_x,xyzend[1],xyzend[2]); //Defining a plane, but need a 3rd point
          std::vector<double> plane_vec = define_plane(tvec_start, tvec_end, pt_3);
          // Stuff plane definition into container
          wire_planar_equations[p].push_back(plane_vec);
          // std::cout <<w << " " << _img_v[p].meta().cols() <<" " <<  plane_vec[0] << " " << plane_vec[1] << " " << plane_vec[2] << " " << plane_vec[3] << "\n";
        }
      }
    }
  }

  TVector3 DQDXBuilder::get_halfway_plane_pt( int w1,  int w2, const TVector3 pt1, const TVector3 pt2, const int plane){
    /*
    This function takes in pt1 and 2 which form a line segment, and w1 and 2
    which give 2 wires. Then if finds the halfway pt between the two driftplanes
    defined by the wires along the line defined by the wiresegment
    */
    // Assure w1 <= w2
    if (w1 > w2){
      int tmp = w1;
      w1 = w2;
      w2 = tmp;
      std::cout << "    Flipping Wire Order\n";
    }
    if (w2 < 0){
      std::cout << "w2 cannot be below 0\n";
      throw  -1;
    }
    else if (w1 > (int)wire_planar_equations[plane].size()-1){
      std::cout << "w1 cannot be above " <<  wire_planar_equations[plane].size()-1 << "\n";
      throw  -1;
    }

    // Both w1 and w2 within bounds, should be almost all cases
    if ((w1 >= 0) && (w2 < (int)wire_planar_equations[plane].size())){


      TVector3 pt_w1 = get_line_plane_intersection_pt(pt1, pt2, wire_planar_equations[plane][w1],false);
      TVector3 pt_w2 = get_line_plane_intersection_pt(pt1, pt2, wire_planar_equations[plane][w2],false);


      double new_x = pt_w1.X() + 0.5*(pt_w2.X()-pt_w1.X());
      double new_y = pt_w1.Y() + 0.5*(pt_w2.Y()-pt_w1.Y());
      double new_z = pt_w1.Z() + 0.5*(pt_w2.Z()-pt_w1.Z());

      TVector3 halfway_pt(new_x,new_y,new_z);

      return halfway_pt;
    }
    else if (w1 < 0){
      throw -1;
      // Take distance between planes for w2 and w2+1 and reverse it
      w1 = w2+1;
      TVector3 pt_w1 = get_line_plane_intersection_pt(pt1, pt2, wire_planar_equations[plane][w1],false);
      TVector3 pt_w2 = get_line_plane_intersection_pt(pt1, pt2, wire_planar_equations[plane][w2],false);

      double new_x = pt_w2.X() - 0.5*(pt_w1.X()-pt_w2.X());
      double new_y = pt_w2.Y() - 0.5*(pt_w1.Y()-pt_w2.Y());
      double new_z = pt_w2.Z() - 0.5*(pt_w1.Z()-pt_w2.Z());

      TVector3 halfway_pt(new_x,new_y,new_z);
      return halfway_pt;
    }
    else{
      throw -1;
      // Take distance between planes for w1 and w1-1 and reverse it
      w2 = w1-1;
      TVector3 pt_w1 = get_line_plane_intersection_pt(pt1, pt2, wire_planar_equations[plane][w1],false);
      TVector3 pt_w2 = get_line_plane_intersection_pt(pt1, pt2, wire_planar_equations[plane][w2],false);

      double new_x = pt_w1.X() - 0.5*(pt_w2.X()-pt_w1.X());
      double new_y = pt_w1.Y() - 0.5*(pt_w2.Y()-pt_w1.Y());
      double new_z = pt_w1.Z() - 0.5*(pt_w2.Z()-pt_w1.Z());

      TVector3 halfway_pt(new_x,new_y,new_z);
      return halfway_pt;
    }
  }

  void DQDXBuilder::Add_Track_To_ResidualDQDX(const larlite::track& reco3d_track, const int plane){
    for (int i=0; i< reco3d_track.NumberTrajectoryPoints();i++){
      residual_dqdx.Fill(reco3d_track.Length(i),reco3d_track.DQdxAtPoint(i,views[plane]));
    }
  }

  void DQDXBuilder::Plot_Orig_New_Raw_Track_Crop(const larlite::track& orig_track, const larlite::track& new_track){
    // std::vector<double> dqdx_test;
    // dqdx_test.reserve(dq_y.size());
    // std::vector<double> s_dist_test_y;
    // s_dist_test_y.reserve(dq_y.size());
    // larcv::Image2D track_img(1008,3456);
    // int min_r = 99999;
    // int max_r = -1;
    // int min_c = 99999;
    // int max_c = -1;
    //
    //
    //
    // std::cout << "Creating Re-Reconstructed Track DQDX for plots\n";
    // std::cout << "index dq dx dqdx col row dist_from_end\n";
    // for (int i = 0;i <= (int)dq_y.size()-1; i++){
    //   // if (i==0) {continue;}
    //   double ddd = dq_y[i]/dx_y[i];
    //   double thresh = 9999999;
    //   int row = calc_row_from_x(new_steps_y[i].x(),_img_v[plane].meta());
    //   int col = geo->NearestWire(new_steps_y[i], plane);
    //   if (row < min_r) min_r=row;
    //   if (row > max_r) max_r=row;
    //   if (col < min_c) min_c=col;
    //   if (col > max_c) max_c=col;
    //   print_tvec(new_steps_y[i]);
    //   std::cout << i << "   "<< dq_y.at(i) << "    " <<
    //     dx_y.at(i) << "    " << ddd << "    " << col << "    "<<
    //     row  << "    " <<
    //     s_distance_y.at(dq_y.size()-1) - s_distance_y.at(i) <<  "\n";
    //   track_img.set_pixel(row,col,ddd);
    //   if (ddd > thresh) ddd=dqdx_test.at(dqdx_test.size()-1);
    //   dqdx_test.push_back(ddd);
    //   s_dist_test_y.push_back(s_distance_y.at(dq_y.size()-1) - s_distance_y.at(i));
    //   if (ddd > residual_dqdx.GetYaxis()->GetXmax()) {
    //     ddd = residual_dqdx.GetYaxis()->GetXmax()-1;
    //   }
    //   residual_dqdx.Fill(s_distance_y.at(dq_y.size()-1) - s_distance_y.at(i), ddd);
    // }
    // min_r -=2;
    // min_c -=2;
    // int num_r = max_r-min_r+2;
    // int num_c = max_c-min_c+2;
    // std::cout << "Make From Revamp\n";
    // larlitecv::make_dqdx_curve(dqdx_test, s_dist_test_y, "dqdx_mine_revamp");
    // larcv::Image2D origtrack_img(1008,3456);
    // for (int i=0;i<(int)reco3d_track.NumberTrajectoryPoints();i++){
    //   TVector3 this_pt = reco3d_track.LocationAtPoint(i);
    //   double dqdx_y = reco3d_track.DQdxAtPoint(i, views[2]);
    //   double wire = geo->WireCoordinate(this_pt, plane);
    //   double row  = calc_row_from_x(this_pt.x(), _img_v[2].meta());
    //   origtrack_img.set_pixel(row,wire,dqdx_y);
    //   // std::cout << wire << " " << row << "\n";
    //
    // }
    //
    // TH2D raw = TH2D("raw","raw ",num_c,min_c,(double)min_c+num_c,num_r,min_r,(double)min_r+num_r);
    // TH2D origtrack_h = TH2D("track_h","track_h ",num_c,min_c,(double)min_c+num_c,num_r,min_r,(double)min_r+num_r);
    // TH2D track_h = TH2D("track_h","track_h ",num_c,min_c,(double)min_c+num_c,num_r,min_r,(double)min_r+num_r);
    // for (int c = min_c; c <=min_c+num_c;c++){
    //   for (int r = min_r; r <=min_r+num_r;r++){
    //     track_h.SetBinContent(c-min_c,r-min_r, track_img.pixel(r,c) );
    //     raw.SetBinContent(c-min_c,r-min_r, _img_v[plane].pixel(r,c) );
    //     origtrack_h.SetBinContent(c-min_c,r-min_r, origtrack_img.pixel(r,c) );
    //   }
    // }
    // larlitecv::make_evdisp(raw, "raw_yplane", "Col (Wire)" , "Row (6xTick)");
    // larlitecv::make_evdisp(track_h, "trackdqdx_yplane");
    // larlitecv::make_evdisp(origtrack_h, "origtrackdqdx_yplane");
    // std::vector<double> dqdx_final(dq_y.size(), -1);
  }

  larlite::track DQDXBuilder::calc_dqdx_track_revamp(const larlite::track& reco3d_track, const int plane) {
    /*
    This is where the magic happens. This function is designed to take
    a reco3d track with 3d points in space, and then create a new track with
    a dqdx vector that isn't charge density, but dqdx at various points
    */
    larlite::track dqdx_track; //Start with empty track

    new_steps_y.clear();
    dq_y.clear();
    dx_y.clear();
    s_distance_y.clear();
    dqdx_track.set_track_id( reco3d_track.ID() );
    dqdx_track.reserve(800); //Reserve a lot of points just in case
    // Just doing a single plane, we'll have to run this function multiple times.
    new_steps_y.reserve(800);
    dq_y.reserve(800);
    dx_y.reserve(800);
    s_distance_y.reserve(800);
    double cumulative_s_distance = 0;
    TVector3 modified_this_pt;
    int i = 1;
    while (i< (int)reco3d_track.NumberTrajectoryPoints()-1){
    // for (int i =1; i<=(int)reco3d_track.NumberTrajectoryPoints()-1; i++){
      TVector3 this_pt = reco3d_track.LocationAtPoint(i-1);
      // If not the first pair, then set this_pt equal to end of last step's dx
      if (i != 1){
        this_pt = modified_this_pt;
      }
      TVector3 next_pt = reco3d_track.LocationAtPoint(i);

      double start_wire_d = geo->WireCoordinate(this_pt, plane);
      double end_wire_d   = geo->WireCoordinate(next_pt, plane);
      int start_wire_i  = -1;
      int end_wire_i    = -1;
      bool wire_is_increasing;
      start_wire_i = (int) std::round(start_wire_d); // round lower wire
      end_wire_i   = (int) std::round(end_wire_d); // round higher wire
      if (start_wire_d <= end_wire_d) {
        wire_is_increasing = true;
      }
      else{
        wire_is_increasing = false;
      }

      // Now we have the wires we are going to take points on for this Reco3D step
      // Each wire will constitute a step in our new track.
      double dq=0;
      double dx;

      int start_row;
      int end_row;
      TVector3 dx_pt1;
      TVector3 dx_pt2;
      TVector3 new_pt_3d;
      if (start_wire_i == end_wire_i){
        if ((this_pt.x() == next_pt.x()) && (this_pt.y() == next_pt.y()) && (this_pt.z() == next_pt.z())){
          i++;
          continue;
        }
        double this_next_distance = distance_between_pt(this_pt,next_pt);

        TVector3 unit_vec = get_unit_vector(this_pt,next_pt);
        // Handle this case, track doesn't jump wires, just a change in ticks
        // Sum over wire the range of ticks branchedd by the projection
        ///////////////////////////////////////////////////////////
        // WARNING NOT GRABBING WIDTH OF THESE PURELY VERTICAL TRACKSTEPS!
        ///////////////////////////////////////////////////////////
        dx_pt1 = get_halfway_plane_pt(start_wire_i-1,start_wire_i,this_pt, next_pt, plane);
        dx_pt2 = get_halfway_plane_pt(start_wire_i,start_wire_i+1,this_pt, next_pt, plane);
        dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
        new_pt_3d = get_line_plane_intersection_pt(dx_pt1, dx_pt2, wire_planar_equations[plane][start_wire_i]);

        if ((dx > 1.0) && (dx  > this_next_distance)){
          if (wire_is_increasing){
            dx_pt2 = next_pt;
            dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
            if ((dx > 1.0) && (dx  > this_next_distance)){
              dx_pt1.SetXYZ(dx_pt2.X() - unit_vec.X(),dx_pt2.Y() - unit_vec.Y(),dx_pt2.Z() - unit_vec.Z());
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
            }
            TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
            new_pt_3d = new_tmp;
          }
          else{
            // TVector3 new_dxpt1(dx_pt2.x()+unit_vec.x(), dx_pt2.y()+unit_vec.y(), dx_pt2.z()+unit_vec.z());
            dx_pt1 = next_pt;
            dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
            if ((dx > 1.0) && (dx  > this_next_distance)){
              dx_pt2.SetXYZ(dx_pt1.X() - unit_vec.X(),dx_pt1.Y() - unit_vec.Y(),dx_pt1.Z() - unit_vec.Z());
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
            }
            TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
            new_pt_3d = new_tmp;

          }
        }
        if (dx < 0.3){
          if (wire_is_increasing){
            dx_pt2.SetXYZ(dx_pt1.x()+unit_vec.x()*0.3, dx_pt1.y()+unit_vec.y()*0.3, dx_pt1.z()+unit_vec.z()*0.3);
            dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
            TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
            new_pt_3d = new_tmp;
          }
          else{
            dx_pt1.SetXYZ(dx_pt2.x()+unit_vec.x()*0.3, dx_pt2.y()+unit_vec.y()*0.3, dx_pt2.z()+unit_vec.z()*0.3);
            dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
            TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
            new_pt_3d = new_tmp;
          }
        }

        // int row_this = calc_row_from_x(new_pt_3d.x(),_img_v[plane].meta());
        // int col_this = geo->NearestWire(new_pt_3d, plane);
        //If step distance is more than 1.0 then lets put a point down early
        //and check for the next reco3d point

        cumulative_s_distance += dx;
        start_row    = calc_row_from_x(dx_pt1.x(),_img_v[plane].meta());
        end_row      = calc_row_from_x(dx_pt2.x(),_img_v[plane].meta());
        //Flip if ordered wrong
        if (start_row > end_row){
          int tmp_row = end_row;
          end_row = start_row;
          start_row = tmp_row;
        }
        dq += sum_charge_dq(start_row, end_row, start_wire_i, plane);
        new_steps_y.push_back(new_pt_3d);
        dq_y.push_back(dq);
        dx_y.push_back(dx);
        s_distance_y.push_back(cumulative_s_distance);
        //Set next reco3d step's this_pt to end of this dx measurement
        if (wire_is_increasing){
          modified_this_pt = dx_pt2 ;
        }
        else{
          modified_this_pt = dx_pt1 ;
        }
        dq=0; //reset dq for next step
      }
      else if (wire_is_increasing){
        //Increasing In Wire
        int col = start_wire_i;
        while (col <= end_wire_i){
          if ((this_pt.x() == next_pt.x()) && (this_pt.y() == next_pt.y()) && (this_pt.z() == next_pt.z())){
            col++;
            continue;
          }
          bool advance_col = true;
          double this_next_distance = distance_between_pt(this_pt,next_pt);
          dx_pt1 = get_halfway_plane_pt(col-1, col   , this_pt, next_pt, plane);
          dx_pt2 = get_halfway_plane_pt(col  , col+1 , this_pt, next_pt, plane);
          TVector3 unit_vec = get_unit_vector(this_pt,next_pt);

          dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
          new_pt_3d = get_line_plane_intersection_pt(dx_pt1, dx_pt2, wire_planar_equations[plane][col]);
          if ((dx > 1.0) && (dx  > this_next_distance)){
            if (wire_is_increasing){
              dx_pt2 = next_pt;
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              if ((dx > 1.0) && (dx  > this_next_distance)){
                dx_pt1.SetXYZ(dx_pt2.X() - unit_vec.X(),dx_pt2.Y() - unit_vec.Y(),dx_pt2.Z() - unit_vec.Z());
                dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              }
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
            else{
              dx_pt1 = next_pt;
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              if ((dx > 1.0) && (dx  > this_next_distance)){
                dx_pt2.SetXYZ(dx_pt1.X() - unit_vec.X(),dx_pt1.Y() - unit_vec.Y(),dx_pt1.Z() - unit_vec.Z());
                dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              }
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
          }
          if (dx < 0.3){
            if (wire_is_increasing){
              dx_pt2.SetXYZ(dx_pt1.x()+unit_vec.x()*0.3, dx_pt1.y()+unit_vec.y()*0.3, dx_pt1.z()+unit_vec.z()*0.3);
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
            else{
              dx_pt1.SetXYZ(dx_pt2.x()+unit_vec.x()*0.3, dx_pt2.y()+unit_vec.y()*0.3, dx_pt2.z()+unit_vec.z()*0.3);
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
          }

          // int row_this = calc_row_from_x(new_pt_3d.x(),_img_v[plane].meta());
          // int col_this = geo->NearestWire(new_pt_3d, plane);
          cumulative_s_distance += dx;
          start_row    = calc_row_from_x(dx_pt1.x(),_img_v[plane].meta());
          end_row      = calc_row_from_x(dx_pt2.x(),_img_v[plane].meta());
          //Flip if ordered wrong
          if (start_row > end_row){
            int tmp_row = end_row;
            end_row = start_row;
            start_row = tmp_row;
          }
          dq += sum_charge_dq(start_row, end_row, col, plane);
          new_steps_y.push_back(new_pt_3d);
          dq_y.push_back(dq);
          dx_y.push_back(dx);
          s_distance_y.push_back(cumulative_s_distance);
          //Set next reco3d step's this_pt to end of this dx measurement
          modified_this_pt = dx_pt2 ;
          this_pt = dx_pt2 ;
          dq=0; //reset dq for next step
          if (advance_col){col++;}
        }
      }
      else if  ( wire_is_increasing != true ) {
        // Decreasing in Wire
        for (int col = start_wire_i; col >= end_wire_i; col--){
          if ((this_pt.x() == next_pt.x()) && (this_pt.y() == next_pt.y()) && (this_pt.z() == next_pt.z())){
            continue;
          }
          double this_next_distance = distance_between_pt(this_pt,next_pt);
          dx_pt1 = get_halfway_plane_pt(col-1, col   , this_pt, next_pt, plane);
          dx_pt2 = get_halfway_plane_pt(col  , col+1 , this_pt, next_pt, plane);
          TVector3 unit_vec = get_unit_vector(this_pt,next_pt);



          dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
          new_pt_3d = get_line_plane_intersection_pt(dx_pt1, dx_pt2, wire_planar_equations[plane][col]);
          if ((dx > 1.0) && (dx  > this_next_distance)){
            if (wire_is_increasing){
              dx_pt2 = next_pt;
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              if ((dx > 1.0) && (dx  > this_next_distance)){
                dx_pt1.SetXYZ(dx_pt2.X() - unit_vec.X(),dx_pt2.Y() - unit_vec.Y(),dx_pt2.Z() - unit_vec.Z());
                dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              }
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
            else{
              dx_pt1 = next_pt;
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              if ((dx > 1.0) && (dx  > this_next_distance)){
                dx_pt2.SetXYZ(dx_pt1.X() - unit_vec.X(),dx_pt1.Y() - unit_vec.Y(),dx_pt1.Z() - unit_vec.Z());
                dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              }
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
          }
          if (dx < 0.3){
            if (wire_is_increasing){
              dx_pt2.SetXYZ(dx_pt1.x()+unit_vec.x()*0.3, dx_pt1.y()+unit_vec.y()*0.3, dx_pt1.z()+unit_vec.z()*0.3);
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
            else{
              dx_pt1.SetXYZ(dx_pt2.x()+unit_vec.x()*0.3, dx_pt2.y()+unit_vec.y()*0.3, dx_pt2.z()+unit_vec.z()*0.3);
              dx = larlitecv::distance_between_pt(dx_pt1,dx_pt2);
              TVector3 new_tmp((dx_pt1.x()+dx_pt2.x())/2, (dx_pt1.y()+dx_pt2.y())/2, (dx_pt1.z()+dx_pt2.z())/2);
              new_pt_3d = new_tmp;
            }
          }

          // int row_this = calc_row_from_x(new_pt_3d.x(),_img_v[plane].meta());
          // int col_this = geo->NearestWire(new_pt_3d, plane);


          cumulative_s_distance += dx;
          start_row    = calc_row_from_x(dx_pt1.x(),_img_v[plane].meta());
          end_row      = calc_row_from_x(dx_pt2.x(),_img_v[plane].meta());
          //Flip if ordered wrong
          if (start_row > end_row){
            int tmp_row = end_row;
            end_row = start_row;
            start_row = tmp_row;
          }
          dq += sum_charge_dq(start_row, end_row, col, plane);
          new_steps_y.push_back(new_pt_3d);
          dq_y.push_back(dq);
          dx_y.push_back(dx);
          s_distance_y.push_back(cumulative_s_distance);
          //Set next reco3d step's this_pt to end of this dx measurement
          modified_this_pt = dx_pt1 ;
          this_pt = dx_pt1 ;

          dq=0; //reset dq for next step
        }
      }
      else{
        std::cout << "This shouldn't happen, wire_is_increasing is either true or false\n";
        throw -1;
      }
      i++;
    }

    std::vector<double> dqdx_final;
    dqdx_final.reserve(new_steps_y.size());
    for (int idx=0;idx < dq_y.size(); idx++){
      TVector3 direction;
      if (idx != dq_y.size()-1){
        direction.SetXYZ(new_steps_y[idx+1].X()-new_steps_y[idx].X() ,
          new_steps_y[idx+1].Y()-new_steps_y[idx].Y(),
          new_steps_y[idx+1].Z()-new_steps_y[idx].Z());
      }
      else{
        // Set last point dir to 0,0,0
        direction.SetXYZ(0.,0.,0.);
      }
      dqdx_track.add_direction(direction);
      dqdx_final.push_back(dq_y[idx]/dx_y[idx]);
      dqdx_track.add_vertex(new_steps_y[idx]);
    }
    //Make all 3 from the y plane for now
    dqdx_track.add_dqdx(dqdx_final);
    dqdx_track.add_dqdx(dqdx_final);
    dqdx_track.add_dqdx(dqdx_final);

    return dqdx_track;
  }

  void DQDXBuilder::test_wirecoord_driftplanes(larlite::track reco3d_track){
    int plane=2;
    // Try comparing the Wire of a 3D point for 2 methods:
    double wire_idx = 2260;
    double wire_start[3];
    double wire_end[3];
    geo->WireEndPoints( (UChar_t)plane, (UInt_t)wire_idx, wire_start, wire_end );
    std::cout << "Wire Start  " << wire_start[0] << " "<<wire_start[1]<<" "<<wire_start[2]<<"\n";
    std::cout << "Wire End  " << wire_end[0] << " "<<wire_end[1]<<" "<<wire_end[2]<<"\n";
    // Make Wire Plane by shifting x coord of end and find plane for those 3 points
    TVector3 wire_start_tvec(wire_start[0], wire_start[1],wire_start[2]);
    TVector3 wire_end_tvec(wire_end[0], wire_end[1],wire_end[2]);
    TVector3 wire_shift_tvec(wire_end[0]+10., wire_end[1],wire_end[2]);
    // Find Vectors from start to end and start to shift
    std::vector<double> vec1 =  {wire_end_tvec.x() - wire_start_tvec.x(), wire_end_tvec.y()-wire_start_tvec.y(), wire_end_tvec.z()-wire_start_tvec.z()};
    std::vector<double> vec2 =  {wire_shift_tvec.x() - wire_start_tvec.x(), wire_shift_tvec.y()-wire_start_tvec.y(), wire_shift_tvec.z()-wire_start_tvec.z()};
    // Cross product components for normal vector, use that and one of the points to get plane equation
    double cross_i = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    double cross_j = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    double cross_k = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    double a = cross_i;
    double b = cross_j;
    double c = cross_k;
    double d = cross_i*(0-wire_start_tvec.x()) + cross_j*(0-wire_start_tvec.y()) + cross_k*(0-wire_start_tvec.z()) ;
    std::vector<double> plane_vec = {a,b,c,d};
    // Now we have a plane for the wire_idx, let's get a point on it:
    TVector3 test_pt(50.,0.,678.25) ; // This is the y plane, so x and y ppoints don't matter
    // Check it's in plane:
    double check_val = a*test_pt.X() + b*test_pt.Y() + c*test_pt.Z() + d;
    std::cout << "check_val " << check_val << "\n"; //If in plane this is 0
    // Now let's check NearestWire for the comparison:
    double wire_exact = geo->WireCoordinate(test_pt, 2);
    std::cout << "wire_exact " << wire_exact << "\n"; //Should be 2260
  }


  int DQDXBuilder::calc_row_from_x(double x1, const larcv::ImageMeta& meta){
    int row;
    float tick = x1/(larutil::LArProperties::GetME()->DriftVelocity()*larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3) + 3200.0;
    if ( tick<meta.min_y() ) {
      row = -1;
	  }
	  else if ( tick>meta.max_y() ) {
      row = -1;
	  }
	  else {
	    // within the image
	    row = meta.row( tick );
	  }
    return row;
  }


  double DQDXBuilder::sum_charge_dq(const int low_row, const int high_row, const int col, const int plane, const int buffer){
    double tot_dq = 0;
    for (int r_idx = low_row-buffer; r_idx < high_row+1+buffer; r_idx++){
      if ((r_idx > 1008) || (r_idx < 0)){
        continue;
      }
      double this_dq = _img_v[plane].pixel(r_idx, col);
      // std::cout << "    R C DQ  " << r_idx << " " << col << "  " << this_dq << "\n";
      tot_dq += this_dq;
    }
    return tot_dq;
  }
}
#endif
