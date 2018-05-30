#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "Constant.h"
#include "Data.h"
#include "Douglas_ADI.h"
#include "Intersections.h"
#include "LOD.h"
#include "Mesh.h"
#include "New_Data_Type.h"
#include "TS.h"
// User-defined surfaces
#include "Surface_Cone.h"
#include "Surface_Cube.h"
#include "Surface_Cylinder.h"
#include "Surface_Dupin_Cyclide.h"
#include "Surface_Ellipsoid.h"
#include "Surface_Heart.h"
#include "Surface_Molecular.h"
#include "Surface_Pile.h"
#include "Surface_Tanglecube.h"
#include "Surface_Torus.h"
// User-defined diffusion coefficient beta
#include "Beta_0.h"
#include "Beta_1.h"
#include "Beta_2.h"
#include "Beta_3.h"
// Users-defined functions
#include "Eq_0.h"
#include "Eq_1.h"
#include "Eq_2.h"
#include "Eq_3.h"
#include "Eq_4.h"
#include "Eq_5.h"
#include "Eq_6.h"
#include "Eq_7.h"

using namespace std;

void ADI_Starting(Int_I, Int_I, Beta &, Mesh &, Intersections &, VecDoub_I &,
                  CubicDoub_I &, ofstream &);
void ADI_Solver(Int_I, Int_I, Beta &, Mesh &, Intersections &, Douglas_ADI &,
                VecDoub_I &, CubicDoub_I &, ofstream &);
void LOD_Starting(Int_I, Beta &, Mesh &, Intersections &, VecDoub_I &,
                  CubicDoub_I &, char &, ofstream &);
void LOD_Solver(Int_I, Beta &, Mesh &, Intersections &, LOD &, VecDoub_I &,
                CubicDoub_I &, char &, ofstream &);
void TS_Starting(Int_I, Beta &, Mesh &, Intersections &, VecDoub_I &,
                 CubicDoub_I &, ofstream &);
void TS_Solver(Int_I, Beta &, Mesh &, Intersections &, TS &, VecDoub_I &,
               CubicDoub_I &, ofstream &);
void Write_txt(Intersections &, Mesh &, Beta &, CubicDoub_I &, Int_I, Doub_I,
               Int_I);

int main(int argc, char *argv[]) {
  ofstream out_file;
  string file_name, current_time, out_file_name;
  VecDoub domain, running_time;
  CubicDoub uh;
  VecInt size;
  Surface_Cartesian *ex_ptr;
  Beta *beta_ptr;
  double t_begin, t_end, t;
  char surface, method;
  int equation, beta_code, accuracy, mib_method;

  // Arguments for a Cube, choose domain [-2,2;-2,2;-2,2], 'C'
  VecDoub arg_cube = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  // Arguments for a Cylinder, choose domain [-4,4;-4,4;-4,4], 'L'
  VecDoub arg_cylinder = {2.0, -3.0, 3.0};
  // Arguments for a Cone, choose domain [-4,4;-4,4;-4,1], 'O'
  VecDoub arg_cone = {1.0, 0.0, 0.0, -3.0, 0.0};
  // Arguments for an Ellipsoid with parameter[2,1,1,1], choose domain
  // [-4,4;-2,2;-2,2], 'E'
  VecDoub arg_ellipsoid = {2.0, 1.0, 1.0, 1.0};
  // Arguments for a Torus with parameter[0.5,0.2], choose domain
  // [-4,4;-4,4;-4,4], 'R'
  VecDoub arg_torus = {1.5, 0.7};
  // Arguments for a Dupin Cyclide with parameter[1,0.3,0.6], choose domain
  // [-5,5;-5,5;-5,5], 'D'
  VecDoub arg_dupin_cyclide = {1.0, 0.3, 0.6};
  // Arguments for a Heart Surface with parameter [1,9/4,1,-1,-1,-9/80], choose
  // domain [-4,4;-2,2;-4,4], 'H'
  VecDoub arg_heart = {1.0, 9.0 / 4.0, 1.0, -1.0, -1.0, -9.0 / 80.0};
  // Arguments for Tanglecube with parameter [1,-5,1,-5,1,-5,10], choose domain
  // [-3.99,3.99;-3.99,3.99;-4.99,4.99], 'T'
  VecDoub arg_tanglecube = {1.0, -5.0, 1.0, -5.0, 1.0, -5.0, 10.0};

  time_t rawtime;
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  current_time = asctime(timeinfo);

  vector<string> files = {"data/data1.txt"};
	// vector<string> files = {"data/data1.txt","data/data2.txt","data/data3.txt","data/data4.txt"};
	// vector<string> files = {"data/data1.txt","data/data2.txt","data/data3.txt","data/data4.txt","data/data5.txt","data/data6.txt","data/data7.txt","data/data8.txt","data/data9.txt","data/data10.txt","data/data11.txt","data/data12.txt"};

  for (int i = 0; i < files.size(); i++) {
		out_file_name = "result/Result_for_" + to_string(i + 1) + ".txt";

    out_file.open(out_file_name, ios::out | ios::app);

    if (out_file.good()) {
      cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% No. " << i + 1
           << " file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      cout << endl;

      out_file << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% No. " << i + 1
               << " file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

      file_name = files[i];

      t_begin = clock();
			
			// Read data from file
      Data data(file_name);

      domain = data.Get_Domain();
      size = data.Get_Size();
      running_time = data.Get_Time();
      beta_code = data.Get_Beta();
      accuracy = data.Get_Accuracy();
      surface = data.Get_Surface();
      method = data.Get_Method();
      mib_method = data.Get_MIB_method();
      equation = data.Get_Equation();

      out_file
          << "------------------------- Basic Info. ------------------------"
          << endl;
      out_file << setiosflags(ios::left) << setw(18) << "Mesh Size"
               << ": " << setw(5) << "NX = " << setiosflags(ios::left)
               << setw(5) << size[0] << setw(5)
               << " NY = " << setiosflags(ios::left) << setw(5) << size[1]
               << setw(5) << " NZ = " << setiosflags(ios::left) << setw(5)
               << size[2] << endl;
      out_file << setprecision(1) << scientific;
      out_file << setiosflags(ios::left) << setw(18) << "Time Step"
               << ": " << running_time[2] << endl;
      out_file << fixed;
      out_file << setiosflags(ios::left) << setw(18) << "Jump"
               << ": " << JP << endl;
      out_file << setiosflags(ios::left) << setw(18) << "Surface"
               << ": " << surface << endl;
      out_file << setiosflags(ios::left) << setw(18) << "Method"
               << ": " << method << endl;
      out_file << setiosflags(ios::left) << setw(18) << "MIB Method"
               << ": L" << mib_method << endl;
      out_file << setiosflags(ios::left) << setw(18) << "Beta No."
               << ": " << beta_code << endl;
      out_file << setiosflags(ios::left) << setw(18) << "Equation No."
               << ": " << equation << endl;
      out_file << setiosflags(ios::left) << setw(18) << "Accuracy"
               << ": " << accuracy << endl
               << endl;

      // Diffusion coefficients initialization
      if (beta_code == 0) {
        Beta_0 beta0;
        beta_ptr = &beta0;
      } else if (beta_code == 1) {
        Beta_1 beta1;
        beta_ptr = &beta1;
      } else if (beta_code == 2) {
        Beta_2 beta2;
        beta_ptr = &beta2;
      } else if (beta_code == 3) {
        Beta_3 beta3;
        beta_ptr = &beta3;
      } else {
        cout << "Beta is not found!" << endl;
        exit(0);
      }
      Beta &beta = *beta_ptr;

      // Surface initialization
      if (surface == 'C') {
        Surface_Cube ex_cube(arg_cube);
        ex_ptr = &ex_cube;
      } else if (surface == 'L') {
        Surface_Cylinder ex_cylinder(arg_cylinder);
        ex_ptr = &ex_cylinder;
      } else if (surface == 'E') {
        Surface_Ellipsoid ex_ellipsoid(arg_ellipsoid);
        ex_ptr = &ex_ellipsoid;
      } else if (surface == 'O') {
        Surface_Cone ex_cone(arg_cone);
        ex_ptr = &ex_cone;
      } else if (surface == 'P') {
        Surface_Pile ex_pile;
        ex_ptr = &ex_pile;
      } else if (surface == 'R') {
        Surface_Torus ex_torus(arg_torus);
        ex_ptr = &ex_torus;
      } else if (surface == 'D') {
        Surface_Dupin_Cyclide ex_dupin_cyclide(arg_dupin_cyclide);
        ex_ptr = &ex_dupin_cyclide;
      } else if (surface == 'M') {
        Surface_Molecular ex_molecular;
        ex_ptr = &ex_molecular;
      } else if (surface == 'H') {
        Surface_Heart ex_heart(arg_heart);
        ex_ptr = &ex_heart;
      } else if (surface == 'T') {
        Surface_Tanglecube ex_tanglecube(arg_tanglecube);
        ex_ptr = &ex_tanglecube;
      } else {
        cout << "Surface is not found!" << endl;
        exit(0);
      }

			// Surface construction
			Surface_Cartesian &ex = *ex_ptr;
			
			// Mesh construction
      Mesh mesh(domain, size, ex);
			
			// Spatial computation
      Intersections inter(ex, mesh, beta, accuracy, mib_method, out_file);

      cout << "Program is running......" << endl << endl;
			
			// Temporal computation, ADI / LOD / Trapezoidal Splitting
      if (method == 'A') {
        ADI_Starting(equation, accuracy, beta, mesh, inter, running_time, uh,
                     out_file);
			} else if (method == 'I' || method == 'C') {
        LOD_Starting(equation, beta, mesh, inter, running_time, uh, method, out_file);
      } else if (method == 'T') {
        TS_Starting(equation, beta, mesh, inter, running_time, uh, out_file);
      } else {
        cout << "Currently, only LOD and Douglas-ADI are applied" << endl;
        exit(0);
      }

      t_end = clock();
      t = ((double)(t_end - t_begin)) / CLOCKS_PER_SEC;

      out_file << setprecision(1) << endl;
      out_file << "CPU time cost: " << t << " seconds" << endl << endl;

			//Write_txt(inter,mesh,beta,uh,equation,running_time[1],i);
    } else {
      cout << "No such file is found!";
      exit(0);
    }

    out_file.close();
  }

  cout << "PROGRAM COMPLETED!" << endl;

  return 1;
}

/*********************************************************************************************
 Initialization of starting solution and do iterations

 INPUT
 equation : code of equation need to be used here
 accuracy : order of accuracy
 beta     : variable coefficients object
 mesh     : mesh object
 inter    : object of all intersection nodes
 time     : vector of 3 double values represent beginning time, finishing time and time step
 uh       : three dimensional uninitialized solution
 *********************************************************************************************/
void ADI_Starting(Int_I equation, Int_I accuracy, Beta &beta, Mesh &mesh,
                  Intersections &inter, VecDoub_I &time, CubicDoub_I &uh,
                  ofstream &out_file) {
  double t_start;
  Equation *eq_ptr;

  t_start = time[0];

  if (equation == 0) {
    Eq_0 eq0(t_start, beta);
    eq_ptr = &eq0;
  } else if (equation == 1) {
    Eq_1 eq1(t_start, beta);
    eq_ptr = &eq1;
  } else if (equation == 2) {
    Eq_2 eq2(t_start, beta);
    eq_ptr = &eq2;
  } else if (equation == 3) {
    Eq_3 eq3(t_start, beta);
    eq_ptr = &eq3;
  } else if (equation == 4) {
    Eq_4 eq4(t_start, beta);
    eq_ptr = &eq4;
  } else if (equation == 5) {
    Eq_5 eq5(t_start, beta);
    eq_ptr = &eq5;
  } else if (equation == 6) {
    Eq_6 eq6(t_start, beta);
    eq_ptr = &eq6;
  } else if (equation == 7) {
    Eq_7 eq7(t_start, beta);
    eq_ptr = &eq7;
  } else {
    cout << "No of equation is not found!" << endl;
    exit(0);
  }

  Equation &eq = *eq_ptr;

  Douglas_ADI adi(inter, mesh, beta, time, accuracy);

  adi.Initialization(eq, uh);

  ADI_Solver(equation, accuracy, beta, mesh, inter, adi, time, uh, out_file);
}

/*********************************************************************************************
 Solve problem in ADI scheme by iterations

 INPUT
 equation : code of equation need to be used here
 accuracy : order of accuracy
 beta     : variable coefficients object
 mesh     : mesh object
 inter    : object of all intersection nodes
 adi      : ADI object
 time     : vector of 3 double values represent beginning time, finishing time and time step
 uh       : three dimensional uninitialized solution
 *********************************************************************************************/
void ADI_Solver(Int_I equation, Int_I accruacy, Beta &beta, Mesh &mesh,
                Intersections &inter, Douglas_ADI &adi, VecDoub_I &time,
                CubicDoub_I &uh, ofstream &out_file) {
  double tnow, dt;
  int loop;
  Equation *eq_dt_ptr;

  out_file << "-------------------- Error of Douglas-ADI with MIB "
              "--------------------"
           << endl;
  tnow = time[0];
  dt = time[2];

  loop = (int)(time[1] / dt) / NPRINT;

  for (int i = 0; i < NPRINT; i++) {
    for (int j = 0; j < loop; j++) {
      tnow += dt;

      if (equation == 0) {
        Eq_0 eq0_dt(tnow, beta);
        eq_dt_ptr = &eq0_dt;
      } else if (equation == 1) {
        Eq_1 eq1_dt(tnow, beta);
        eq_dt_ptr = &eq1_dt;
      } else if (equation == 2) {
        Eq_2 eq2_dt(tnow, beta);
        eq_dt_ptr = &eq2_dt;
      } else if (equation == 3) {
        Eq_3 eq3_dt(tnow, beta);
        eq_dt_ptr = &eq3_dt;
      } else if (equation == 4) {
        Eq_4 eq4_dt(tnow, beta);
        eq_dt_ptr = &eq4_dt;
      } else if (equation == 5) {
        Eq_5 eq5_dt(tnow, beta);
        eq_dt_ptr = &eq5_dt;
      } else if (equation == 6) {
        Eq_6 eq6(tnow, beta);
        eq_dt_ptr = &eq6;
      } else if (equation == 7) {
        Eq_7 eq7(tnow, beta);
        eq_dt_ptr = &eq7;
      } else {
        cout << "No of equation is not found!" << endl;
        exit(0);
      }

      Equation &eq_dt = *eq_dt_ptr;

      adi.Solve_2nd(eq_dt, inter, uh, beta);
    }

    out_file << setprecision(1) << scientific << "T = " << tnow << endl;
    out_file << fixed;

    Equation &eq_now = *eq_dt_ptr;

    inter.Refresh_Fp(eq_now);
    inter.Error_Fp(out_file);

    inter.Refresh_Jump(eq_now, uh, beta);
    inter.Error_Jump(out_file);

    adi.Error(eq_now, uh, out_file);
  }
}

/*********************************************************************************************
 Initialization of starting solution and do iterations

 INPUT
 equation : code of equation need to be used here
 accuracy : order of accuracy
 beta     : variable coefficients object
 mesh     : mesh object
 inter    : object of all intersection nodes
 time     : vector of 3 double values represent beginning time, finishing time and time step
 uh       : three dimensional uninitialized solution
 *********************************************************************************************/
void LOD_Starting(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                  VecDoub_I &time, CubicDoub_I &uh, char &method, ofstream &out_file) {
  double t_start;
  Equation *eq_ptr;

  t_start = time[0];

  if (equation == 0) {
    Eq_0 eq0(t_start, beta);
    eq_ptr = &eq0;
  } else if (equation == 1) {
    Eq_1 eq1(t_start, beta);
    eq_ptr = &eq1;
  } else if (equation == 2) {
    Eq_2 eq2(t_start, beta);
    eq_ptr = &eq2;
  } else if (equation == 3) {
    Eq_3 eq3(t_start, beta);
    eq_ptr = &eq3;
  } else if (equation == 4) {
    Eq_4 eq4(t_start, beta);
    eq_ptr = &eq4;
  } else if (equation == 5) {
    Eq_5 eq5(t_start, beta);
    eq_ptr = &eq5;
  } else if (equation == 6) {
    Eq_6 eq6(t_start, beta);
    eq_ptr = &eq6;
  } else if (equation == 7) {
    Eq_7 eq7(t_start, beta);
    eq_ptr = &eq7;
  } else {
    cout << "No of equation is not found!" << endl;
    exit(0);
  }

  Equation &eq = *eq_ptr;

  LOD lod(inter, mesh, beta, time);

  lod.Initialization(eq, uh);

  LOD_Solver(equation, beta, mesh, inter, lod, time, uh, method, out_file);
}

/*********************************************************************************************
 Solve problem in LOD scheme by iterations

 INPUT
 equation : code of equation need to be used here
 accuracy : order of accuracy
 beta     : variable coefficients object
 mesh     : mesh object
 inter    : object of all intersection nodes
 lod      : LOD object
 time     : vector of 3 double values represent beginning time, finishing time and time step
 uh       : three dimensional uninitialized solution
 *********************************************************************************************/
void LOD_Solver(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                LOD &lod, VecDoub_I &time, CubicDoub_I &uh, char &method,
                ofstream &out_file) {
  double tnow, dt;
  int loop;
  Equation *eq_dt_ptr;

  out_file << "-------------------- Error of LOD with MIB--------------------"
           << endl;
  tnow = time[0];
  dt = time[2];

  loop = (int)(time[1] / dt) / NPRINT;

  for (int i = 0; i < NPRINT; i++) {
    for (int j = 0; j < loop; j++) {
      tnow += dt;

      if (equation == 0) {
        Eq_0 eq0_dt(tnow, beta);
        eq_dt_ptr = &eq0_dt;
      } else if (equation == 1) {
        Eq_1 eq1_dt(tnow, beta);
        eq_dt_ptr = &eq1_dt;
      } else if (equation == 2) {
        Eq_2 eq2_dt(tnow, beta);
        eq_dt_ptr = &eq2_dt;
      } else if (equation == 3) {
        Eq_3 eq3_dt(tnow, beta);
        eq_dt_ptr = &eq3_dt;
      } else if (equation == 4) {
        Eq_4 eq4_dt(tnow, beta);
        eq_dt_ptr = &eq4_dt;
      } else if (equation == 5) {
        Eq_5 eq5_dt(tnow, beta);
        eq_dt_ptr = &eq5_dt;
      } else if (equation == 6) {
        Eq_6 eq6(tnow, beta);
        eq_dt_ptr = &eq6;
      } else if (equation == 7) {
        Eq_7 eq7(tnow, beta);
        eq_dt_ptr = &eq7;
      } else {
        cout << "No of equation is not found!" << endl;
        exit(0);
      }

      Equation &eq_dt = *eq_dt_ptr;

      lod.Solve_2nd(eq_dt, inter, uh, beta, method);
    }

    out_file << setprecision(1) << scientific << "T = " << tnow << endl;
    out_file << fixed;

    Equation &eq_now = *eq_dt_ptr;

    inter.Refresh_Fp(eq_now);
    inter.Error_Fp(out_file);

    inter.Refresh_Jump(eq_now, uh, beta);
    inter.Error_Jump(out_file);

    lod.Error(eq_now, uh, out_file);
  }
}

/*********************************************************************************************
 Initialization of starting solution and do iterations

 INPUT
 equation : code of equation need to be used here
 accuracy : order of accuracy
 beta     : variable coefficients object
 mesh     : mesh object
 inter    : object of all intersection nodes
 time     : vector of 3 double values represent beginning time, finishing time and time step
 uh       : three dimensional uninitialized solution
 *********************************************************************************************/
void TS_Starting(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                 VecDoub_I &time, CubicDoub_I &uh, ofstream &out_file) {
  double t_start;
  Equation *eq_ptr;

  t_start = time[0];

  if (equation == 0) {
    Eq_0 eq0(t_start, beta);
    eq_ptr = &eq0;
  } else if (equation == 1) {
    Eq_1 eq1(t_start, beta);
    eq_ptr = &eq1;
  } else if (equation == 2) {
    Eq_2 eq2(t_start, beta);
    eq_ptr = &eq2;
  } else if (equation == 3) {
    Eq_3 eq3(t_start, beta);
    eq_ptr = &eq3;
  } else if (equation == 4) {
    Eq_4 eq4(t_start, beta);
    eq_ptr = &eq4;
  } else if (equation == 5) {
    Eq_5 eq5(t_start, beta);
    eq_ptr = &eq5;
  } else if (equation == 6) {
    Eq_6 eq6(t_start, beta);
    eq_ptr = &eq6;
  } else if (equation == 7) {
    Eq_7 eq7(t_start, beta);
    eq_ptr = &eq7;
  } else {
    cout << "No of equation is not found!" << endl;
    exit(0);
  }

  Equation &eq = *eq_ptr;

  TS ts(inter, mesh, beta, time);

  ts.Initialization(eq, uh);

  TS_Solver(equation, beta, mesh, inter, ts, time, uh, out_file);
}

/*********************************************************************************************
 Solve problem in Trapezoidal Splitting scheme by iterations

 INPUT
 equation : code of equation need to be used here
 accuracy : order of accuracy
 beta     : variable coefficients object
 mesh     : mesh object
 inter    : object of all intersection nodes
 ts       : Trapezoidal splitting object
 time     : vector of 3 double values represent beginning time, finishing time and time step
 uh       : three dimensional uninitialized solution
 *********************************************************************************************/
void TS_Solver(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
               TS &ts, VecDoub_I &time, CubicDoub_I &uh, ofstream &out_file) {
  double tnow, dt;
  int loop;
  Equation *eq_ptr, *eq_half_dt_ptr, *eq_dt_ptr;

  out_file << "-------------------- Error of Trapezoidal Splitting with MIB "
              "--------------------"
           << endl;
  tnow = time[0];
  dt = time[2];

  loop = (int)(time[1] / dt) / NPRINT;

  for (int i = 0; i < NPRINT; i++) {
    for (int j = 0; j < loop; j++) {
      tnow += dt;

      if (equation == 0) {
        Eq_0 eq0(tnow - dt, beta);
        Eq_0 eq0_half_dt(tnow - 0.5 * dt, beta);
        Eq_0 eq0_dt(tnow, beta);
        eq_ptr = &eq0;
        eq_half_dt_ptr = &eq0_half_dt;
        eq_dt_ptr = &eq0_dt;
      } else if (equation == 1) {
        Eq_1 eq1(tnow - dt, beta);
        Eq_1 eq1_half_dt(tnow - 0.5 * dt, beta);
        Eq_1 eq1_dt(tnow, beta);
        eq_ptr = &eq1;
        eq_half_dt_ptr = &eq1_half_dt;
        eq_dt_ptr = &eq1_dt;
      } else if (equation == 2) {
        Eq_2 eq2(tnow - dt, beta);
        Eq_2 eq2_half_dt(tnow - 0.5 * dt, beta);
        Eq_2 eq2_dt(tnow, beta);
        eq_ptr = &eq2;
        eq_half_dt_ptr = &eq2_half_dt;
        eq_dt_ptr = &eq2_dt;
      } else if (equation == 3) {
        Eq_3 eq3(tnow - dt, beta);
        Eq_3 eq3_half_dt(tnow - 0.5 * dt, beta);
        Eq_3 eq3_dt(tnow, beta);
        eq_ptr = &eq3;
        eq_half_dt_ptr = &eq3_half_dt;
        eq_dt_ptr = &eq3_dt;
      } else if (equation == 4) {
        Eq_4 eq4(tnow - dt, beta);
        Eq_4 eq4_half_dt(tnow - 0.5 * dt, beta);
        Eq_4 eq4_dt(tnow, beta);
        eq_ptr = &eq4;
        eq_half_dt_ptr = &eq4_half_dt;
        eq_dt_ptr = &eq4_dt;
      } else if (equation == 5) {
        Eq_5 eq5(tnow - dt, beta);
        Eq_5 eq5_half_dt(tnow - 0.5 * dt, beta);
        Eq_5 eq5_dt(tnow, beta);
        eq_ptr = &eq5;
        eq_half_dt_ptr = &eq5_half_dt;
        eq_dt_ptr = &eq5_dt;
      } else if (equation == 6) {
        Eq_6 eq6(tnow - dt, beta);
        Eq_6 eq6_half_dt(tnow - 0.5 * dt, beta);
        Eq_6 eq6_dt(tnow, beta);
        eq_ptr = &eq6;
        eq_half_dt_ptr = &eq6_half_dt;
        eq_dt_ptr = &eq6_dt;
      } else if (equation == 7) {
        Eq_7 eq7(tnow - dt, beta);
        Eq_7 eq7_half_dt(tnow - 0.5 * dt, beta);
        Eq_7 eq7_dt(tnow, beta);
        eq_ptr = &eq7;
        eq_half_dt_ptr = &eq7_half_dt;
        eq_dt_ptr = &eq7_dt;
      } else {
        cout << "No of equation is not found!" << endl;
        exit(0);
      }

      Equation &eq = *eq_ptr;
      Equation &eq_half_dt = *eq_half_dt_ptr;
      Equation &eq_dt = *eq_dt_ptr;

      ts.Solve_2nd(eq, eq_half_dt, eq_dt, inter, uh, beta);
    }

    out_file << setprecision(1) << scientific << "T = " << tnow << endl;
    out_file << fixed;

    Equation &eq_now = *eq_dt_ptr;

    inter.Refresh_Fp(eq_now);
    inter.Error_Fp(out_file);

    inter.Refresh_Jump(eq_now, uh, beta);
    inter.Error_Jump(out_file);

    ts.Error(eq_now, uh, out_file);
  }
}

/*********************************************************************************************
 Write the function value for closest outside node around interface

 INPUT
 inter    : object of all intersection nodes
 mesh     : mesh object
 beta     : variable coefficients object
 uh       : three dimensional uninitialized solution
 equation : code of equation need to be used here
 t        : current time
 i        : file number
 *********************************************************************************************/
void Write_txt(Intersections &inter, Mesh &mesh, Beta &beta, CubicDoub_I &uh,
               Int_I equation, Doub_I t, Int_I i) {
  ofstream err_txt_file, sol_txt_file;
  string err_txt_file_name, sol_txt_file_name;
  Equation *eq_ptr;

  int ix, iy, iz;
  double error;

  sol_txt_file_name = "result/Numerical Solutions<" + to_string(i) + ">.txt";
  err_txt_file_name = "result/Numerical Error<" + to_string(i) + ">.txt";

  sol_txt_file.open(sol_txt_file_name, ios::out);
  err_txt_file.open(err_txt_file_name, ios::out);

  if (equation == 0) {
    Eq_0 eq0(t, beta);
    eq_ptr = &eq0;
  } else if (equation == 1) {
    Eq_1 eq1(t, beta);
    eq_ptr = &eq1;
  } else if (equation == 2) {
    Eq_2 eq2(t, beta);
    eq_ptr = &eq2;
  } else if (equation == 3) {
    Eq_3 eq3(t, beta);
    eq_ptr = &eq3;
  } else if (equation == 4) {
    Eq_4 eq4(t, beta);
    eq_ptr = &eq4;
  } else if (equation == 5) {
    Eq_5 eq5(t, beta);
    eq_ptr = &eq5;
  } else if (equation == 6) {
    Eq_6 eq6(t, beta);
    eq_ptr = &eq6;
  } else if (equation == 7) {
    Eq_7 eq7(t, beta);
    eq_ptr = &eq7;
  } else {
    cout << "No of equation is not found!" << endl;
    exit(0);
  }

  Equation &eq = *eq_ptr;

  if (sol_txt_file.good() && err_txt_file.good()) {
    sol_txt_file << setiosflags(ios::right) << setprecision(6) << fixed;
    err_txt_file << setiosflags(ios::right) << setprecision(6) << fixed;

    // Z-direction
    for (int ix = 1; ix < mesh.nx - 1; ix++) {
      for (int iy = 1; iy < mesh.ny - 1; iy++) {
        for (int ip = 0; ip < inter.ifpz[ix][iy].size(); ip++) {
          // Always use clostest outside grid node value
          if (abs(inter.ifpz[ix][iy][ip].ID) % 2 == 1) {
            iz = inter.ifpz[ix][iy][ip].left_loc;
          } else if (abs(inter.ifpz[ix][iy][ip].ID) % 2 == 0) {
            iz = inter.ifpz[ix][iy][ip].left_loc + 1;

          } else {
            cout << "X Wrong ID is found" << endl;
            exit(0);
          }

          error = abs(uh[ix][iy][iz] -
                      eq.Outer_u(mesh.xi[ix], mesh.yi[iy], mesh.zi[iz]));

          sol_txt_file << setw(5) << ix << setw(5) << iy << setw(5) << iz
                       << setw(20) << uh[ix][iy][iz] << endl;
          err_txt_file << setw(5) << ix << setw(5) << iy << setw(5) << iz
                       << setw(20) << error << endl;
        }
      }
    }

    // Y-direction
    for (int ix = 1; ix < mesh.nx - 1; ix++) {
      for (int iz = 1; iz < mesh.nz - 1; iz++) {
        for (int ip = 0; ip < inter.ifpy[ix][iz].size(); ip++) {
          // Always use clostest outside grid node value
          if (abs(inter.ifpy[ix][iz][ip].ID) % 2 == 1) {
            iy = inter.ifpy[ix][iz][ip].left_loc;
          } else if (abs(inter.ifpy[ix][iz][ip].ID) % 2 == 0) {
            iy = inter.ifpy[ix][iz][ip].left_loc + 1;
          } else {
            cout << "Y Wrong ID is found" << endl;
            exit(0);
          }

          error = abs(uh[ix][iy][iz] -
                      eq.Outer_u(mesh.xi[ix], mesh.yi[iy], mesh.zi[iz]));

          sol_txt_file << setw(5) << ix << setw(5) << iy << setw(5) << iz
                       << setw(20) << uh[ix][iy][iz] << endl;
          err_txt_file << setw(5) << ix << setw(5) << iy << setw(5) << iz
                       << setw(20) << error << endl;
        }
      }
    }

    // X-direction
    for (int iy = 1; iy < mesh.ny - 1; iy++) {
      for (int iz = 1; iz < mesh.nz - 1; iz++) {
        for (int ip = 0; ip < inter.ifpx[iy][iz].size(); ip++) {
          // Always use clostest outside grid node value
          if (abs(inter.ifpx[iy][iz][ip].ID) % 2 == 1) {
            ix = inter.ifpx[iy][iz][ip].left_loc;
          } else if (abs(inter.ifpx[iy][iz][ip].ID) % 2 == 0) {
            ix = inter.ifpx[iy][iz][ip].left_loc + 1;
          } else {
            cout << "Wrong ID is found" << endl;
            exit(0);
          }

          error = abs(uh[ix][iy][iz] -
                      eq.Outer_u(mesh.xi[ix], mesh.yi[iy], mesh.zi[iz]));

          sol_txt_file << setw(5) << ix << setw(5) << iy << setw(5) << iz
                       << setw(20) << uh[ix][iy][iz] << endl;
          err_txt_file << setw(5) << ix << setw(5) << iy << setw(5) << iz
                       << setw(20) << error << endl;
        }
      }
    }

  } else {
    cout << "No such file is found!";
    exit(0);
  }

  sol_txt_file.close();
  err_txt_file.close();
}