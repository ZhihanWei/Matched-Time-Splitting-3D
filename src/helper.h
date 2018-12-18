#pragma once

#include <iomanip>
#include <memory>
#include "configuration/read_config.h"
#include "constant.h"
#include "mesh/mesh.h"
#include "spatial/intersection.h"
#include "temporal/douglas_adi.h"
#include "temporal/lod.h"
#include "temporal/trapezoidal.h"
// User-defined diffusion coefficient beta
#include "diffusion/beta.h"
#include "diffusion/beta_0.h"
#include "diffusion/beta_1.h"
#include "diffusion/beta_2.h"
#include "diffusion/beta_3.h"
#include "diffusion/beta_4.h"
// Users-defined functions
#include "equation/eq_0.h"
#include "equation/eq_1.h"
#include "equation/eq_2.h"
#include "equation/eq_3.h"
#include "equation/eq_4.h"
#include "equation/eq_5.h"
#include "equation/eq_6.h"
#include "equation/eq_7.h"
#include "equation/equation.h"
// User-defined surfaces
#include "surface/surface_cartesian.h"
#include "surface/surface_cone.h"
#include "surface/surface_cube.h"
#include "surface/surface_cylinder.h"
#include "surface/surface_dupin_cyclide.h"
#include "surface/surface_ellipsoid.h"
#include "surface/surface_heart.h"
#include "surface/surface_molecular.h"
#include "surface/surface_pile.h"
#include "surface/surface_tanglecube.h"
#include "surface/surface_torus.h"

#define LOG_FATAL(msg)                                                                         \
  std::cerr << "@ " << __FILE__ << " line: " << __LINE__ << "; [ERROR]: " << msg << std::endl; \
  exit(EXIT_FAILURE);

using namespace std;

inline void print_help(char *binary_name) {
  cout << "Usage: " << binary_name << " -i input [-o output] \n \n";
  cout << "-i, --input [required] \n"
          "      input configuration file, sample can be found in directory \"example\" \n"
          "-o, --output [required] \n"
          "      output file \n";
}

inline unique_ptr<Equation> find_equation(double t, int equation_type, Beta &beta) {
  unique_ptr<Equation> ptr;

  switch (equation_type) {
    case 0:
      ptr = make_unique<Eq_0>(t, beta);
      break;
    case 1:
      ptr = make_unique<Eq_1>(t, beta);
      break;
    case 2:
      ptr = make_unique<Eq_2>(t, beta);
      break;
    case 3:
      ptr = make_unique<Eq_3>(t, beta);
      break;
    case 4:
      ptr = make_unique<Eq_4>(t, beta);
      break;
    case 5:
      ptr = make_unique<Eq_5>(t, beta);
      break;
    case 6:
      ptr = make_unique<Eq_6>(t, beta);
      break;
    case 7:
      ptr = make_unique<Eq_7>(t, beta);
      break;
    default:
      LOG_FATAL("Provided equation type doesn't exist, check configuration");
      break;
  }

  return ptr;
}

inline unique_ptr<Beta> find_diffusion_coefficient(int diffusion_coeff_type) {
  unique_ptr<Beta> ptr;

  switch (diffusion_coeff_type) {
    case 0:
      ptr = make_unique<Beta_0>();
      break;
    case 1:
      ptr = make_unique<Beta_1>();
      break;
    case 2:
      ptr = make_unique<Beta_2>();
      break;
    case 3:
      ptr = make_unique<Beta_3>();
      break;
    case 4:
      ptr = make_unique<Beta_4>();
      break;
    default:
      LOG_FATAL("Provided diffusion coefficient type doesn't exist, check configuration!");
      break;
  }
  return ptr;
}

inline unique_ptr<Surface_Cartesian> find_surface(Surface_Type surface) {
  unique_ptr<Surface_Cartesian> ptr;

  switch (surface) {
    case Surface_Type::TANGLECUBE: {
      // suggested domain: [-4,4;-4,4;-5,5]
      VecDoub args{1.0, -5.0, 1.0, -5.0, 1.0, -5.0, 10.0};
      ptr = make_unique<Surface_Tanglecube>(args);
      break;
    }
    case Surface_Type::CUBE: {
      // suggested domain [-2,2;-2,2;-2,2]
      VecDoub args{-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
      ptr = make_unique<Surface_Cube>(args);
      break;
    }
    case Surface_Type::CYLINDER: {
      // suggested domain: [-4,4;-4,4;-4,4]
      VecDoub args{2.0, -3.0, 3.0};
      ptr = make_unique<Surface_Cylinder>(args);
      break;
    }
    case Surface_Type::CONE: {
      // suggested domain: [-4,4;-4,4;-4,1]
      VecDoub args{1.0, 0.0, 0.0, -3.0, 0.0};
      ptr = make_unique<Surface_Cone>(args);
      break;
    }
    case Surface_Type::ELLIPSOID: {
      // suggested domain: [-4,4;-2,2;-2,2]
      VecDoub args{2.0, 1.0, 1.0, 1.0};
      ptr = make_unique<Surface_Ellipsoid>(args);
      break;
    }
    case Surface_Type::TORUS: {
      // suggested domain: [-4,4;-4,4;-4,4]
      VecDoub args{1.5, 0.7};
      ptr = make_unique<Surface_Torus>(args);
      break;
    }
    case Surface_Type::DUPIN_CYCLIDE: {
      // suggested domain: [-5,5;-5,5;-5,5]
      VecDoub args{1.0, 0.3, 0.6};
      ptr = make_unique<Surface_Dupin_Cyclide>(args);
      break;
    }
    case Surface_Type::HEART: {
      // domain [-4,4;-2,2;-4,4], 'H'
      VecDoub args{1.0, 9.0 / 4.0, 1.0, -1.0, -1.0, -9.0 / 80.0};
      ptr = make_unique<Surface_Heart>(args);
      break;
    }
    case Surface_Type::MOLECULAR: {
      ptr = make_unique<Surface_Molecular>();
      break;
    }
    case Surface_Type::PILE: {
      ptr = make_unique<Surface_Pile>();
      break;
    }
    default:
      LOG_FATAL("Provided surface type doesn't exist, check configuration!");
      break;
  }

  return ptr;
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
 time     : vector of 3 double values represent beginning time, finishing time
 and time step uh       : three dimensional uninitialized solution
 *********************************************************************************************/
inline void ADI_Solver(Int_I equation, Int_I accruacy, Beta &beta, Mesh &mesh,
                       Intersections &inter, Douglas_ADI &adi, VecDoub_I &time,
                       CubicDoub_I &uh, ofstream &out_stream) {
  double dt = time[2], t_now = time[0], t_finish = time[1];
  int loop;

  out_stream << "----------- Error of Douglas-ADI with MIB -----------" << endl;

  loop = (int)(t_finish / dt) / NPRINT;

  for (int i = 0; i < NPRINT; i++) {
    for (int j = 0; j < loop; j++) {
      t_now += dt;

      unique_ptr<Equation> eq_ptr = move(find_equation(t_now, equation, beta));
      if (eq_ptr.get() == nullptr) {
        LOG_FATAL("Equation pointer is NULL");
      }

      adi.Solve_2nd(*eq_ptr, inter, uh, beta);
    }

    out_stream << setprecision(1) << scientific << "T = " << t_now << endl;
    out_stream << fixed;

    unique_ptr<Equation> eq_ptr = move(find_equation(t_now, equation, beta));
    if (eq_ptr.get() == nullptr) {
      LOG_FATAL("Equation pointer is NULL");
    }

    inter.Refresh_Fp(*eq_ptr);
    inter.Error_Fp(out_stream);

    inter.Refresh_Jump(*eq_ptr, uh, beta);
    inter.Error_Jump(out_stream);

    adi.Error(*eq_ptr, uh, out_stream);
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
 time     : vector of 3 double values represent beginning time, finishing time
 and time step uh       : three dimensional uninitialized solution
 *********************************************************************************************/
inline void Call_ADI(Int_I equation, Int_I accuracy, Beta &beta, Mesh &mesh,
                     Intersections &inter, VecDoub_I &time, CubicDoub_I &uh,
                     ofstream &out_stream) {
  unique_ptr<Equation> eq_ptr = move(find_equation(time[0], equation, beta));
  if (eq_ptr.get() == nullptr) {
    LOG_FATAL("Equation pointer is NULL")
  }

  Douglas_ADI adi(inter, mesh, beta, time, accuracy);

  adi.Initialization(*eq_ptr, uh);

  ADI_Solver(equation, accuracy, beta, mesh, inter, adi, time, uh, out_stream);
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
 time     : vector of 3 double values represent beginning time, finishing time
 and time step uh       : three dimensional uninitialized solution
 *********************************************************************************************/
inline void LOD_Solver(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                       LOD &lod, VecDoub_I &time, CubicDoub_I &uh, char &method,
                       ofstream &out_stream) {
  double dt = time[2], t_now = time[0], t_finish = time[1];
  int loop;

  out_stream << "------------ Error of LOD with MIB -------------" << endl;

  loop = (int)(t_finish / dt) / NPRINT;

  for (int i = 0; i < NPRINT; i++) {
    for (int j = 0; j < loop; j++) {
      t_now += dt;

      unique_ptr<Equation> eq_ptr = move(find_equation(t_now, equation, beta));
      if (eq_ptr.get() == nullptr) {
        LOG_FATAL("Equation pointer is NULL");
      }

      lod.Solve_2nd(*eq_ptr, inter, uh, beta, method);
    }

    out_stream << setprecision(1) << scientific << "T = " << t_now << endl;
    out_stream << fixed;

    unique_ptr<Equation> eq_ptr = move(find_equation(t_now, equation, beta));
    if (eq_ptr.get() == nullptr) {
      LOG_FATAL("Equation pointer is NULL");
    }

    inter.Refresh_Fp(*eq_ptr);
    inter.Error_Fp(out_stream);

    inter.Refresh_Jump(*eq_ptr, uh, beta);
    inter.Error_Jump(out_stream);

    lod.Error(*eq_ptr, uh, out_stream);
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
 time     : vector of 3 double values represent beginning time, finishing time
 and time step uh       : three dimensional uninitialized solution
 *********************************************************************************************/
inline void Call_LOD(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                     VecDoub_I &time, CubicDoub_I &uh, char &method,
                     ofstream &out_stream) {
  unique_ptr<Equation> eq_ptr = move(find_equation(time[0], equation, beta));
  if (eq_ptr.get() == nullptr) {
    LOG_FATAL("Equation pointer is NULL");
  }

  LOD lod(inter, mesh, beta, time);

  lod.Initialization(*eq_ptr, uh);

  LOD_Solver(equation, beta, mesh, inter, lod, time, uh, method, out_stream);
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
 time     : vector of 3 double values represent beginning time, finishing time
 and time step uh       : three dimensional uninitialized solution
 *********************************************************************************************/
inline void TS_Solver(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                      TS &ts, VecDoub_I &time, CubicDoub_I &uh, ofstream &out_stream) {
  double dt = time[2], t_now = time[0], t_finish = time[1];
  int loop;

  out_stream << "--------------- Error of TS with MIB ---------------" << endl;

  loop = (int)(t_finish / dt) / NPRINT;

  for (int i = 0; i < NPRINT; i++) {
    for (int j = 0; j < loop; j++) {
      t_now += dt;

      unique_ptr<Equation> eq_ptr = move(find_equation(t_now - dt, equation, beta));
      unique_ptr<Equation> eq_ptr_half_dt = move(find_equation(t_now - 0.5 * dt, equation, beta));
      unique_ptr<Equation> eq_ptr_dt = move(find_equation(t_now, equation, beta));
      if (eq_ptr.get() == nullptr || eq_ptr_half_dt.get() == nullptr || eq_ptr_dt.get() == nullptr) {
        LOG_FATAL("Equation pointer is NULL");
      }

      ts.Solve_2nd(*eq_ptr, *eq_ptr_half_dt, *eq_ptr_dt, inter, uh, beta);
    }

    out_stream << setprecision(1) << scientific << "T = " << t_now << endl;
    out_stream << fixed;

    unique_ptr<Equation> eq_ptr_dt = move(find_equation(t_now, equation, beta));
    if (eq_ptr_dt.get() == nullptr) {
      LOG_FATAL("Equation pointer is NULL");
    }

    inter.Refresh_Fp(*eq_ptr_dt);
    inter.Error_Fp(out_stream);

    inter.Refresh_Jump(*eq_ptr_dt, uh, beta);
    inter.Error_Jump(out_stream);

    ts.Error(*eq_ptr_dt, uh, out_stream);
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
 time     : vector of 3 double values represent beginning time, finishing time
 and time step uh       : three dimensional uninitialized solution
 *********************************************************************************************/
inline void Call_TS(Int_I equation, Beta &beta, Mesh &mesh, Intersections &inter,
                    VecDoub_I &time, CubicDoub_I &uh, ofstream &out_stream) {
  unique_ptr<Equation> eq_ptr = move(find_equation(time[0], equation, beta));
  if (eq_ptr.get() == nullptr) {
    LOG_FATAL("Equation pointer is NULL");
  }

  TS ts(inter, mesh, beta, time);

  ts.Initialization(*eq_ptr, uh);

  TS_Solver(equation, beta, mesh, inter, ts, time, uh, out_stream);
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

/*
inline void Write_txt(Intersections &inter, Mesh &mesh, Beta &beta, CubicDoub_I &uh,
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
        for (size_t ip = 0; ip < inter.ifpz[ix][iy].size(); ip++) {
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
        for (size_t ip = 0; ip < inter.ifpy[ix][iz].size(); ip++) {
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
        for (size_t ip = 0; ip < inter.ifpx[iy][iz].size(); ip++) {
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
*/
