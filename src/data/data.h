#pragma once

#include <string>
#include "constant.h"
#include "data/data_type.h"

class Data {
private:
  enum Prt_name {
    e_xl,
    e_xr,
    e_yl,
    e_yr,
    e_zl,
    e_zr,
    e_t_start,
    e_t_finish,
    e_t_step,
    e_nx,
    e_ny,
    e_nz,
    e_beta_inside,
    e_beta_outside,
    e_surface,
    e_equation,
    e_method,
    e_mib_method,
    e_accuracy,
    e_beta
  };

  double xl, xr, yl, yr, zl, zr;
  double t_start, t_finish, t_step;
  int nx, ny, nz;
  char surface;
  char method;
  int mib_method;
  int equation;
  int accuracy;
  int beta;

  Prt_name Translation(const string &in_string);

public:
  Data(const string &file_name);

  void Display();

  VecDoub Get_Domain() const;
  VecInt Get_Size() const;
  VecDoub Get_Time() const;
  char Get_Method() const;
  char Get_Surface() const;
  int Get_Beta() const;
  int Get_Accuracy() const;
  int Get_MIB_method() const;
  int Get_Equation() const;
  double Get_Tol() const;
};
