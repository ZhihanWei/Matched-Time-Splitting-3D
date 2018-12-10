#pragma once

#include <string>
#include "constant.h"
#include "data/data_type.h"

class Data {
private:
  enum Config {
    MAX_X,
    MIN_X,
    MAX_Y,
    MIN_Y,
    MAX_Z,
    MIN_Z,
    TIME_START,
    TIME_TERMINATE,
    TIME_STEP,
    NX,
    NY,
    NZ,
    SURFACE,
    EQUATION,
    TEMPORAL_METHOD,
    SPATIAL_METHOD,
    SPATIAL_ACCURACY,
    DIFFUSION_COEFFICIENT,
    COMMENT,
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

  Config translate(const string &in_string);

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
