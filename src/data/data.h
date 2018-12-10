#pragma once

#include <string>
#include "constant.h"
#include "data/data_type.h"

class Data {
private:
  double x_max, x_min, y_max, y_min, z_max, z_min;
  double t_start, t_terminate, t_step;
  int nx, ny, nz;
  string surface, temporal_method, spatial_method;
  int spatial_accuracy, equation, diffusion_coef;

  Config Translate(const string &);
  string ParseSurface(const string &);
  string ParseTemporalMethod(const string &);
  string ParseSpatialMethod(const string &);
  int ParseEquationMethod(const string &);
  int ParseDiffusionCoefficient(const string &);
  int ParseSpatialAccuracy(const string &);

public:
  Data(const string &file_name);

  void Display();
  VecDoub GetDomain() const;
  VecInt GetMesh() const;
  VecDoub GetTime() const;
  char GetTemporalMethod() const;
  char GetSurface() const;
  int GetDiffusionCoeff() const;
  int GetSpatialAccuracy() const;
  int GetSpatialMethod() const;
  int GetEquation() const;
};
