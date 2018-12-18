#pragma once

#include <string>
#include "constant.h"

class ReadConfig {
private:
  double x_max, x_min, y_max, y_min, z_max, z_min;
  double t_start, t_terminate, t_step;
  int nx, ny, nz;
  int spatial_accuracy, equation, diffusion_coef;
  Surface_Type surface;
  Temporal_Method_Type temporal_method;
  Spatial_Method_Type spatial_method;

  Config Translate(const string &);
  Surface_Type ParseSurface(const string &);
  Temporal_Method_Type ParseTemporalMethod(const string &);
  Spatial_Method_Type ParseSpatialMethod(const string &);
  int ParseEquationMethod(const string &);
  int ParseDiffusionCoefficient(const string &);
  int ParseSpatialAccuracy(const string &);

public:
  ReadConfig(const string &file_name);

  void Display();
  VecDoub GetDomain() const;
  VecInt GetMesh() const;
  VecDoub GetTime() const;
  int GetEquation() const;
  int GetDiffusionCoef() const;
  int GetSpatialAccuracy() const;
  Spatial_Method_Type GetSpatialMethod() const;
  Temporal_Method_Type GetTemporalMethod() const;
  Spatial_Method_Type GetSpatialmethod() const;
  Surface_Type GetSurface() const;
};
