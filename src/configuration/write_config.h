#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "configuration/read_config.h"
#include "constant.h"
#include "helper.h"

string write_surface(const Surface_Type s) {
  string r;

  switch (s) {
    case Surface_Type::TANGLECUBE:
      r = "tanglecube";
      break;
    case Surface_Type::CUBE:
      r = "cube";
      break;
    case Surface_Type::CYLINDER:
      r = "cylinder";
      break;
    case Surface_Type::ELLIPSOID:
      r = "ellipsoid";
      break;
    case Surface_Type::CONE:
      r = "cone";
      break;
    case Surface_Type::PILE:
      r = "pile";
      break;
    case Surface_Type::TORUS:
      r = "torus";
      break;
    case Surface_Type::DUPIN_CYCLIDE:
      r = "dupin_cyclide";
      break;
    case Surface_Type::MOLECULAR:
      r = "molecular";
      break;
    case Surface_Type::HEART:
      r = "heart";
      break;
    default:
      LOG_FATAL("Corresponding string of given surface doesn't exist, cannot be writen to stream");
      break;
  }

  return r;
}

string write_spatial_method(const Spatial_Method_Type s) {
  string r;

  switch (s) {
    case Spatial_Method_Type::MIB_V1:
      r = "mib_v1";
      break;
    case Spatial_Method_Type::MIB_V2:
      r = "mib_v2";
      break;
    default:
      LOG_FATAL("Corresponding string of given spatial method doesn't exist, cannot be writen to stream");
      break;
  }

  return r;
}

string write_temporal_method(const Temporal_Method_Type s) {
  string r;

  switch (s) {
    case Temporal_Method_Type::ADI:
      r = "ADI";
      break;
    case Temporal_Method_Type::LOD_IE:
      r = "LOD-IE";
      break;
    case Temporal_Method_Type::LOD_CN:
      r = "LOD-CN";
      break;
    case Temporal_Method_Type::TS:
      r = "TS";
      break;
    default:
      LOG_FATAL("Corresponding string of given temporal method doesn't exist, cannot be writen to stream");
      break;
  }

  return r;
}

void write_configuration(ofstream& out_stream, ReadConfig& config) {
  const VecDoub domain = config.GetDomain();
  const VecInt size = config.GetMesh();
  const VecDoub time_info = config.GetTime();
  const int diffusion_coef_type = config.GetDiffusionCoef();
  const int equation_type = config.GetEquation();
  const int spatial_accuracy_type = config.GetSpatialAccuracy();
  const Surface_Type surface = config.GetSurface();
  const Temporal_Method_Type temporal_method = config.GetTemporalMethod();
  const Spatial_Method_Type spatial_method = config.GetSpatialmethod();

  out_stream << "------------ Configuration informations ------------" << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Mesh Size"
             << ": " << setw(5) << "NX = " << setiosflags(ios::left)
             << setw(5) << size.at(0) << setw(5)
             << " NY = " << setiosflags(ios::left) << setw(5) << size.at(1)
             << setw(5) << " NZ = " << setiosflags(ios::left) << setw(5)
             << size.at(2) << endl;
  out_stream << setprecision(1) << scientific;
  out_stream << setiosflags(ios::left) << setw(18) << "Time Step"
             << ": " << time_info.at(2) << endl;
  out_stream << fixed;
  out_stream << setiosflags(ios::left) << setw(18) << "Jump"
             << ": " << JP << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Diffusion Coefficient Type "
             << ": " << diffusion_coef_type << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Equation Type"
             << ": " << equation_type << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Order of Spatial Accuracy"
             << ": " << spatial_accuracy_type << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Surface"
             << ": " << write_surface(surface) << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Temporal Method"
             << ": " << write_temporal_method(temporal_method) << endl;
  out_stream << setiosflags(ios::left) << setw(18) << "Spatial Method"
             << ": " << write_spatial_method(spatial_method) << endl
             << endl;
}
