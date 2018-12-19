#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include "configuration/read_config.h"
#include "constant.h"
#include "helper.h"

string write_surface(const Surface_Type s) {
  string r;

  switch (s) {
    case Surface_Type::TANGLECUBE:
      r = "Tanglecube";
      break;
    case Surface_Type::CUBE:
      r = "Cube";
      break;
    case Surface_Type::CYLINDER:
      r = "Cylinder";
      break;
    case Surface_Type::ELLIPSOID:
      r = "Ellipsoid";
      break;
    case Surface_Type::CONE:
      r = "Cone";
      break;
    case Surface_Type::PILE:
      r = "Pile";
      break;
    case Surface_Type::TORUS:
      r = "Torus";
      break;
    case Surface_Type::DUPIN_CYCLIDE:
      r = "Dupin-Cyclide";
      break;
    case Surface_Type::MOLECULAR:
      r = "Molecular";
      break;
    case Surface_Type::HEART:
      r = "Heart";
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
      r = "MIB-V1";
      break;
    case Spatial_Method_Type::MIB_V2:
      r = "MIB-V2";
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

  out_stream << "----------------- Configuration informations -----------------" << endl;
  out_stream << "Mesh Size                          :   "
             << "NX = " << size.at(0) << "; NY = " << size.at(1)
             << "; NZ = " << size.at(2) << endl;
  out_stream << "Time Step                          :   "
             << setprecision(1) << scientific << time_info.at(2) << endl;
  string jump_type = JP == 'r' ? "real" : "approximate";
  out_stream << "Jump                               :   "
             << fixed << jump_type << endl;
  out_stream << "Diffusion Coefficient Type         :   "
             << diffusion_coef_type << endl;
  out_stream << "Equation Type                      :   "
             << equation_type << endl;
  out_stream << "Order of Spatial Accuracy          :   "
             << spatial_accuracy_type << endl;
  out_stream << "Surface                            :   "
             << write_surface(surface) << endl;
  out_stream << "Temporal Method                    :   "
             << write_temporal_method(temporal_method) << endl;
  out_stream << "Spatial Method                     :   "
             << write_spatial_method(spatial_method) << endl
             << endl;
}
