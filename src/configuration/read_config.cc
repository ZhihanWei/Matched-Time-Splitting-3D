#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>

#include "configuration/read_config.h"
#include "helper.h"

using namespace std;

/***********************************************************************
      read configuration parameters from input configuration file
                          Constructor

 INPUT
 file_name : name of txt file
 **********************************************************************/
ReadConfig::ReadConfig(const string &file) {
  ifstream read_stream;

  read_stream.open(file, ios::in);
  if (!read_stream.good()) {
    LOG_FATAL("Read configuration file failed, check path of " + file);
  }

  string line;

  while (getline(read_stream, line)) {
    if (line.find("#") != string::npos || line.empty()) {
      continue;
    }

    string config = line.substr(0, line.find("=") - 1);
    config.erase(remove(config.begin(), config.end(), ' '), config.end());
    string e = line.substr(line.find("=") + 1);
    e.erase(remove(e.begin(), e.end(), ' '), e.end());

    switch (Translate(config)) {
      case MAX_X:
        x_max = stoi(e);
        break;
      case MIN_X:
        x_min = stoi(e);
        break;
      case MAX_Y:
        y_max = stoi(e);
        break;
      case MIN_Y:
        y_min = stoi(e);
        break;
      case MAX_Z:
        z_max = stoi(e);
        break;
      case MIN_Z:
        z_min = stoi(e);
        break;
      case TIME_START:
        t_start = stod(e);
        break;
      case TIME_TERMINATE:
        t_terminate = stod(e);
        break;
      case TIME_STEP:
        t_step = stod(e);
        break;
      case NX:
        nx = stoi(e);
        break;
      case NY:
        ny = stoi(e);
        break;
      case NZ:
        nz = stoi(e);
        break;
      case SURFACE:
        surface = ParseSurface(e);
        break;
      case TEMPORAL_METHOD:
        temporal_method = ParseTemporalMethod(e);
        break;
      case EQUATION:
        equation = ParseEquationMethod(e);
        break;
      case SPATIAL_METHOD:
        spatial_method = ParseSpatialMethod(e);
        break;
      case DIFFUSION_COEFFICIENT:
        diffusion_coef = ParseDiffusionCoefficient(e);
        break;
      case SPATIAL_ACCURACY:
        spatial_accuracy = ParseSpatialAccuracy(e);
        break;
    }
  }

  if ((spatial_accuracy == 4) && (surface != Surface_Type::CUBE)) {
    LOG_FATAL("Currently, 4th order of spatial convergence rate only suppport surface \"cube\", please modify configuration");
  }

  read_stream.close();
}

/***********************************************************************************************
                            Translate a string to a enumeration member
                                        called by Constructor

 INPUT
 arg : string name

 OUTPUT
 enumeration member
 ***********************************************************************************************/
Config ReadConfig::Translate(const string &arg) {
  if (arg == "x-max") {
    return MAX_X;
  } else if (arg == "x-min") {
    return MIN_X;
  } else if (arg == "y-max") {
    return MAX_Y;
  } else if (arg == "y-min") {
    return MIN_Y;
  } else if (arg == "z-max") {
    return MAX_Z;
  } else if (arg == "z-min") {
    return MIN_Z;
  } else if (arg == "time-start") {
    return TIME_START;
  } else if (arg == "time-finish") {
    return TIME_TERMINATE;
  } else if (arg == "time-step") {
    return TIME_STEP;
  } else if (arg == "x-grid") {
    return NX;
  } else if (arg == "y-grid") {
    return NY;
  } else if (arg == "z-grid") {
    return NZ;
  } else if (arg == "surface") {
    return SURFACE;
  } else if (arg == "temporal-method") {
    return TEMPORAL_METHOD;
  } else if (arg == "spatial-method") {
    return SPATIAL_METHOD;
  } else if (arg == "equation") {
    return EQUATION;
  } else if (arg == "beta") {
    return DIFFUSION_COEFFICIENT;
  } else if (arg == "spatial-accuracy") {
    return SPATIAL_ACCURACY;
  } else {
    LOG_FATAL("Bad element \"" + arg + "\" in configuration file, check configuration file example at \"example/config.txt\"");
  }
}

Temporal_Method_Type ReadConfig::ParseTemporalMethod(const string &in_str) {
  if (in_str == "adi") {
    return Temporal_Method_Type::ADI;
  } else if (in_str == "lod-ie") {
    return Temporal_Method_Type::LOD_IE;
  } else if (in_str == "lod-cn") {
    return Temporal_Method_Type::LOD_CN;
  } else if (in_str == "ts") {
    return Temporal_Method_Type::TS;
  } else {
    LOG_FATAL("Temporal-method type not found! Currently support: \"adi\", \"lod-ie\", \"lod-cn\", \"ts\", check configuration file");
  }
}

Surface_Type ReadConfig::ParseSurface(const string &in_str) {
  if (in_str == "tanglecube") {
    return Surface_Type::TANGLECUBE;
  } else if (in_str == "cube") {
    return Surface_Type::CUBE;
  } else if (in_str == "cylinder") {
    return Surface_Type::CYLINDER;
  } else if (in_str == "ellipsoid") {
    return Surface_Type::ELLIPSOID;
  } else if (in_str == "cone") {
    return Surface_Type::CONE;
  } else if (in_str == "pile") {
    return Surface_Type::PILE;
  } else if (in_str == "torus") {
    return Surface_Type::TORUS;
  } else if (in_str == "dupin_cyclide") {
    return Surface_Type::DUPIN_CYCLIDE;
  } else if (in_str == "molecular") {
    return Surface_Type::MOLECULAR;
  } else if (in_str == "heart") {
    return Surface_Type::HEART;
  } else {
    LOG_FATAL(
        "Surface type not found! Currently support: \"tanglecube\", \"cube\", \"cylinder\", \"ellipsoid\", \"cone\", \"pile\", \
        \"torus\", \"dupin_cylide\", \"molecular\", \"heart\", check configuration file");
  }
}

Spatial_Method_Type ReadConfig::ParseSpatialMethod(const string &in_str) {
  if (in_str == "mib-v1") {
    return Spatial_Method_Type::MIB_V1;
  } else if (in_str == "mib-v2") {
    return Spatial_Method_Type::MIB_V2;
  } else {
    LOG_FATAL("Spatial method type not found! Currently support: \"mib-v1\", \"mib-v2\", check configuration file");
  }
}

int ReadConfig::ParseDiffusionCoefficient(const string &in_str) {
  int v = stoi(in_str);
  std::set<int> set_of_diffusion_coeff_type = {0, 1, 2, 3, 4};
  auto search = set_of_diffusion_coeff_type.find(v);
  if (search == set_of_diffusion_coeff_type.end()) {
    LOG_FATAL("Diffusion coefficient type not found! Currently support: 0 - 4, check configuration file");
  }
  return v;
}

int ReadConfig::ParseEquationMethod(const string &in_str) {
  int v = stoi(in_str);
  std::set<int> set_of_diffusion_coeff_type = {0, 1, 2, 3, 4, 5, 6, 7};
  auto search = set_of_diffusion_coeff_type.find(v);
  if (search == set_of_diffusion_coeff_type.end()) {
    LOG_FATAL("Equation type not found! Currently support: 0 - 7, check configuration file");
  }
  return v;
}

int ReadConfig::ParseSpatialAccuracy(const string &in_str) {
  int v = stoi(in_str);
  std::set<int> set_of_diffusion_coeff_type = {2, 4};
  auto search = set_of_diffusion_coeff_type.find(v);
  if (search == set_of_diffusion_coeff_type.end()) {
    LOG_FATAL("Spatial accuracy type not found! Currently support: 2 and 4, check configuration file");
  }
  return v;
}

/***********************************************************************************************
                              Display part of configuration parameters
***********************************************************************************************/
void ReadConfig::Display() {
  cout << "X boundary: [" << x_min << ", " << x_max << "] \n"
       << "Y boundary: [" << y_min << ", " << y_max << "] \n"
       << "Z boundary: [" << z_min << ", " << z_max << "] \n"
       << "Mesh size: [" << nx << ", " << ny << ", " << nz << "] \n"
       << "Time info: [" << t_start << ", " << t_terminate << ", " << t_step << "] \n";
}

VecDoub ReadConfig::GetDomain() const {
  VecDoub domain;

  domain.push_back(x_max);
  domain.push_back(x_min);
  domain.push_back(y_max);
  domain.push_back(y_min);
  domain.push_back(z_max);
  domain.push_back(z_min);

  return domain;
}

VecInt ReadConfig::GetMesh() const {
  VecInt size;

  size.push_back(nx);
  size.push_back(ny);
  size.push_back(nz);

  return size;
}

VecDoub ReadConfig::GetTime() const {
  VecDoub time;

  time.push_back(t_start);
  time.push_back(t_terminate);
  time.push_back(t_step);

  return time;
}

Surface_Type ReadConfig::GetSurface() const {
  return surface;
}

Temporal_Method_Type ReadConfig::GetTemporalMethod() const {
  return temporal_method;
}

int ReadConfig::GetDiffusionCoef() const {
  return diffusion_coef;
}

int ReadConfig::GetSpatialAccuracy() const {
  return spatial_accuracy;
}

int ReadConfig::GetEquation() const {
  return equation;
}

Spatial_Method_Type ReadConfig::GetSpatialmethod() const {
  return spatial_method;
}
