#include <fstream>
#include <iostream>

#include "data/data.h"
#include "helper.h"

using namespace std;

/***********************************************************************
      read configuration parameters from input configuration file
                          Constructor

 INPUT
 file_name : name of txt file
 **********************************************************************/
Data::Data(const string &file) {
  ifstream read_stream;

  read_stream.open(file_name, ios::in);
  if (!read_stream.good()) {
    LOG_FATAL("read configuration file failed, check path of " + file);
  }

  char tmp;  // for '=' in configuration file
  string config;

  while (read_stream >> config >> tmp) {
    string e;

    switch (Translate(config)) {
      case MAX_X:
        read_stream >> e;
        x_max = atoi(e);
        break;
      case MIN_X:
        read_stream >> e;
        x_min = atoi(e);
        break;
      case MAX_Y:
        read_stream >> e;
        y_max = atoi(e);
        break;
      case MIN_Y:
        read_stream >> e;
        y_min = atoi(e);
        break;
      case MAX_Z:
        read_stream >> e;
        z_max = atoi(e);
        break;
      case MIN_Z:
        read_stream >> e;
        z_min = atoi(e);
        break;
      case TIME_START:
        read_stream >> e;
        t_start = atod(e);
        break;
      case TIME_TERMINATE:
        read_stream >> e;
        t_terminate = atod(e);
        break;
      case TIME_STEP:
        read_stream >> e;
        t_step = atod(e);
        break;
      case NX:
        read_stream >> e;
        nx = atoi(e);
        break;
      case NY:
        read_stream >> e;
        ny = atoi(e);
        break;
      case NZ:
        read_stream >> e;
        nz = atoi(e);
        break;
      case SURFACE:
        read_stream >> e;
        surface = ParseSurface(e);
        break;
      case TEMPORAL_METHOD:
        read_stream >> e;
        temporal_method = ParseTemporalMethod(e);
        break;
      case EQUATION:
        read_stream >> e;
        equation = ParseEquationMethod(e);
        break;
      case SPATIAL_METHOD:
        read_stream >> e;
        spatial_method = ParseTemporalMethod(e);
        break;
      case DIFFUSION_COEFFICIENT:
        read_stream >> e;
        diffusion_coef = ParseDiffusionCoefficient(e);
        break;
      case SPATIAL_ACCURACY:
        read_stream >> e;
        accuracy = ParseSpatialAccuracy(e);
        break;
      case COMMENT:
        break;
    }
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
Data::Config Data::Translate(const string &arg) {
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
  } else if (arg.find("#") == 0) {
    return COMMENT;
  } else {
    LOG_FATAL("bad element in configuration file, check configuration file example at \"example/config.txt\"");
  }
}

string ParseTemporalMethod(const string &in_str) {
  return in_str;
}

string ParseSurface(const string &in_str) {
  return in_str;
}

string ParseSpatialMethod(const string &in_str) {
  return in_str;
}

int ParseDiffusionCoefficient(const string &in_str) {
  return atoi(in_str);
}

int ParseEquationMethod(const string &in_str) {
  return atoi(in_str);
}

int ParseSpatialMethod(const string &in_str) {
  return atoi(in_str);
}

/***********************************************************************************************
                              Display part of configuration parameters
***********************************************************************************************/
void Data::Display(){
    cout << "X boundary: [" << x_min << ", " << x_max << "] \n"
         << "Y boundary: [" << y_min << ", " << y_max << "] \n"
         << "Z boundary: [" << z_min << ", " << z_max << "] \n"
         << "Mesh size: [" << nx << ", " << ny << ", " << nz << "] \n"
         << "Time info: [" << t_start << ", " << t_terminate << ", " << t_step << "] \n"}

VecDoub Data::GetDomain() const {
  VecDoub domain;

  domain.push_back(x_max);
  domain.push_back(x_min);
  domain.push_back(y_max);
  domain.push_back(y_min);
  domain.push_back(z_max);
  domain.push_back(z_min);

  return domain;
}

VecInt Data::GetMesh() const {
  VecInt size;

  size.push_back(nx);
  size.push_back(ny);
  size.push_back(nz);

  return size;
}

VecDoub Data::GetTime() const {
  VecDoub time;

  time.push_back(t_start);
  time.push_back(t_terminate);
  time.push_back(t_step);

  return time;
}

string &Data::GetSurface() const {
  return surface;
}

string &Data::GetTemporalMethod() const {
  return temporal_method;
}

int Data::GetDiffusionCoeff() const {
  return diffusion_coeff;
}

int Data::GetSpatialAccuracy() const {
  return accuracy;
}

int Data::GetEquation() const {
  return equation;
}

string &Data::GetSpatialmethod() const {
  return spatial_method;
}
