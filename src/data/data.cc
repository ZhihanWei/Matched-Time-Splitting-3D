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
Data::Data(const string &file_name) {
  ifstream read_file;

  read_file.open(file_name, ios::in);
  if (read_file.good()) {
    char tmp;  // for '=' in configuration file
    string config;

    while (read_file >> config >> tmp) {
      int parameter_value_1;
      double parameter_value_2;
      char parameter_value_3;

      switch (translate(config)) {
        case MAX_X:
          read_file >> parameter_value_2;
          xl = parameter_value_2;
          break;
        case MIN_X:
          read_file >> parameter_value_2;
          xr = parameter_value_2;
          break;
        case MAX_Y:
          read_file >> parameter_value_2;
          yl = parameter_value_2;
          break;
        case MIN_Y:
          read_file >> parameter_value_2;
          yr = parameter_value_2;
          break;
        case MAX_Z:
          read_file >> parameter_value_2;
          zl = parameter_value_2;
          break;
        case MIN_Z:
          read_file >> parameter_value_2;
          zr = parameter_value_2;
          break;
        case TIME_START:
          read_file >> parameter_value_2;
          t_start = parameter_value_2;
          break;
        case TIME_TERMINATE:
          read_file >> parameter_value_2;
          t_finish = parameter_value_2;
          break;
        case TIME_STEP:
          read_file >> parameter_value_2;
          t_step = parameter_value_2;
          break;
        case NX:
          read_file >> parameter_value_1;
          nx = parameter_value_1;
          break;
        case NY:
          read_file >> parameter_value_1;
          ny = parameter_value_1;
          break;
        case NZ:
          read_file >> parameter_value_1;
          nz = parameter_value_1;
          break;
        case SURFACE:
          read_file >> parameter_value_3;
          surface = parameter_value_3;
          break;
        case TEMPORAL_METHOD:
          read_file >> parameter_value_3;
          method = parameter_value_3;
          break;
        case EQUATION:
          read_file >> parameter_value_1;
          equation = parameter_value_1;
          break;
        case SPATIAL_METHOD:
          read_file >> parameter_value_1;
          mib_method = parameter_value_1;
          break;
        case DIFFUSION_COEFFICIENT:
          read_file >> parameter_value_1;
          beta = parameter_value_1;
          break;
        case SPATIAL_ACCURACY:
          read_file >> parameter_value_1;
          accuracy = parameter_value_1;
          break;
          if (accuracy % 2 == 1) {
            LOG_FATAL("[ERROR]: currently only support 2nd or 4th order of spatial accuracy");
          }
          break;
        case COMMENT:
          break;
        default:
          LOG_FATAL("[ERROR]: bad element in configuration file, check configuration file example at \"example/config.txt\"");
      }
    }
  } else {
    cout << "Data file is not opened!";
    exit(0);
  }
  read_file.close();
}

/***********************************************************************************************
                            Translate a string to a enumeration member
                                        called by Constructor

 INPUT
 in_string : string name

 OUTPUT
 enumeration member
 ***********************************************************************************************/
Data::Config Data::translate(const string &in_string) {
  if (in_string == "x-max") {
    return MAX_X;
  } else if (in_string == "x-min") {
    return MIN_X;
  } else if (in_string == "y-max") {
    return MAX_Y;
  } else if (in_string == "y-min") {
    return MIN_Y;
  } else if (in_string == "z-max") {
    return MAX_Z;
  } else if (in_string == "z-min") {
    return MIN_Z;
  } else if (in_string == "time-start") {
    return TIME_START;
  } else if (in_string == "time-finish") {
    return TIME_TERMINATE;
  } else if (in_string == "time-step") {
    return TIME_STEP;
  } else if (in_string == "x-grid") {
    return NX;
  } else if (in_string == "y-grid") {
    return NY;
  } else if (in_string == "z-grid") {
    return NZ;
  } else if (in_string == "surface") {
    return SURFACE;
  } else if (in_string == "temporal-method") {
    return TEMPORAL_METHOD;
  } else if (in_string == "spatial-method") {
    return SPATIAL_METHOD;
  } else if (in_string == "equation") {
    return EQUATION;
  } else if (in_string == "beta") {
    return DIFFUSION_COEFFICIENT;
  } else if (in_string == "spatial-accuracy") {
    return SPATIAL_ACCURACY;
  } else if (in_string.find("#") == 0) {
    return COMMENT;
  } else {
    LOG_FATAL("[ERROR]: bad element in configuration file, check configuration file example at \"example/config.txt\"");
  }
}

/***********************************************************************************************
                                    Print data
 ***********************************************************************************************/
void Data::Display() {
  cout << "xl = " << xl << " xr = " << xr << " yl = " << yl << " yr = " << yr
       << " zl = " << zl << " zr = " << zr << endl;
  cout << "t_start = " << t_start << " t_finish = " << t_finish
       << " t_step = " << t_step << endl;
  cout << "nx = " << nx << " ny = " << ny << " nz = " << nz << endl;
}

/***********************************************************************************************
                        Get data given by order [xl,xr,yl,yr,zl,zr]
 ***********************************************************************************************/
VecDoub Data::Get_Domain() const {
  VecDoub domain;

  domain.push_back(xl);
  domain.push_back(xr);
  domain.push_back(yl);
  domain.push_back(yr);
  domain.push_back(zl);
  domain.push_back(zr);

  return domain;
}

/***********************************************************************************************
                                Get data given by order [nx,ny,nz]
 ***********************************************************************************************/
VecInt Data::Get_Size() const {
  VecInt size;

  size.push_back(nx);
  size.push_back(ny);
  size.push_back(nz);

  return size;
}

/***********************************************************************************************
                            Get data given by order [t_start,t_finish,t_step]
 ***********************************************************************************************/
VecDoub Data::Get_Time() const {
  VecDoub time;

  time.push_back(t_start);
  time.push_back(t_finish);
  time.push_back(t_step);

  return time;
}

/***********************************************************************************************
                                    Get data Surface
 ***********************************************************************************************/
char Data::Get_Surface() const { return surface; }

/***********************************************************************************************
                                    Get data Method
 ***********************************************************************************************/
char Data::Get_Method() const { return method; }

/***********************************************************************************************
                                    Get data Beta
 ***********************************************************************************************/
int Data::Get_Beta() const { return beta; }

/***********************************************************************************************
                                   Get data Accuracy
 ***********************************************************************************************/
int Data::Get_Accuracy() const { return accuracy; }

/***********************************************************************************************
                                    Get data Equation
 ***********************************************************************************************/
int Data::Get_Equation() const { return equation; }

/***********************************************************************************************
                                    Get data MIB method
 ***********************************************************************************************/
int Data::Get_MIB_method() const { return mib_method; }
