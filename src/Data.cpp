#include "Data.h"
#include <fstream>
#include <iostream>

using namespace std;

/*********************************************************************
                        read data from data.txt file
                                Constructor

 INPUT
 file_name : name of txt file
 **********************************************************************/
Data::Data(const string& file_name) {
  ifstream read_file;

  read_file.open(file_name, ios::in);
  if (read_file.good()) {
    char tmp;
    string parameter_name;

    while (read_file >> parameter_name >> tmp) {
      string a;
      int parameter_value_1;
      double parameter_value_2;
      char parameter_value_3;

      switch (Translation(parameter_name)) {
        case e_xl:
          read_file >> parameter_value_2;
          xl = parameter_value_2;
          break;
        case e_xr:
          read_file >> parameter_value_2;
          xr = parameter_value_2;
          break;
        case e_yl:
          read_file >> parameter_value_2;
          yl = parameter_value_2;
          break;
        case e_yr:
          read_file >> parameter_value_2;
          yr = parameter_value_2;
          break;
        case e_zl:
          read_file >> parameter_value_2;
          zl = parameter_value_2;
          break;
        case e_zr:
          read_file >> parameter_value_2;
          zr = parameter_value_2;
          break;
        case e_t_start:
          read_file >> parameter_value_2;
          t_start = parameter_value_2;
          break;
        case e_t_finish:
          read_file >> parameter_value_2;
          t_finish = parameter_value_2;
          break;
        case e_t_step:
          read_file >> parameter_value_2;
          t_step = parameter_value_2;
          break;
        case e_nx:
          read_file >> parameter_value_1;
          nx = parameter_value_1;
          break;
        case e_ny:
          read_file >> parameter_value_1;
          ny = parameter_value_1;
          break;
        case e_nz:
          read_file >> parameter_value_1;
          nz = parameter_value_1;
          break;
        case e_surface:
          read_file >> parameter_value_3;
          surface = parameter_value_3;
          break;
        case e_method:
          read_file >> parameter_value_3;
          method = parameter_value_3;
          break;
        case e_equation:
          read_file >> parameter_value_1;
          equation = parameter_value_1;
          break;
        case e_mib_method:
          read_file >> parameter_value_1;
          mib_method = parameter_value_1;
          break;
        case e_beta:
          read_file >> parameter_value_1;
          beta = parameter_value_1;
          break;
        case e_accuracy:
          read_file >> parameter_value_1;
          accuracy = parameter_value_1;
          break;
          if (accuracy % 2 == 1) {
            cout << "Order of accuracy has to be an even order" << endl;
            exit(0);
          }
          break;
        default:
          cout << "No corresponding parameter found, check data.txt file"
               << endl;
          exit(0);
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
Data::Prt_name Data::Translation(const string& in_string) {
  if (in_string == "xl") {
    return e_xl;
  } else if (in_string == "xr") {
    return e_xr;
  } else if (in_string == "yl") {
    return e_yl;
  } else if (in_string == "yr") {
    return e_yr;
  } else if (in_string == "zl") {
    return e_zl;
  } else if (in_string == "zr") {
    return e_zr;
  } else if (in_string == "t_start") {
    return e_t_start;
  } else if (in_string == "t_finish") {
    return e_t_finish;
  } else if (in_string == "t_step") {
    return e_t_step;
  } else if (in_string == "nx") {
    return e_nx;
  } else if (in_string == "ny") {
    return e_ny;
  } else if (in_string == "nz") {
    return e_nz;
  } else if (in_string == "surface") {
    return e_surface;
  } else if (in_string == "method") {
    return e_method;
  } else if (in_string == "mib_method") {
    return e_mib_method;
  } else if (in_string == "equation") {
    return e_equation;
  } else if (in_string == "beta") {
    return e_beta;
  } else if (in_string == "accuracy") {
    return e_accuracy;
  } else {
    cout << "No corresponding parameter found, check data.txt file" << endl;
    exit(0);
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
