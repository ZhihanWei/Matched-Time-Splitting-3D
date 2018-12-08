#include <iostream>
#include "surface/surface_cylinder.h"

using namespace std;
/****************************************************
 Constructor

 INPUT
 arg : arguments of Cube object
 *****************************************************/
Surface_Cylinder::Surface_Cylinder(VecDoub_I arg) {
  r = arg[0];
  z_lower = arg[1];
  z_upper = arg[2];
}

/**********************************************************************
 Implicit function with little perturbation

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 positive when outside; negative when inside
 ***********************************************************************/
double Surface_Cylinder::Set_Node(Doub_I x, Doub_I y, Doub_I z) const {
  double diff;

  if ((z > z_lower) && (z < z_upper)) {
    if ((x * x + y * y) < r * r) {
      diff = -1 + TOL_SETUP;
    } else if ((x * x + y * y) == r * r) {
      cout << "Node is on grid" << endl;
      exit(0);
    } else {
      diff = 1 + TOL_SETUP;
    }
  } else if ((z == z_lower) || (z == z_upper)) {
    cout << "Node is on grid" << endl;
    exit(0);
  } else {
    diff = 1 + TOL_SETUP;
  }

  return diff;
}

/*********************************************************************************
 Private member function to get value of implicit function without perturbation

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of function
 **********************************************************************************/
double Surface_Cylinder::func(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  if ((z > z_lower) && (z < z_upper)) {
    temp = x * x + y * y - r * r;
  } else {
    cout << "Function is not defined!" << endl;
    exit(0);
  }

  return temp;
}

/**********************************************************************
 Calculate derivative of function in x-direction

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 derivative of function in x-direction
 ***********************************************************************/
double Surface_Cylinder::dfx(const Doub_I x, const Doub_I y,
                             const Doub_I z) const {
  double temp;

  if ((z > z_lower) && (z < z_upper)) {
    temp = 2 * x;
  } else {
    cout << "Function is not defined!" << endl;
    exit(0);
  }

  return temp;
}

/**********************************************************************
 Calculate derivative of function in y-direction

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 derivative of function in y-direction
 ***********************************************************************/
double Surface_Cylinder::dfy(const Doub_I x, const Doub_I y,
                             const Doub_I z) const {
  double temp;

  if ((z > z_lower) && (z < z_upper)) {
    temp = 2 * y;
  } else {
    cout << "Function is not defined!" << endl;
    exit(0);
  }

  return temp;
}

/***********************************************************************************************
 Calculate the on-grid interface point, bd1 can be either smaller or greater
 than bd2

 INPUT
 x  : x coordinate of the given point
 y  : y coordinate of the given point
 z1 : z coordinate of ond boundary node for the given point
 z2 : z coordinate of the other boundary node for the given point

 OUTPUT
 z coordinate of on-grid intersection point
 ***********************************************************************************************/
double Surface_Cylinder::Gamma_z(Doub_I x, Doub_I y, Doub_I z1, Doub_I z2) {
  double temp, up_bd, lw_bd;

  if (z1 > z2) {
    up_bd = z1;
    lw_bd = z2;
  } else if (z1 < z2) {
    up_bd = z2;
    lw_bd = z1;
  } else {
    cout << "Two bounds are too close!";
    exit(0);
  }

  if ((z_lower > lw_bd) && (z_lower < up_bd)) {
    temp = z_lower;
  } else if ((z_upper > lw_bd) && (z_upper < up_bd)) {
    temp = z_upper;
  } else {
    cout << "No interface is found in y-direction" << endl;
  }

  return temp;
}

/**********************************************************************
 Get normal direction at given point

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 vector of normal direction with 3 components
 ***********************************************************************/
VecDoub Surface_Cylinder::normal(Doub_I x, Doub_I y, Doub_I z) {
  VecDoub normal;
  double total, temp;

  if (z == z_lower) {
    normal.push_back(0);
    normal.push_back(0);
    normal.push_back(-1);
  } else if (z == z_upper) {
    normal.push_back(0);
    normal.push_back(0);
    normal.push_back(1);
  } else {
    total = sqrt(pow(2 * x, 2) + pow(2 * y, 2));

    temp = 2 * x / total;
    normal.push_back(temp);
    temp = 2 * y / total;
    normal.push_back(temp);
    normal.push_back(0);
  }

  return normal;
}

/**********************************************************************
 Show the quation of implicit function
 ***********************************************************************/
void Surface_Cylinder::display() {
  cout << "The equation of Cylinder surface is:" << endl;
  cout << "f(x) = x^2 + y^2 - " << r << "^2 = 0, where z in [" << z_lower << ","
       << z_upper << "]" << endl;
}

/*****************************************************************************************
 Public member function to get the value of implicit function without
 perturbation

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of implicit function without perturbation
 ******************************************************************************************/
double Surface_Cylinder::check(Doub_I x, Doub_I y, Doub_I z) {
  double temp;

  if ((z > z_lower) && (z < z_upper)) {
    temp = x * x + y * y - r * r;
  } else {
    temp = 1;
  }

  return temp;
}

/******************************************************************************
 The following member function should not be defined for cube interface
 *******************************************************************************/
double Surface_Cylinder::dfz(const Doub_I x, const Doub_I y,
                             const Doub_I z) const {
  cout << "No implicit function defined for cylinder surface" << endl;
  exit(0);

  return (0.0);
}
