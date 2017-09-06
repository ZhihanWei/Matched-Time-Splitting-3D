#include "Surface_Ellipsoid.h"
#include <cmath>
#include <iostream>
#include "Constant.h"

using namespace std;

/****************************************************
                       Constructor

 INPUT
 arg : arguments of Ellipsoid object
 *****************************************************/
Surface_Ellipsoid::Surface_Ellipsoid(VecDoub_I arg) {
  a = arg[0];
  b = arg[1];
  c = arg[2];
  radius = arg[3];
}

/**********************************************************************
                   Implicit function with little perturbation

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of implicit function with perturbation
 ***********************************************************************/
double Surface_Ellipsoid::Set_Node(Doub_I x, Doub_I y, Doub_I z) const {
  double diff;

  diff = x * x / (a * a) + y * y / (b * b) + z * z / (c * c) - radius * radius +
         TOL_SETUP;

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
double Surface_Ellipsoid::func(const Doub_I x, const Doub_I y,
                               const Doub_I z) const {
  double temp;

  temp = x * x / (a * a) + y * y / (b * b) + z * z / (c * c) - radius * radius;

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
double Surface_Ellipsoid::dfx(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  temp = 2 * x / (a * a);

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
double Surface_Ellipsoid::dfy(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  temp = 2 * y / (b * b);

  return temp;
}

/**********************************************************************
 Calculate derivative of function in z-direction

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 derivative of function in z-direction
 ***********************************************************************/
double Surface_Ellipsoid::dfz(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  temp = 2 * z / (c * c);

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
VecDoub Surface_Ellipsoid::normal(Doub_I x, Doub_I y, Doub_I z) {
  VecDoub normal;
  double total, temp;

  total = sqrt(pow(2 * x / (a * a), 2) + pow(2 * y / (b * b), 2) +
               pow(2 * z / (c * c), 2));

  temp = 2 * x / (a * a) / total;
  normal.push_back(temp);
  temp = 2 * y / (b * b) / total;
  normal.push_back(temp);
  temp = 2 * z / (c * c) / total;
  normal.push_back(temp);

  return normal;
}

/**********************************************************************
 Show the quation of implicit function
 ***********************************************************************/
void Surface_Ellipsoid::display() {
  cout << "The equation of Ellipsoid surface is:" << endl;
  cout << "f(x) = 1/" << a * a << "*x^2+"
       << "1/" << b * b << "*y^2+"
       << "1/" << c * c << "*z^2-(" << radius << ")^2 = 0" << endl;
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
double Surface_Ellipsoid::check(Doub_I x, Doub_I y, Doub_I z) {
  double temp;

  temp = x * x / (a * a) + y * y / (b * b) + z * z / (c * c) - radius * radius;

  return temp;
}