#include "Surface_Molecular.h"
#include <cmath>
#include <iostream>
#include "Constant.h"

using namespace std;

/****************************************************
 Constructor(default)
 *****************************************************/
Surface_Molecular::Surface_Molecular() {}

/**********************************************************************
 Implicit function with little perturbation

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of implicit function with perturbation
 ***********************************************************************/
double Surface_Molecular::Set_Node(Doub_I x, Doub_I y, Doub_I z) const {
  double diff;

  diff = pow((x * x + y * y + z * z + 0.6), 2) - 3.5 * y * y - 0.6 + TOL_SETUP;

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
double Surface_Molecular::func(const Doub_I x, const Doub_I y,
                               const Doub_I z) const {
  double temp;

  temp = pow((x * x + y * y + z * z + 0.6), 2) - 3.5 * y * y - 0.6;

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
double Surface_Molecular::dfx(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  temp = 4 * x * (x * x + y * y + z * z + 0.6);

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
double Surface_Molecular::dfy(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  temp = 4 * y * (x * x + y * y + z * z + 0.6) - 7 * y;

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
double Surface_Molecular::dfz(const Doub_I x, const Doub_I y,
                              const Doub_I z) const {
  double temp;

  temp = 4 * z * (x * x + y * y + z * z + 0.6);

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
VecDoub Surface_Molecular::normal(Doub_I x, Doub_I y, Doub_I z) {
  VecDoub normal;
  double total, temp;

  total = sqrt(pow(4 * x * (x * x + y * y + z * z + 0.6), 2) +
               pow((4 * y * (x * x + y * y + z * z + 0.6) - 7 * y), 2) +
               pow(4 * z * (x * x + y * y + z * z + 0.6), 2));

  temp = (4 * x * (x * x + y * y + z * z + 0.6)) / total;
  normal.push_back(temp);
  temp = (4 * y * (x * x + y * y + z * z + 0.6) - 7 * y) / total;
  normal.push_back(temp);
  temp = (4 * z * (x * x + y * y + z * z + 0.6)) / total;
  normal.push_back(temp);

  return normal;
}

/**********************************************************************
 Show the quation of implicit function
 ***********************************************************************/
void Surface_Molecular::display() {
  cout << "The equation of Molecular surface is:" << endl;
  cout << "f(x) = (x^2+y^2+5*z^2+0.6)^2-3.5*y*y-0.6 = 0" << endl;
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
double Surface_Molecular::check(Doub_I x, Doub_I y, Doub_I z) {
  double temp;

  temp = pow((x * x + y * y + z * z + 0.6), 2) - 3.5 * y * y - 0.6;

  return temp;
}