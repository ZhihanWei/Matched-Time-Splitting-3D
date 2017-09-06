#include "Surface_Pile.h"
#include <cmath>
#include <iostream>
#include "Constant.h"

using namespace std;

/****************************************************
 Constructor(default)
 *****************************************************/
Surface_Pile::Surface_Pile() {}

/**********************************************************************
 Implicit function with little perturbation

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of implicit function with perturbation
 ***********************************************************************/
double Surface_Pile::Set_Node(Doub_I x, Doub_I y, Doub_I z) const {
  double diff;

  diff = 5 * pow((x * x + y * y + 5 * z * z - 10), 3) +
         5 * pow((x * x + 60 * y * y + 3 * z * z - 2), 2) + TOL_SETUP;

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
double Surface_Pile::func(const Doub_I x, const Doub_I y,
                          const Doub_I z) const {
  double temp;

  temp = 5 * pow((x * x + y * y + 5 * z * z - 10), 3) +
         5 * pow((x * x + 60 * y * y + 3 * z * z - 2), 2);

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
double Surface_Pile::dfx(const Doub_I x, const Doub_I y, const Doub_I z) const {
  double temp;

  temp = 30 * x * pow((x * x + y * y + 5 * z * z - 10), 2) -
         20 * x * (x * x + 60 * y * y + 3 * z * z - 2);

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
double Surface_Pile::dfy(const Doub_I x, const Doub_I y, const Doub_I z) const {
  double temp;

  temp = 30 * y * pow((x * x + y * y + 5 * z * z - 10), 2) -
         1200 * y * (x * x + 60 * y * y + 3 * z * z - 2);

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
double Surface_Pile::dfz(const Doub_I x, const Doub_I y, const Doub_I z) const {
  double temp;

  temp = 150 * z * pow((x * x + y * y + 5 * z * z - 10), 2) -
         60 * z * (x * x + 60 * y * y + 3 * z * z - 2);

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
VecDoub Surface_Pile::normal(Doub_I x, Doub_I y, Doub_I z) {
  VecDoub normal;
  double total, temp;

  total = sqrt(pow(30 * x * pow((x * x + y * y + 5 * z * z - 10), 2) -
                       20 * x * (x * x + 60 * y * y + 3 * z * z - 2),
                   2) +
               pow(30 * y * pow((x * x + y * y + 5 * z * z - 10), 2) -
                       1200 * y * (x * x + 60 * y * y + 3 * z * z - 2),
                   2) +
               pow(150 * z * pow((x * x + y * y + 5 * z * z - 10), 2) -
                       60 * z * (x * x + 60 * y * y + 3 * z * z - 2),
                   2));

  temp = (30 * x * pow((x * x + y * y + 5 * z * z - 10), 2) -
          20 * x * (x * x + 60 * y * y + 3 * z * z - 2)) /
         total;
  normal.push_back(temp);
  temp = (30 * y * pow((x * x + y * y + 5 * z * z - 10), 2) -
          1200 * y * (x * x + 60 * y * y + 3 * z * z - 2)) /
         total;
  normal.push_back(temp);
  temp = (150 * z * pow((x * x + y * y + 5 * z * z - 10), 2) -
          60 * z * (x * x + 60 * y * y + 3 * z * z - 2)) /
         total;
  normal.push_back(temp);

  return normal;
}

/**********************************************************************
 Show the quation of implicit function
 ***********************************************************************/
void Surface_Pile::display() {
  cout << "The equation of Pile surface is:" << endl;
  cout << "f(x) = 5*(x^2+y^2+5*z^2-10)^3+5*(x^2+60*y^2+3*z^2-2)^2 = 0" << endl;
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
double Surface_Pile::check(Doub_I x, Doub_I y, Doub_I z) {
  double temp;

  temp = 5 * pow((x * x + y * y + 5 * z * z - 10), 3) +
         5 * pow((x * x + 60 * y * y + 3 * z * z - 2), 2);

  return temp;
}