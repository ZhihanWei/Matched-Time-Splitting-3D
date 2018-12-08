#include "equation/equation.h"

/****************************************************
            Analytical jump [u]

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 analutical of jump [u]
 *****************************************************/
double Equation::Jump_u(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = Outer_u(x, y, z) - Inner_u(x, y, z);

  return temp;
}

/***************************************************************
 Analytical jump [beta_u_xi]

 INPUT
 x  : x coordinate of given point
 y  : y coordinate of given point
 z  : z coordinate of given point
 vx : x direction component of transformation vector
 vy : y direction component of transformation vector
 vz : z direction component of transformation vector

 OUTPUT
 analutical of jump [beta_u_xi]
 ***************************************************************/
double Equation::Jump_betau_xi(Doub_I x, Doub_I y, Doub_I z, Doub_I vx,
                               Doub_I vy, Doub_I vz) const {
  double temp;

  temp = vx * Jump_betau_x(x, y, z) + vy * Jump_betau_y(x, y, z) +
         vz * Jump_betau_z(x, y, z);

  return temp;
}

/***************************************************************
 Analytical jump [u_eta]

 INPUT
 x  : x coordinate of given point
 y  : y coordinate of given point
 z  : z coordinate of given point
 vx : x direction component of transformation vector
 vy : y direction component of transformation vector
 vz : z direction component of transformation vector

 OUTPUT
 analutical of jump [u_eta]
 ***************************************************************/
double Equation::Jump_u_eta(Doub_I x, Doub_I y, Doub_I z, Doub_I vx, Doub_I vy,
                            Doub_I vz) const {
  double temp;

  temp = vx * (Outer_dux(x, y, z) - Inner_dux(x, y, z)) +
         vy * (Outer_duy(x, y, z) - Inner_duy(x, y, z)) +
         vz * (Outer_duz(x, y, z) - Inner_duz(x, y, z));

  return temp;
}

/***************************************************************
 Analytical jump [u_tau]

 INPUT
 x  : x coordinate of given point
 y  : y coordinate of given point
 z  : z coordinate of given point
 vx : x direction component of transformation vector
 vy : y direction component of transformation vector
 vz : z direction component of transformation vector

 OUTPUT
 analutical of jump [u_tau]
 ***************************************************************/
double Equation::Jump_u_tau(Doub_I x, Doub_I y, Doub_I z, Doub_I vx, Doub_I vy,
                            Doub_I vz) const {
  double temp;

  temp = vx * (Outer_dux(x, y, z) - Inner_dux(x, y, z)) +
         vy * (Outer_duy(x, y, z) - Inner_duy(x, y, z)) +
         vz * (Outer_duz(x, y, z) - Inner_duz(x, y, z));

  return temp;
}

/***************************************************************
 Analytical jump [beta_u_eta]

 INPUT
 x  : x coordinate of given point
 y  : y coordinate of given point
 z  : z coordinate of given point
 vx : x direction component of transformation vector
 vy : y direction component of transformation vector
 vz : z direction component of transformation vector

 OUTPUT
 analutical of jump [beta_u_eta]
 ***************************************************************/
double Equation::Jump_betau_eta(Doub_I x, Doub_I y, Doub_I z, Doub_I vx,
                                Doub_I vy, Doub_I vz) const {
  double temp;

  temp = vx * Jump_betau_x(x, y, z) + vy * Jump_betau_y(x, y, z) +
         vz * Jump_betau_z(x, y, z);

  return temp;
}

/***************************************************************
 Analytical jump [beta_u_tau]

 INPUT
 x  : x coordinate of given point
 y  : y coordinate of given point
 z  : z coordinate of given point
 vx : x direction component of transformation vector
 vy : y direction component of transformation vector
 vz : z direction component of transformation vector

 OUTPUT
 analutical of jump [beta_u_tau]
 ***************************************************************/
double Equation::Jump_betau_tau(Doub_I x, Doub_I y, Doub_I z, Doub_I vx,
                                Doub_I vy, Doub_I vz) const {
  double temp;

  temp = vx * Jump_betau_x(x, y, z) + vy * Jump_betau_y(x, y, z) +
         vz * Jump_betau_z(x, y, z);

  return temp;
}
