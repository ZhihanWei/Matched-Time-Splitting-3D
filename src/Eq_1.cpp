#include "Eq_1.h"
#include <cmath>

/******************************************************************
                            Constructor

 INPUT
 ti   : current time
 beta : vector of 2 double values represent beta^{-} and beta^{+}
 ******************************************************************/
Eq_1::Eq_1(Doub_I ti, Beta& in_beta) : beta(in_beta) { t = ti; }

/****************************************************
 Analytical solution in Omega^{-}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of analytical solution at given node
 *****************************************************/
double Eq_1::Inner_u(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = exp(0.5*z)*exp(0.25*t);
  temp = exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/****************************************************
 Analytical solution in Omega^{+}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of analytical solution at given node
 *****************************************************/
double Eq_1::Outer_u(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = exp(0.5*z)*exp(0.25*t);
  temp = exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/****************************************************
 Source term in Omega^{+}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of source term at given node
 *****************************************************/
double Eq_1::Inner_f(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0.75 * Inner_u(x, y, z) -
         beta.Inside(x, y, z) * 0.75 * Inner_u(x, y, z) -
         beta.Inside_Dx(x, y, z) * Inner_dux(x, y, z) -
         beta.Inside_Dy(x, y, z) * Inner_duy(x, y, z) -
         beta.Inside_Dz(x, y, z) * Inner_duz(x, y, z);

  return temp;
}

/****************************************************
 Source term in Omega^{+}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of source term at given node
 *****************************************************/
double Eq_1::Outer_f(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0.75 * Outer_u(x, y, z) -
         beta.Outside(x, y, z) * 0.75 * Outer_u(x, y, z) -
         beta.Outside_Dx(x, y, z) * Outer_dux(x, y, z) -
         beta.Outside_Dy(x, y, z) * Outer_duy(x, y, z) -
         beta.Outside_Dz(x, y, z) * Outer_duz(x, y, z);

  return temp;
}

/**********************************************************************
 1st derivative in x-direction of analytical solution in Omega^{-}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of 1st derivative in x-direction at given node
 ***********************************************************************/
double Eq_1::Inner_dux(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0;
  temp = 0.5 * exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/**********************************************************************
 1st derivative in x-direction of analytical solution in Omega^{+}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of 1st derivative in x-direction at given node
 ***********************************************************************/
double Eq_1::Outer_dux(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0;
  temp = 0.5 * exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);
  ;

  return temp;
}

/**********************************************************************
 1st derivative in y-direction of analytical solution in Omega^{-}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of 1st derivative in y-direction at given node
 ***********************************************************************/
double Eq_1::Inner_duy(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0;
  temp = 0.5 * exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/**********************************************************************
 1st derivative in y-direction of analytical solution in Omega^{+}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of 1st derivative in y-direction at given node
 ***********************************************************************/
double Eq_1::Outer_duy(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0;
  temp = 0.5 * exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/**********************************************************************
 1st derivative in z-direction of analytical solution in Omega^{-}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of 1st derivative in z-direction at given node
 ***********************************************************************/
double Eq_1::Inner_duz(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0.5*exp(0.5*z)*exp(0.25*t);
  temp = 0.5 * exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/**********************************************************************
 1st derivative in z-direction of analytical solution in Omega^{+}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 value of 1st derivative in z-direction at given node
 ***********************************************************************/
double Eq_1::Outer_duz(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0.5*exp(0.5*z)*exp(0.25*t);
  temp = 0.5 * exp(0.5 * x) * exp(0.5 * y) * exp(0.5 * z) * exp(0.75 * t);

  return temp;
}

/***************************************************************
 Analytical jump [beta_u] in x-direction

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 analutical of jump [beta_u_x]
 ***************************************************************/
double Eq_1::Jump_betau_x(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  // temp = 0;
  temp = beta.Outside(x, y, z) * Outer_dux(x, y, z) -
         beta.Inside(x, y, z) * Inner_dux(x, y, z);

  return temp;
}

/***************************************************************
 Analytical jump [beta_u] in y-direction

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 analutical of jump [beta_u_y]
 ***************************************************************/
double Eq_1::Jump_betau_y(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = beta.Outside(x, y, z) * Outer_duy(x, y, z) -
         beta.Inside(x, y, z) * Inner_duy(x, y, z);

  return temp;
}

/***************************************************************
 Analytical jump [beta_u] in z-direction

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 analutical of jump [beta_u_z]
 ***************************************************************/
double Eq_1::Jump_betau_z(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = beta.Outside(x, y, z) * Outer_duz(x, y, z) -
         beta.Inside(x, y, z) * Inner_duz(x, y, z);

  return temp;
}