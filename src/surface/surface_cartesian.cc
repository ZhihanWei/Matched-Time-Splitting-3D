#include <cmath>
#include <iostream>

#include "surface/surface_cartesian.h"
#include "constant.h"

using namespace std;

/***********************************************************************************************
 Calculate the on-grid interface point, bd1 can be either smaller or greater
 than bd2

 INPUT
 x1 : x coordinate of one boundary node for the given point
 x2 : x coordinate of the other boundary node for the given point
 y  : y coordinate of the given point
 z  : z coordinate of the given point

 OUTPUT
 x coordinate of on-grid intersection point
 ***********************************************************************************************/
double Surface_Cartesian::Gamma_x(Doub_I x1, Doub_I x2, Doub_I y, Doub_I z) {
  double temp, up_bd, lw_bd;

  if (x1 > x2) {
    up_bd = x1;
    lw_bd = x2;
  } else if (x1 < x2) {
    up_bd = x2;
    lw_bd = x1;
  } else {
    cout << "Two bounds are too close!";
    exit(0);
  }

  if (abs(x1 - x2) < TOL_RTSAFE) {
    temp = (x1 + x2) / 2;
  } else {
    temp = Rtsafe_x(x1, x2, y, z);
  }

  return temp;
}

/***********************************************************************************************
 Calculate the on-grid interface point, bd1 can be either smaller or greater
 than bd2

 INPUT
 x  : x coordinate of the given point
 y1 : y coordinate of one boundary node for the given point
 y2 : y coordinate of the other boundary node for the given point
 z  : z coordinate of the given point

 OUTPUT
 y coordinate of on-grid intersection point
 ***********************************************************************************************/
double Surface_Cartesian::Gamma_y(Doub_I x, Doub_I y1, Doub_I y2, Doub_I z) {
  double temp, up_bd, lw_bd;

  if (y1 > y2) {
    up_bd = y1;
    lw_bd = y2;
  } else if (y1 < y2) {
    up_bd = y2;
    lw_bd = y1;
  } else {
    cout << "Two bounds are too close!";
    exit(0);
  }

  if (abs(y1 - y2) < TOL_RTSAFE) {
    temp = (y1 + y2) / 2;
  } else {
    temp = Rtsafe_y(x, y1, y2, z);
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
double Surface_Cartesian::Gamma_z(Doub_I x, Doub_I y, Doub_I z1, Doub_I z2) {
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

  if (abs(z1 - z2) < TOL_RTSAFE) {
    temp = (z1 + z2) / 2;
  } else {
    temp = Rtsafe_z(x, y, z1, z2);
  }

  return temp;
}

/*******************************************************************************************************************
 Using a combination of Newton-Raphson and bisection, return the root of a
 function bracketed between x1 and x2.
 The root will be refined until its accuracy is known within xacc. funcd is a
 user-supplied struct that returns
 the function value as a functor and the first derivative of the function at the
 point x as the function df

 INPUT
 x1 : x coordinate of left node for the given point
 x2 : x coordinate of right node for the given point
 y  : y coordinate of the given point
 z  : z coordinate of the given point

 OUTPUT
 root of function for x-direction
 *******************************************************************************************************************/
double Surface_Cartesian::Rtsafe_x(const Doub_I x1, const Doub_I x2,
                                   const Doub_I y, const Doub_I z) {
  const int MAXIT = 100; // Maximum allowed number of iterations
  double fl, fh, xh, xl, rts, dxold, dx, f, df, temp;
  int j;

  fl = func(x1, y, z);
  fh = func(x2, y, z);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    throw("Root must be bracketed in rtsafe");
  }
  if (fl == 0.0) {
    return x1;
  }
  if (fh == 0.0) {
    return x2;
  }
  if (fl < 0.0) // Orient the search so that f(x1)<0
  {
    xl = x1;
    xh = x2;
  } else {
    xh = x1;
    xl = x2;
  }
  rts = 0.5 * (x1 + x2); // Initialize the guess for root
  dxold = abs(x2 - x1);  // the "stepsize before last"
  dx = dxold;            // and the last step
  f = func(rts, y, z);
  df = dfx(rts, y, z);
  // Bisect if Newton out of range or not decreasing fast enough
  for (j = 0; j < MAXIT; j++) // Loop over allowed iterations
  {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
        (abs(2.0 * f) > abs(dxold * df))) {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts) // Change in root is negligible
      {
        return rts;
      }
    } else // Newton step acceptable, take it
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts) {
        return rts;
      }
    }
    if (abs(dx) < TOL_RTSAFE) // Convergence criterion
    {
      return rts;
    }
    f = func(rts, y, z);
    df = dfx(rts, y, z);
    // The one new functon evaluation per iteration
    if (f < 0.0) // Maintain the bracket on the root
    {
      xl = rts;
    } else {
      xh = rts;
    }
  }
  throw("Maximum number of iterations exceeded in rtsafe");
}

/*******************************************************************************************************************
 Using a combination of Newton-Raphson and bisection, return the root of a
 function bracketed between x1 and x2.
 The root will be refined until its accuracy is known within xacc. funcd is a
 user-supplied struct that returns
 the function value as a functor and the first derivative of the function at the
 point x as the function df

 INPUT
 x  : x coordinate of the given point
 y1 : y coordinate of left node for the given point
 y2 : y coordinate of right node for the given point
 z  : z coordinate of the given point

 OUTPUT
 root of function for y-direction
 *******************************************************************************************************************/
double Surface_Cartesian::Rtsafe_y(const Doub_I x, const Doub_I y1,
                                   const Doub_I y2, const Doub_I z) {
  const int MAXIT = 100; // Maximum allowed number of iterations
  double fl, fh, xh, xl, rts, dxold, dx, f, df, temp;
  int j;

  fl = func(x, y1, z);
  fh = func(x, y2, z);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    throw("Root must be bracketed in rtsafe");
  }
  if (fl == 0.0) {
    return y1;
  }
  if (fh == 0.0) {
    return y2;
  }
  if (fl < 0.0) // Orient the search so that f(x1)<0
  {
    xl = y1;
    xh = y2;
  } else {
    xh = y1;
    xl = y2;
  }
  rts = 0.5 * (y1 + y2); // Initialize the guess for root
  dxold = abs(y2 - y1);  // the "stepsize before last"
  dx = dxold;            // and the last step
  f = func(x, rts, z);
  df = dfy(x, rts, z);
  // Bisect if Newton out of range or not decreasing fast enough
  for (j = 0; j < MAXIT; j++) // Loop over allowed iterations
  {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
        (abs(2.0 * f) > abs(dxold * df))) {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts) // Change in root is negligible
      {
        return rts;
      }
    } else // Newton step acceptable, take it
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts) {
        return rts;
      }
    }
    if (abs(dx) < TOL_RTSAFE) // Convergence criterion
    {
      return rts;
    }
    f = func(x, rts, z);
    df = dfy(x, rts, z);
    // The one new functon evaluation per iteration
    if (f < 0.0) // Maintain the bracket on the root
    {
      xl = rts;
    } else {
      xh = rts;
    }
  }
  throw("Maximum number of iterations exceeded in rtsafe");
}

/*******************************************************************************************************************
 Using a combination of Newton-Raphson and bisection, return the root of a
 function bracketed between x1 and x2.
 The root will be refined until its accuracy is known within xacc. funcd is a
 user-supplied struct that returns
 the function value as a functor and the first derivative of the function at the
 point x as the function df

 INPUT
 x  : x coordinate of the given point
 y  : y coordinate of the given point
 z1 : z coordinate of left node for the given point
 z2 : z coordinate of right node for the given point

 OUTPUT
 root of function for z-direction
 *******************************************************************************************************************/
double Surface_Cartesian::Rtsafe_z(const Doub_I x, const Doub_I y,
                                   const Doub_I z1, const Doub_I z2) {
  const int MAXIT = 100; // Maximum allowed number of iterations
  double fl, fh, xh, xl, rts, dxold, dx, f, df, temp;
  int j;

  fl = func(x, y, z1);
  fh = func(x, y, z2);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    throw("Root must be bracketed in rtsafe");
  }
  if (fl == 0.0) {
    return z1;
  }
  if (fh == 0.0) {
    return z2;
  }
  if (fl < 0.0) // Orient the search so that f(x1)<0
  {
    xl = z1;
    xh = z2;
  } else {
    xh = z1;
    xl = z2;
  }
  rts = 0.5 * (z1 + z2); // Initialize the guess for root
  dxold = abs(z2 - z1);  // the "stepsize before last"
  dx = dxold;            // and the last step
  f = func(x, y, rts);
  df = dfz(x, y, rts);
  // Bisect if Newton out of range or not decreasing fast enough
  for (j = 0; j < MAXIT; j++) // Loop over allowed iterations
  {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
        (abs(2.0 * f) > abs(dxold * df))) {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts) // Change in root is negligible
      {
        return rts;
      }
    } else // Newton step acceptable, take it
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts) {
        return rts;
      }
    }
    if (abs(dx) < TOL_RTSAFE) // Convergence criterion
    {
      return rts;
    }
    f = func(x, y, rts);
    df = dfz(x, y, rts);
    // The one new functon evaluation per iteration
    if (f < 0.0) // Maintain the bracket on the root
    {
      xl = rts;
    } else {
      xh = rts;
    }
  }
  throw("Maximum number of iterations exceeded in rtsafe");
}
