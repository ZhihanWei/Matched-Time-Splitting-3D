#include "surface/surface_cube.h"
#include "constant.h"
#include <cmath>
#include <iostream>

using namespace std;

/****************************************************
                        Constructor

 INPUT
 arg : arguments of Cube object
 *****************************************************/
Surface_Cube::Surface_Cube(VecDoub_I arg) {
  xl = arg[0];
  xr = arg[1];
  yl = arg[2];
  yr = arg[3];
  zl = arg[4];
  zr = arg[5];
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
double Surface_Cube::Set_Node(Doub_I x, Doub_I y, Doub_I z) const {
  double diff;

  if ((x > xl) && (x < xr) && (y > yl) && (y < yr) && (z > zl) && (z < zr)) {
    diff = -1 + TOL_SETUP;
  } else {
    diff = 1 + TOL_SETUP;
  }

  return diff;
}

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
double Surface_Cube::Gamma_x(Doub_I x1, Doub_I x2, Doub_I y, Doub_I z) {
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

  if ((xl > lw_bd) && (xl < up_bd)) {
    temp = xl;
  } else if ((xr > lw_bd) && (xr < up_bd)) {
    temp = xr;
  } else {
    cout << "No interface is found in x-direction" << endl;
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
double Surface_Cube::Gamma_y(Doub_I x, Doub_I y1, Doub_I y2, Doub_I z) {
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

  if ((yl > lw_bd) && (yl < up_bd)) {
    temp = yl;
  } else if ((yr > lw_bd) && (yr < up_bd)) {
    temp = yr;
  } else {
    cout << "No interface is found in y-direction" << endl;
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
double Surface_Cube::Gamma_z(Doub_I x, Doub_I y, Doub_I z1, Doub_I z2) {
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

  if ((zl > lw_bd) && (zl < up_bd)) {
    temp = zl;
  } else if ((zr > lw_bd) && (zr < up_bd)) {
    temp = zr;
  } else {
    cout << "No interface is found in z-direction" << endl;
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
VecDoub Surface_Cube::normal(Doub_I x, Doub_I y, Doub_I z) {
  VecDoub normal;

  if (x == xl) {
    if ((y == yl) || (y == yr) || (z == zl) || (z == zr)) {
      cout << "Mesh lines are on the interface of cube, reset grid number"
           << endl;
      exit(0);
    }
    normal.push_back(-1);
    normal.push_back(0);
    normal.push_back(0);
  } else if (x == xr) {
    if ((y == yl) || (y == yr) || (z == zl) || (z == zr)) {
      cout << "Mesh lines are on the interface of cube, reset grid number"
           << endl;
      exit(0);
    }
    normal.push_back(1);
    normal.push_back(0);
    normal.push_back(0);
  } else if (y == yl) {
    if ((x == xl) || (x == xr) || (z == zl) || (z == zr)) {
      cout << "Mesh lines are on the interface of cube, reset grid number"
           << endl;
      exit(0);
    }
    normal.push_back(0);
    normal.push_back(-1);
    normal.push_back(0);
  } else if (y == yr) {
    if ((x == xl) || (x == xr) || (z == zl) || (z == zr)) {
      cout << "Mesh lines are on the interface of cube, reset grid number"
           << endl;
      exit(0);
    }
    normal.push_back(0);
    normal.push_back(1);
    normal.push_back(0);
  } else if (z == zl) {
    if ((x == xl) || (x == xr) || (y == yl) || (y == yr)) {
      cout << "Mesh lines are on the interface of cube, reset grid number"
           << endl;
      exit(0);
    }
    normal.push_back(0);
    normal.push_back(0);
    normal.push_back(-1);
  } else if (z == zr) {
    if ((x == xl) || (x == xr) || (y == yl) || (y == yr)) {
      cout << "Mesh lines are on the interface of cube, reset grid number"
           << endl;
      exit(0);
    }
    normal.push_back(0);
    normal.push_back(0);
    normal.push_back(1);
  } else {
    cout << "The point is not an intersection for cube" << endl;
    exit(0);
  }

  return normal;
}

/**********************************************************************
 Show the parameters of cube
 ***********************************************************************/
void Surface_Cube::display() {
  cout << "The surface is a cube with parameters: " << endl;
  cout << "left x coordinate: " << xl << " right x coordinate: " << xr << endl;
  cout << "left y coordinate: " << yl << " right y coordinate: " << yr << endl;
  cout << "left z coordinate: " << zl << " right z coordinate: " << zr << endl;
}

/******************************************************************************
 The following member function should not be defined for cube interface
 *******************************************************************************/
double Surface_Cube::func(const Doub_I x, const Doub_I y,
                          const Doub_I z) const {
  cout << "No implicit function defined for cube surface" << endl;
  exit(0);

  return (0.0);
}

double Surface_Cube::dfx(const Doub_I x, const Doub_I y, const Doub_I z) const {
  cout << "No implicit function defined for cube surface" << endl;
  exit(0);

  return (0.0);
}

double Surface_Cube::dfy(const Doub_I x, const Doub_I y, const Doub_I z) const {
  cout << "No implicit function defined for cube surface" << endl;
  exit(0);

  return (0.0);
}

double Surface_Cube::dfz(const Doub_I x, const Doub_I y, const Doub_I z) const {
  cout << "No implicit function defined for cube surface" << endl;
  exit(0);

  return (0.0);
}

double Surface_Cube::check(Doub_I x, Doub_I y, Doub_I z) {
  cout << "No implicit function defined for cube surface" << endl;
  exit(0);

  return (0.0);
}
