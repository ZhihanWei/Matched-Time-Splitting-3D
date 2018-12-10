#pragma once

#include <string>
#include <vector>
#include "constant.h"
#include "surface/surface_cartesian.h"

using namespace std;

class Mesh {
private:
  int To1d(int, int, int);

public:
  double xl, xr, yl, yr, zl, zr, dx, dy, dz;
  int nx, ny, nz;

  double *xi, *yi, *zi, *mesh_value;

  Mesh(const VecDoub_I, const VecInt_I, Surface_Cartesian &);

  //~Mesh();

  void display();
};
