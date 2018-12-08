#pragma once

#include <fstream>
#include <vector>

#include "diffusion/beta.h"
#include "constant.h"
#include "equation/equation.h"
#include "spatial/intersection.h"
#include "mesh/mesh.h"
#include "data/data_type.h"
#include "surface/surface_cartesian.h"

using namespace std;

class TS {
private:
  int nx, ny, nz;
  double dx, dy, dz, dt;
  double *xi, *yi, *zi, *indicator;

  void D_xx_r(Equation &, Intersections &, CubicDoub &, CubicDoub &, Beta &);
  void D_yy_r(Equation &, Intersections &, CubicDoub &, CubicDoub &, Beta &);
  void D_zz_r(Equation &, Intersections &, CubicDoub &, CubicDoub &, Beta &);

  void D_xx_l(Int_I, Int_I, Equation &, Intersections &, CubicDoub &, Beta &,
              VecDoub_O &, VecDoub_O &, VecDoub_O &, VecDoub_O &);
  void D_yy_l(Int_I, Int_I, Equation &, Intersections &, CubicDoub &, Beta &,
              VecDoub_O &, VecDoub_O &, VecDoub_O &, VecDoub_O &);
  void D_zz_l(Int_I, Int_I, Equation &, Intersections &, CubicDoub &, Beta &,
              VecDoub_O &, VecDoub_O &, VecDoub_O &, VecDoub_O &);

  void Convert2Tri_irr(Int_I, MatrixDoub_I &, VecDoub_O &, VecDoub_O &,
                       VecDoub_O &, VecDoub_O &);
  void Convert2Tri_cor(Int_I, MatrixDoub_I &, VecDoub_O &, VecDoub_O &,
                       VecDoub_O &, VecDoub_O &);

  void Operator_weights_x(Beta &, VecDoub_O &, Int_I, Int_I, Int_I, Doub_I);
  void Operator_weights_y(Beta &, VecDoub_O &, Int_I, Int_I, Int_I, Doub_I);
  void Operator_weights_z(Beta &, VecDoub_O &, Int_I, Int_I, Int_I, Doub_I);

  void TDMA(VecDoub_I &, VecDoub_I &, VecDoub_I &, VecDoub_I &, VecDoub_O &);
  void Weights(Doub_I, VecDoub_I &, Int_I, Int_I, MatrixDoub_O &);
  int To1d(Int_I, Int_I, Int_I);
  void Src_2nd(Equation &, Intersections &, CubicDoub &);
  void Set_bc(Equation &, CubicDoub &);

public:
  TS(Intersections &, Mesh &, Beta &, VecDoub_I);

  void Solve_2nd(Equation &, Equation &, Equation &, Intersections &,
                 CubicDoub_I &, Beta &);
  void Initialization(Equation &, CubicDoub_I &);
  void Error(Equation &, CubicDoub_I &, ofstream &);
};
