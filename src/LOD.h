#ifndef __LOD_H__
#define __LOD_H__

#include <fstream>
#include <vector>
#include "Beta.h"
#include "Constant.h"
#include "Equation.h"
#include "Intersections.h"
#include "Mesh.h"
#include "New_Data_Type.h"
#include "Surface_Cartesian.h"

using namespace std;

class LOD {
 private:
  int nx, ny, nz;
  double dx, dy, dz, dt;
  double *xi, *yi, *zi, *indicator;

  void D_xx_r_ie(CubicDoub &, CubicDoub &);
  void D_yy_r_ie(CubicDoub &, CubicDoub &);
  void D_zz_r_ie(CubicDoub &, CubicDoub &);
	
	void D_xx_r_cn(Equation &, Intersections &, CubicDoub &, CubicDoub &, Beta &);
  void D_yy_r_cn(Equation &, Intersections &, CubicDoub &, CubicDoub &, Beta &);
  void D_zz_r_cn(Equation &, Intersections &, CubicDoub &, CubicDoub &, Beta &);
	
	void D_xx_l(Int_I, Int_I, Equation &, Intersections &, CubicDoub &,
							Beta &, VecDoub_O &, VecDoub_O &, VecDoub_O &, VecDoub_O &);
  void D_yy_l(Int_I, Int_I, Equation &, Intersections &, CubicDoub &,
							Beta &, VecDoub_O &, VecDoub_O &, VecDoub_O &, VecDoub_O &);
  void D_zz_l(Int_I, Int_I, Equation &, Intersections &, CubicDoub &,
							Beta &, VecDoub_O &, VecDoub_O &, VecDoub_O &, VecDoub_O &);
	
	void Convert2Tri_irr(Int_I, MatrixDoub_I &, VecDoub_O &,
											 VecDoub_O &, VecDoub_O &, VecDoub_O &);
  void Convert2Tri_cor(Int_I, MatrixDoub_I &, VecDoub_O &,
                       VecDoub_O &, VecDoub_O &, VecDoub_O &);
	
  void Operator_weights_x(Beta &, VecDoub_O &, Int_I, Int_I, Int_I, Doub_I);
  void Operator_weights_y(Beta &, VecDoub_O &, Int_I, Int_I, Int_I, Doub_I);
  void Operator_weights_z(Beta &, VecDoub_O &, Int_I, Int_I, Int_I, Doub_I);
	
	void Src_2nd(Equation &, Intersections &, CubicDoub &);
	void Set_bc(Equation &, CubicDoub &);
	void TDMA(VecDoub_I &, VecDoub_I &, VecDoub_I &, VecDoub_I &, VecDoub_O &);
  void Weights(Doub_I, VecDoub_I &, Int_I, Int_I, MatrixDoub_O &);
  int To1d(Int_I, Int_I, Int_I);

 public:
  LOD(Intersections &, Mesh &, Beta &, VecDoub_I);
	
	void Solve_2nd(Equation &, Intersections &, CubicDoub_I &, Beta &);
  void Initialization(Equation &, CubicDoub_I &);
  void Error(Equation &, CubicDoub_I &, ofstream &);
};

#endif
