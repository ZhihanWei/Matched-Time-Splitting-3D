#ifndef __Interface_H__
#define __Interface_H__

#include <fstream>
#include <string>
#include <vector>
#include "Beta.h"
#include "Constant.h"
#include "Equation.h"
#include "Mesh.h"
#include "New_Data_Type.h"
#include "Surface_Cartesian.h"

using namespace std;

class Intersections {
 private:
  double dx, dy, dz;
  int nx, ny, nz, accuracy;
  double *xi, *yi, *zi, *indicator;

  //***************** Auxilary functions ***********************
	int To1d(Int_I, Int_I, Int_I);
  void Weights(Doub_I, VecDoub_I&, Int_I, Int_I, MatrixDoub_O&);
	void Search_indx(Char_I, Int_I, Int_I, Int_I, VecInt_O&, VecInt_O&);
	void Intersection_display(Intersection_Data, int);
	
	//***************** Initialization functions ***********************
	void Setup_Intersections(Surface_Cartesian&, Beta&, ofstream&);
	void Getdata_irr_z(Int_I, Int_I, Int_I, Surface_Cartesian&,
										 Intersection_Data&, Beta&);
	void Getdata_irr_x(Int_I, Int_I, Int_I, Surface_Cartesian&,
										 Intersection_Data&, Beta&);
	void Getdata_irr_y(Int_I, Int_I, Int_I, Surface_Cartesian&,
										 Intersection_Data&, Beta&);
	void Getdata_cor_z(Int_I, Int_I, Int_I, Surface_Cartesian&,
										 Intersection_Data&, Intersection_Data&, Beta&);
	void Getdata_cor_x(Int_I, Int_I, Int_I, Surface_Cartesian&,
										 Intersection_Data&, Intersection_Data&, Beta&);
	void Getdata_cor_y(Int_I, Int_I, Int_I, Surface_Cartesian&,
										 Intersection_Data&, Intersection_Data&, Beta&);
	
	//***************** MIB functions ***********************
	//----------------- MIB-L2 related ----------------------
	void Setup_MIB_L2(Beta&);
	void Irregular_MIB_L2(Intersection_Data&, Beta&, Doub_I);
	void Irregular_MIB_L2_Recursive(Intersection_Data&, Beta&,
																	Doub_I, Int_I,Int_I, VecDoub_O&);
	void Irregular_MIB_L2_2nd(Intersection_Data&, Doub_I, Int_I);
	void Corner_MIB_L2_2nd(Intersection_Data&, Intersection_Data&, Beta&, Doub_I);
	void Get_irr_weights(Doub_I, Doub_I, Int_I, MatrixDoub_O&, MatrixDoub_O&);
	void Get_cor_weights(Doub_I, Doub_I, Doub_I, MatrixDoub_O&,
											 MatrixDoub_O&, MatrixDoub_O&, MatrixDoub_O&);
	//----------------- MIB-L2 related ----------------------
	void Setup_MIB_L1(Beta&);
	void MIB_L1(Intersection_Data&, Beta&, Doub_I, VecDoub&);
	void Irregular_MIB_L1(Intersection_Data&, Beta&, Doub_I);
	void Corner_MIB_L1(Intersection_Data&, Intersection_Data&, Beta&, Doub_I);
	
	// Coordinator transformation
  void Pmatrix_Setup_x(Intersection_Data&);
  void Pmatrix_Setup_y(Intersection_Data&);
  void Pmatrix_Setup_z(Intersection_Data&);

	
	//***************** Jump functions ***********************
  void Auxiliary_eta_x(Intersection_Data&);
  void Auxiliary_eta_y(Intersection_Data&);
  void Auxiliary_eta_z(Intersection_Data&);
  void Auxiliary_tau_x(Intersection_Data&);
  void Auxiliary_tau_y(Intersection_Data&);
  void Auxiliary_tau_z(Intersection_Data&);

  double Eta_x(Intersection_Data&, CubicDoub&, Beta&);
  double Tau_x(Intersection_Data&, CubicDoub&, Beta&);
  double Eta_y(Intersection_Data&, CubicDoub&, Beta&);
  double Tau_y(Intersection_Data&, CubicDoub&, Beta&);
  double Eta_z(Intersection_Data&, CubicDoub&, Beta&);
  double Tau_z(Intersection_Data&, CubicDoub&, Beta&);

  double Eta_Plus(Intersection_Data&, Beta&, Doub_I);
  double Tau_Plus(Intersection_Data&, Beta&, Doub_I);
  double Eta_Minus(Intersection_Data&, Beta&, Doub_I);
  double Tau_Minus(Intersection_Data&, Beta&, Doub_I);

 public:
  typedef vector<vector<vector<Intersection_Data> > > Cubic_Intersection_Data;

  // Main product of the class with all the informations for an interface
  Cubic_Intersection_Data ifpx, ifpy, ifpz;

  Intersections(Surface_Cartesian&, Mesh&, Beta&, Int_I, Int_I, ofstream&);

	//***************** Auxilary functions ***********************
  void Smallest_gamma_x();

  void Check_size_xy();
  void Check_size_xz();
  void Check_size_yz();

  void Refresh_Fp(Equation&);
  void Error_Fp(ofstream&);
  void Refresh_Jump(Equation&, CubicDoub&, Beta& beta);
  void Error_Jump(ofstream&);
  void Display();
  void Check_coord(Surface_Cartesian&);
};

#endif
