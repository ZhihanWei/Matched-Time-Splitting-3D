#ifndef __New_Data_Type_H__
#define __New_Data_Type_H__

#include <vector>
#include "Constant.h"

using namespace std;

// Three coordinate of interface points
struct Coord {
  double x_value;
  double y_value;
  double z_value;
};
// Indxes for intersects plane, X-direction: indx1 -> iy, indx2 ->iz;
// Y-direction: indx1 -> ix, indx2 -> iz; Z-direction: indx1 -> iy, indx2 -> iz;
struct Line {
  int indx1;
  int indx2;
};
// Jumps on the interface
struct Jump {
  double u;
  double betau_xi;
  double u_eta;
  double u_tau;
  double u_dir;  // Jump on x or y or z direction, u_x or u_y or u_z
  double err;

  double betau_eta;
  double eta_err;
  double betau_tau;
  double tau_err;
};
// Normal direction
struct Local_Coord {
  VecDoub normal;
};
// Outside and inside diffusion coefficients for all used grid nodes ordered
// from small index to large index
struct Diff_Coef {
  double inside;
  double outside;
};
// Zenith and Azimuth of interface point
struct Angle {
  double zenith;
  double azimuth;
};
// Weights for fictitious points
struct Wei {
  MatrixDoub weil;
  MatrixDoub weir;
};
// Error for each interfaces with two FPs
struct Err {
  VecDoub errl;
  VecDoub errr;
};
// All data used to approximate Eta and Tau jumps
struct Approximated_Jumps {
	// on which axis to approximate auxiliary points coordinate location of auxiliary points;
  VecChar auxlaxis;
	VecDoub auxl;           // size: 2; 1st for upper, 2nd for lower
  VecInt uin_auxlnodes;   // indices of upper inside auxiliary nodes
  VecInt uout_auxlnodes;  // indices of upper outside auxiliary nodes
  VecInt lin_auxlnodes;   // indices of lower inside auxiliary nodes
  VecInt lout_auxlnodes;  // indices of lower outside auxiliary nodes
	// approximation region in Omega{+} or Omege{-};
	char region;            // 'o' for Omega{+}, 'i' for Omega{-}
};
struct Intersection_Data {
  Coord coord;
  Line line;
  Wei wei;
  Err err;
  Jump jump;
  Diff_Coef diffcoef;
  Local_Coord local;
  MatrixDoub p;
  Approximated_Jumps eta;
  Approximated_Jumps tau;
  char dir;            // x, y or z direction
  int ID;							 // number of intersection on the grid line 
  int left_loc;        // index of left node of interface intersection
  double gamma;        // distance between left node and interface intersection
};

#endif
