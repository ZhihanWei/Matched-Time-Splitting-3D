#pragma once

#include <cmath>
#include <vector>

using namespace std;

const double PI = acos(-1);
const int NPRINT = 1;
const double TOL_SETUP = 1.0e-12;
const double TOL_RTSAFE = 1.0e-14;

// region of approximated jump condition ('o': for outside only; others)
const char REG = 'i';
// jump condition ('r': for real jump only; others)
const char JP = 'r';

const bool DEBUGFLAG = false;
const bool WARNINGS = false;

typedef int Int_I, Int_O, Int_IO;
typedef double Doub_I, Doub_O, Doub_IO;
typedef char Char_I, Char_O, Char_IO;
typedef vector<int> VecInt, VecInt_I, VecInt_O, VecInt_IO;
typedef vector<char> VecChar, VecChar_I, VecChar_O, VecChar_IO;
typedef vector<double> VecDoub, VecDoub_I, VecDoub_O, VecDoub_IO;
typedef vector<vector<double>> MatrixDoub, MatrixDoub_I, MatrixDoub_O,
    MatrixDoub_IO;
typedef vector<vector<vector<double>>> CubicDoub, CubicDoub_I, CubicDoub_O,
    CubicDoub_IO;

enum Config {
  MAX_X,
  MIN_X,
  MAX_Y,
  MIN_Y,
  MAX_Z,
  MIN_Z,
  TIME_START,
  TIME_TERMINATE,
  TIME_STEP,
  NX,
  NY,
  NZ,
  SURFACE,
  EQUATION,
  TEMPORAL_METHOD,
  SPATIAL_METHOD,
  SPATIAL_ACCURACY,
  DIFFUSION_COEFFICIENT,
};

enum class Temporal_Method_Type {
  ADI,
  LOD_IE,
  LOD_CN,
  TS,
};

enum class Surface_Type {
  TANGLECUBE,
  CUBE,
  CYLINDER,
  ELLIPSOID,
  CONE,
  PILE,
  TORUS,
  DUPIN_CYCLIDE,
  MOLECULAR,
  HEART,
};

enum class Spatial_Method_Type {
  MIB_V1,
  MIB_V2
};

// Three coordinate of interface points
struct Coord {
  double x_value;
  double y_value;
  double z_value;
};
// Indxes for intersects plane,
// X-direction: indx1 -> iy, indx2 ->iz;
// Y-direction: indx1 -> ix, indx2 -> iz;
// Z-direction: indx1 -> ix, indx2 -> iy;
struct Line {
  int indx1;
  int indx2;
};
// Jumps on the interface
struct Jump {
  double u;
  double betau_xi;  // Normal direction jump
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
  // on which axis to approximate auxiliary points coordinate location of
  // auxiliary points;
  VecChar auxlaxis;
  VecDoub auxl;           // size: 2; 1st for upper, 2nd for lower
  VecInt uin_auxlnodes;   // indices of upper inside auxiliary nodes
  VecInt uout_auxlnodes;  // indices of upper outside auxiliary nodes
  VecInt lin_auxlnodes;   // indices of lower inside auxiliary nodes
  VecInt lout_auxlnodes;  // indices of lower outside auxiliary nodes
  // approximation region in Omega{+} or Omege{-};
  char region;  // 'o' for Omega{+}, 'i' for Omega{-}
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
  char dir;      // x, y or z direction
  int ID;        // number of intersection on the grid line
  int left_loc;  // index of left node of interface intersection
  double gamma;  // distance between left node and interface intersection
};
