#ifndef __LUDCMP_H__
#define __LUDCMP_H__

#include <vector>
#include "Constant.h"

class LU {
 private:
  int n;
  MatrixDoub lu;       // stores the decomposition
  VecInt indx;         // stores the permutation
  double d;            // used by det
  MatrixDoub_I &aref;  // used only by mprove
 public:
  LU(MatrixDoub_I &a);  // constructor
  void solve(VecDoub_I &b, VecDoub_O &x);        // solve for a single vector
  void solve(MatrixDoub_I &b, MatrixDoub_O &x);  // solve for a matrix
	// calculate matrix inverse A^-1 return determinant of A
	void inverse(MatrixDoub_IO &ainv);
  double det();
	// iterative improvement of a solution
  void mprove(VecDoub_I &b, VecDoub_IO &x);
};

#endif