#include "LU.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "Constant.h"

using namespace std;

/*Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU
 decomposition of a
 rowwise permutation of itself. a is input. On output, indx[0..n-1] is an output
 vector
 that records the row permutation effected by the partial pivoting; d is output
 as +/-1
 depending on whether the number of row interchanges was even or odd,
 respectively.
 This routine is used in combination with solve to solve linear equations or
 invert a matrix.*/
LU::LU(MatrixDoub_I &a) : n((int)a.size()), lu(a), aref(a), indx(n) {
  const double TINY = 1.0e-40;
  int i, imax, j, k;
  double big, temp;
  VecDoub vv(n);  // vv stores the imiplicit scaling of each row

  d = 1.0;  // no row interchanges yet
  imax = 0;
 
  // loop over rows to get the implicit scaling information
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++) {
      if ((temp = abs(lu[i][j])) > big) {
        big = temp;
      }
    }
    if (big == 0.0) {
      cout << "Singular matrix in LUdcmp" << endl;
      exit(0);
    };
    // No nonzero largest element
    vv[i] = 1.0 / big;  // save the scaling
  }
  // this is the outermost kth loop
  for (k = 0; k < n; k++) {
    big = 0.0;  // initialize for the search for largest pivot element
    for (i = k; i < n; i++) {
      temp = vv[i] * abs(lu[i][k]);
      if (temp > big) {
        big = temp;
        imax = i;
      }
    }
    // interchange rows
    if (k != imax) {
      for (j = 0; j < n; j++) {
        temp = lu[imax][j];
        lu[imax][j] = lu[k][j];
        lu[k][j] = temp;
      }
      d = -d;            // change the parity of d
      vv[imax] = vv[k];  // interchange the scale factor
    }
    indx[k] = imax;
    if (lu[k][k] == 0.0) {
      lu[k][k] = TINY;  // if the pivot element is zer, the matrix is singular
    }
    for (i = k + 1; i < n; i++) {
      temp = lu[i][k] /= lu[k][k];  // divide by the pivot element
      for (j = k + 1; j < n; j++)  // innermost loop: reduce remaining submatrix
      {
        lu[i][j] -= temp * lu[k][j];
      }
    }
  }
}

/*Solves the set of n linear equations Ax=b using the stored LU decomposition of
 A.
 b[0..n-1] is input as the right-hand side vector b, while x returns the
 solution vector x;
 b and x may reference the same vector, in which case the solution overwrites
 the input.
 This routine takes into account the possibility that b will begin with many
 zero elements,
 so it is efficient for use in matrix inversion.*/
void LU::solve(VecDoub_I &b, VecDoub_O &x) {
  int i, ii, ip, j;
  double sum;

  ii = 0;
  if (b.size() != n || x.size() != n) {
    cout << "LUdcmp::solve bad sizes" << endl;
    exit(0);
  }
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }
  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = x[ip];
    x[ip] = x[i];
    if (ii != 0) {
      for (j = ii - 1; j < i; j++) {
        sum -= lu[i][j] * x[j];
      }
    } else if (sum != 0.0) {  // a nonzero element was encountered, do the sums in the loop above
      ii = i + 1;
    }
    x[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = x[i];
    for (j = i + 1; j < n; j++) {
      sum -= lu[i][j] * x[j];
    }
    x[i] = sum / lu[i][i];  // store a component of the solution vector X
  }
}

/*Solves m sets of n linear equations AX=B using the stored LU decomposition of
 A.
 The matrix b[0..n-1][0..m-1] inputs the right-hand sides, while
 x[0..n-1][0..m-1]
 returns the solution A-1*B. b and x may reference the same matrix, in which
 case
 the solution overwrites the input.*/
void LU::solve(MatrixDoub_I &b, MatrixDoub_O &x) {
  int i, j, m;
  VecDoub xx(n);

  m = (int)b[0].size();
  if ((int)b.size() != n || (int)x.size() != n ||
      (int)b[0].size() != (int)x[0].size()) {
    cout << "LUdcmp::solve bad size" << endl;
    exit(0);
  }
    for (j = 0; j < m; j++) {  // copy and solve each column in turn
      for (i = 0; i < n; i++) {
        xx[i] = b[i][j];
      }
      solve(xx, xx);
      for (i = 0; i < n; i++) {
        x[i][j] = xx[i];
     }
  }
}

/*Using the stored LU decomposition, return in ainv the matrix inverse A^-1*/
void LU::inverse(MatrixDoub_IO &ainv) {
  int i, j;

  ainv.resize(n);
  for (i = 0; i < n; i++) {
    ainv[i].resize(n);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      ainv[i][j] = 0.0;
    }
    ainv[i][i] = 1.0;
  }
  solve(ainv, ainv);
}

/*Using the stored LU decomposition, return the determinant of the matrix A.*/
double LU::det() {
  int i;
  double dd;

  dd = d;
  for (i = 0; i < n; i++) {
    dd *= lu[i][i];
  }
  return dd;
}

/*Improves a solution vector x[0..n-1] of the linear set of equations Ax=b.
 The vectors b[0..n-1] and x[0..n-1] are input. On output, x[0..n-1] is modified,
 to an improved set of values.*/
void LU::mprove(VecDoub_I &b, VecDoub_IO &x) {
  int i, j;
  long double sdp;
  VecDoub r(n);

  for (i = 0; i < n; i++) {
    sdp = -b[i];
    for (j = 0; j < n; j++) {
      sdp += (long double)aref[i][j] * (long double)x[j];
    }
    r[i] = sdp;
  }
  solve(r, r);             // solve the error term
  for (i = 0; i < n; i++)  // subtract the error term from old solution
  {
    x[i] -= r[i];
  }
}
