#ifndef __Constant_H__
#define __Constant_H__

#include <vector>
#include <cmath>

using namespace std;

const double PI = acos(-1);
const int NPRINT = 1;
const double TOL_SETUP = 1.0e-12;
const double TOL_RTSAFE = 1.0e-14;
//region of approximated jump condition ('o': for outside only; others)
const char REG = 'i';
//indicator of jump condition (approximated or analytical; 'r': for real; others)
const char JP = 'r';
const bool DEBUGFLAG = false;
const bool WARNINGS = false;

typedef int Int_I,Int_O,Int_IO;
typedef double Doub_I,Doub_O,Doub_IO;
typedef char Char_I,Char_O,Char_IO;
typedef vector<int> VecInt,VecInt_I,VecInt_O,VecInt_IO;
typedef vector<char> VecChar,VecChar_I,VecChar_O,VecChar_IO;
typedef vector<double> VecDoub,VecDoub_I,VecDoub_O,VecDoub_IO;
typedef vector<vector<double> > MatrixDoub,MatrixDoub_I,MatrixDoub_O,MatrixDoub_IO;
typedef vector<vector<vector<double> > > CubicDoub,CubicDoub_I,CubicDoub_O,CubicDoub_IO;

#endif
