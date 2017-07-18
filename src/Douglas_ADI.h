#ifndef __Douglas_ADI_H__
#define __Douglas_ADI_H__

#include <vector>
#include <fstream>
#include "Constant.h"
#include "New_Data_Type.h"
#include "Equation.h"
#include "Surface_Cartesian.h"
#include "Intersections.h"
#include "Mesh.h"
#include "Beta.h"

using namespace std;

class Douglas_ADI
{
private:
    int nx, ny, nz;
    double dx, dy, dz, dt;
    double *xi, *yi, *zi, *indicator;
    
    double jump_u, jump_ul, jump_ur;
    double jump_betaux, jump_betauxl, jump_betauxr;
    double jump_betauy, jump_betauyl, jump_betauyr;
    double jump_betauz, jump_betauzl, jump_betauzr;
    
    CubicDoub uhs, uhss, uhs2, uh1;
    CubicDoub src;
    
    void Weights(Doub_I, VecDoub_I&, Int_I, Int_I, MatrixDoub_O&);
    void TDMA(VecDoub_I&, VecDoub_I&, VecDoub_I&, VecDoub_I&, VecDoub_O&);
    void PDMA(VecDoub_I&, VecDoub_I&, VecDoub_I&, VecDoub_I&, VecDoub_I&, VecDoub_I&, VecDoub_O&);
    void Diagonal(MatrixDoub_IO&, VecDoub_IO&);
    void Src_2nd(Equation&, Intersections&, Beta&, CubicDoub_I&);
    void Convert2Tri_irr(Int_I, MatrixDoub_I&, VecDoub_O&, VecDoub_O&, VecDoub_O&, VecDoub_O&);
    void Convert2Tri_cor(Int_I, MatrixDoub_I&, VecDoub_O&, VecDoub_O&, VecDoub_O&, VecDoub_O&);
    void Operator_weights_x(Beta&, VecDoub_O&, Int_I, Int_I, Int_I, Doub_I);
    void Operator_weights_y(Beta&, VecDoub_O&, Int_I, Int_I, Int_I, Doub_I);
    void Operator_weights_z(Beta&, VecDoub_O&, Int_I, Int_I, Int_I, Doub_I);
    int To1d(Int_I, Int_I, Int_I);
public:
    Douglas_ADI(Intersections&, Mesh&, Beta&, VecDoub_I, Int_I);
    void Solve_2nd(Equation&, Intersections&, CubicDoub_I&, Beta&);
    void Initialization(Equation&, CubicDoub_I&);
    void Error(Equation&, CubicDoub_I&, ofstream&);
};

#endif
