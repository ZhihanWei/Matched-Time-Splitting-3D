#ifndef __Mesh_Initial_H__
#define __Mesh_Initial_H__

#include <string>
#include <vector>
#include "Constant.h"
#include "Surface_Cartesian.h"

using namespace std;

class Mesh
{
private:
    int To1d(int, int, int);
    
public:
    double xl,xr,yl,yr,zl,zr,dx,dy,dz;
    int nx,ny,nz;
    
    double *xi, *yi, *zi, *mesh_value;
    
    Mesh(const VecDoub_I, const VecInt_I, Surface_Cartesian&);
    
    //~Mesh();
    
    void display();
};

#endif
