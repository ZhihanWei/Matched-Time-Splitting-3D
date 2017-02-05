#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "Constant.h"
#include "Data.h"
#include "Surface_Cartesian.h"
#include "Mesh.h"

using namespace std;

/**********************************************************************
                               Constructor
 
 INPUT
 domain : domain of 3D coordinate
 size   : grid size in each direction
 ex     : cartesian surface object
 **********************************************************************/
Mesh::Mesh(const VecDoub_I domain, const VecInt_I size, Surface_Cartesian& ex)
{
    int indx;
    
    xl = domain[0];
    xr = domain[1];
    yl = domain[2];
    yr = domain[3];
    zl = domain[4];
    zr = domain[5];
    nx = size[0];
    ny = size[1];
    nz = size[2]; 
    
    xi = new double[nx];
    yi = new double[ny];
    zi = new double[nz];
    mesh_value = new double[nx*ny*nz];
    
    dx = (xr - xl) / (nx - 1.0);
    dy = (yr - yl) / (ny - 1.0);
    dz = (zr - zl) / (nz - 1.0);
    
    for(int ix = 0; ix < nx; ix++)
    {
        xi[ix] = xl + ix * dx;
    }
    
    for(int iy = 0; iy < ny; iy++)
    {
        yi[iy] = yl + iy * dy;
    }
    
    for(int iz = 0; iz < nz; iz++)
    {
        zi[iz] = zl + iz * dz;
    }
    
    for(int ix = 0; ix < nx; ix++)
    {
        for(int iy = 0; iy < ny; iy++)
        {
            for(int iz = 0; iz < nz; iz++)
            {
                indx = To1d(ix,iy,iz);
                
                if(ex.Set_Node(xi[ix],yi[iy],zi[iz]) > TOL_SETUP)
                {
                    mesh_value[indx] = 1;  //Outside
                }
                else if(ex.Set_Node(xi[ix],yi[iy],zi[iz]) < TOL_SETUP)
                {
                    mesh_value[indx] = -1;  //Inside
                }
                else
                {
                    cout << "Interface on grid, reset the mesh!";
                    exit(0);
                }
            }
        }
    }
}

/*************************
Mesh::~Mesh()
{
    delete []mesh_value;
}
************************/

/**********************************************************
 Convert 3D coefficients to 1D coefficients
 
 INPUT
 ix : index of grid node on x-axis
 iy : index of grid node on y-axis
 iz : index of grid node on z-axis
 
 OUTPUT
 1D coefficients converted from 3D
 **********************************************************/
int Mesh::To1d(int ix, int iy, int iz)
{
    int temp;
    
    temp = (iz*nx*ny) + (iy*nx) + ix;
    
    return temp;
}

/**********************************************************
 Show informations of mesh 
 **********************************************************/
void Mesh::display()
{
    int indx;
    
    cout << "dx = " << dx << "dy = " << dy << "dz = " << dz << endl;
    cout << "nx = " << nx << "ny = " << ny << "nz = " << nz << endl;
    
    for (int ix = 0; ix < nx; ix++)
        for (int iy = 0; iy < ny; iy++)
            for (int iz = 0; iz < nz; iz++)
            {
                indx = To1d(ix,iy,iz);
                
                cout << ix << " " << iy << " " << iz << " " << endl;
                cout << setprecision(5);
                cout << xi[ix] << " " << yi[iy] << " " << zi[iz] << endl;
                cout << mesh_value[indx] << endl;
            }
}