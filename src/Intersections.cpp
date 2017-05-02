#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "Intersections.h"
#include "LU.h"

using namespace std;

/**************************************************************************************************************
 Constructor
 
 INPUT
 ex          : cartesian surface object
 mesh        : mesh object (at least contains: nx, dx, grids location xi[n] and grid indicator mesh_value[i])
 beta        : diffusion coefficient object representing beta^{-} and beta^{+}
 TOL_ITYPE   : tolerance of itype
 accuracy    : accuracy of scheme order
 **************************************************************************************************************/
Intersections::Intersections(Surface_Cartesian& ex, Mesh& mesh, Beta& beta, Doub_I TOL_ITYPE, Int_I in_accuracy, ofstream& out_file)
{
    tol_type = TOL_ITYPE;
    
    nx = mesh.nx;
    ny = mesh.ny;
    nz = mesh.nz;
    dx = mesh.dx;
    dy = mesh.dy;
    dz = mesh.dz;
    
    xi = &mesh.xi[0];
    yi = &mesh.yi[0];
    zi = &mesh.zi[0];
    indicator = &mesh.mesh_value[0];
    
    accuracy = in_accuracy;
    
    Setup_Intersections(ex,beta,out_file);
    
    Setup_MIB(beta);
}

/**********************************************************
 Convert 3D coefficients to 1D coefficients
 
 INPUT
 ix : index of grid node on x-axis
 iy : index of grid node on y-axis
 iz : index of grid node on z-axis
 
 OUTPUT
 1D coefficients converted from 3D
 **********************************************************/
int Intersections::To1d(Int_I ix, Int_I iy, Int_I iz)
{
    int temp;
    
    temp = (iz*nx*ny) + (iy*nx) + ix;
    
    return temp;
}

/**********************************************************
 Setup initial values for the interseciton node
 
 INPUT
 ex : cartesian surface object
 **********************************************************/
void Intersections::Setup_Intersections(Surface_Cartesian& ex, Beta& beta, ofstream& out_file)
{
    int ifp,indx,indx1,indx2,indx3,indx4;
    Intersection_Data inter_node,inter_node_left,inter_node_right;
    
    out_file << "Z-direction computation: " << endl;
    //For the line of intersect on X-plane and Y-plane, check in Z-direction
    ifpz.resize(nx);
    for(int ix = 1; ix < nx-1; ix++)
    {
        ifpz[ix].resize(ny);
        for(int iy = 1; iy < ny-1; iy++)
        {
            ifp = 0;
            
            int iz = 1;
            while(iz < nz-1)
            {
                indx = To1d(ix,iy,iz);
                indx1 = To1d(ix,iy,iz+1);
                indx2 = To1d(ix,iy,iz+2);
                indx3 = To1d(ix,iy,iz+3);
                indx4 = To1d(ix,iy,iz+4);
                
                if(indicator[indx] * indicator[indx1] < 0)
                {
                    //Irregular interface point
                    if(indicator[indx1] * indicator[indx2] > 0)
                    {
                        if(accuracy == 4)
                        {
                            if(indicator[indx1] * indicator[indx4] < 0)
                            {
                                cout << "4th order cannot deal with corner interface" << endl;
                                exit(0);
                            }
                        }
                        
                        ifp += 1;
                        inter_node.ID = ifp;
                        Getdata_irr_z(ix,iy,iz,ex,inter_node,beta);
                        
                        ifpz[ix][iy].push_back(inter_node);
                        
                        iz += 2;
                    }
                    //Corner interface point
                    else
                    {
                        if(indicator[indx2] * indicator[indx3] < 0)
                        {
                            cout << "Two corner interfaces coherent in Z-direction at " << "IX: " << ix << " IY: " << iy << endl;
                            exit(0);
                        }
                        else
                        {
                            out_file << "Corner interface is found in Z-direction at " << "IX: " << ix << " IY: " << iy << endl;
                            
                            ifp += 1;
                            inter_node_left.ID = -ifp;
                            ifp += 1;
                            inter_node_right.ID = -ifp;
                            Getdata_cor_z(ix,iy,iz,ex,inter_node_left,inter_node_right,beta);
                            
                            ifpz[ix][iy].push_back(inter_node_left);
                            ifpz[ix][iy].push_back(inter_node_right);
                            
                            iz += 3;
                        }
                    }
                }
                else
                {
                    iz += 1;
                }
            }
            if(ifp%2 == 1)
            {
                cout << "It's an open surface in Z-direction, change the size of domain" << endl;
                exit(0);
            }
        }
    }
    
    out_file << endl << "Y-direction computation: " << endl;
    //For the line of intersect on X-plane and Z-plane, check in Y-direction
    ifpy.resize(nx);
    for(int ix = 1; ix < nx-1; ix++)
    {
        ifpy[ix].resize(nz);
        for(int iz = 1; iz < nz-1; iz++)
        {
            ifp = 0;
            
            int iy = 1;
            while(iy < ny-1)
            {
                indx = To1d(ix,iy,iz);
                indx1 = To1d(ix,iy+1,iz);
                indx2 = To1d(ix,iy+2,iz);
                indx3 = To1d(ix,iy+3,iz);
                indx4 = To1d(ix,iy+4,iz);
                
                if(indicator[indx] * indicator[indx1] < 0)
                {
                    //Irregular interface point
                    if(indicator[indx1] * indicator[indx2] > 0)
                    {
                        if(accuracy == 4)
                        {
                            if(indicator[indx1] * indicator[indx4] < 0)
                            {
                                cout << "4th order cannot deal with corner interface" << endl;
                                exit(0);
                            }
                        }
                        ifp += 1;
                        inter_node.ID = ifp;
                        Getdata_irr_y(ix,iy,iz,ex,inter_node,beta);
                        
                        ifpy[ix][iz].push_back(inter_node);
                        
                        iy += 2;
                    }
                    //Corner interface point
                    else
                    {
                        if(indicator[indx2] * indicator[indx3] < 0)
                        {
                            cout << "Two corner interfaces coherent in Y-direction at " << "IX: " << ix << " IZ: " << iz << endl;
                            exit(0);
                        }
                        else
                        {
                            out_file << "Corner interface is found in Y-direction at " << "IX: " << ix << " IZ: " << iz << endl;
                            
                            ifp += 1;
                            inter_node_left.ID = -ifp;
                            ifp += 1;
                            inter_node_right.ID = -ifp;
                            
                            Getdata_cor_y(ix,iy,iz,ex,inter_node_left,inter_node_right,beta);
                            
                            ifpy[ix][iz].push_back(inter_node_left);
                            ifpy[ix][iz].push_back(inter_node_right);
                            
                            iy += 3;
                        }
                    }
                }
                else
                {
                    iy += 1;
                }
            }
            if(ifp%2 == 1)
            {
                cout << "It's an open surface in Y-direction, change the size of domain" << endl;
                exit(0);
            }
        }
    }
    
    out_file << endl << "X-direction computation: " << endl;
    //For the line of intersect on Y-plane and Z-plane, check in X-direction
    ifpx.resize(ny);
    for(int iy = 1; iy < ny-1; iy++)
    {
        ifpx[iy].resize(nz);
        for(int iz = 1; iz < nz-1; iz++)
        {
            ifp = 0;
            
            int ix = 1;
            while(ix < nx-1)
            {
                indx = To1d(ix,iy,iz);
                indx1 = To1d(ix+1,iy,iz);
                indx2 = To1d(ix+2,iy,iz);
                indx3 = To1d(ix+3,iy,iz);
                indx4 = To1d(ix+4,iy,iz);
                
                if(indicator[indx] * indicator[indx1] < 0)
                {
                    //Irregular interface point
                    if(indicator[indx1] * indicator[indx2] > 0)
                    {
                        if(accuracy == 4)
                        {
                            if(indicator[indx1] * indicator[indx4] < 0)
                            {
                                cout << "4th order cannot deal with corner interface" << endl;
                                exit(0);
                            }
                        }
                        ifp += 1;
                        inter_node.ID = ifp;
                        Getdata_irr_x(ix,iy,iz,ex,inter_node,beta);
                        
                        ifpx[iy][iz].push_back(inter_node);
                        
                        ix += 2;
                    }
                    //Corner interface point
                    else
                    {
                        if(indicator[indx2] * indicator[indx3] < 0)
                        {
                            cout << "Two corner interfaces coherent in X-direction at " << "IY: " << iy << " IZ: " << iz << endl;
                            exit(0);
                        }
                        else
                        {
                            out_file << "Corner interface is found in X-direction at " << "IY: " << iy << " IZ: " << iz << endl;
                            
                            ifp += 1;
                            inter_node_left.ID = -ifp;
                            ifp += 1;
                            inter_node_right.ID = -ifp;
                            Getdata_cor_x(ix,iy,iz,ex,inter_node_left,inter_node_right,beta);
                            
                            ifpx[iy][iz].push_back(inter_node_left);
                            ifpx[iy][iz].push_back(inter_node_right);
                            
                            ix += 3;
                        }
                    }
                }
                else
                {
                    ix += 1;
                }
            }
            if(ifp%2 == 1)
            {
                cout << "It's an open surface in X-direction, change the size of domain" << endl;
                exit(0);
            }
        }
    }
    
    out_file << endl;
}

/*************************************************************************************************
 Initialize irregular interface node on xy plane in z-direction without MIB weights and error
 
 INPUT
 ix         : index of grid node on x-axis
 iy         : index of grid node on y-axis
 iz         : index of closest left grid node on z-axis
 ex         : cartesian surface object
 inter_node : intersection node need to be initialized
 *************************************************************************************************/
void Intersections::Getdata_irr_z(Int_I ix, Int_I iy, Int_I iz, Surface_Cartesian& ex, Intersection_Data& inter_node, Beta& beta)
{
    inter_node.dir = 'z';
    
    inter_node.err.errr.resize(accuracy/2);
    inter_node.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node.err.errr[i] = 0;
        inter_node.err.errl[i] = 0;
    }
    
    inter_node.coord.x_value = xi[ix];
    inter_node.coord.y_value = yi[iy];
    inter_node.coord.z_value = ex.Gamma_z(xi[ix],yi[iy],zi[iz],zi[iz+1]);
    
    inter_node.line.indx1 = ix;
    inter_node.line.indx2 = iy;
    
    inter_node.left_loc = iz;
    
    inter_node.gamma = inter_node.coord.z_value - zi[iz];
    
    inter_node.jump.u = 0;
    inter_node.jump.betau_xi = 0;
    inter_node.jump.u_eta = 0;
    inter_node.jump.u_tau = 0;
    inter_node.jump.u_dir = 0;
    inter_node.jump.err = 0;
    
    inter_node.diffcoef.inside = beta.Inside(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    inter_node.diffcoef.outside = beta.Outside(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    
    inter_node.local.normal = ex.normal(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    
    Pmatrix_Setup_z(inter_node);
    Auxiliary_eta_z(inter_node);
    Auxiliary_tau_z(inter_node);
}

/*************************************************************************************************
 Initialize irregular interface node on xz plane in y-direction without MIB weights and error
 
 INPUT
 ix         : index of grid node on x-axis
 iy         : index of closest left grid node on y-axis
 iz         : index of grid node on z-axis
 ex         : cartesian surface object
 inter_node : intersection node need to be initialized
 *************************************************************************************************/
void Intersections::Getdata_irr_y(Int_I ix, Int_I iy, Int_I iz, Surface_Cartesian& ex, Intersection_Data& inter_node, Beta& beta)
{
    inter_node.dir = 'y';
    
    inter_node.err.errr.resize(accuracy/2);
    inter_node.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node.err.errr[i] = 0;
        inter_node.err.errl[i] = 0;
    }
    
    inter_node.coord.x_value = xi[ix];
    inter_node.coord.y_value = ex.Gamma_y(xi[ix],yi[iy],yi[iy+1],zi[iz]);
    inter_node.coord.z_value = zi[iz];
    
    inter_node.line.indx1 = ix;
    inter_node.line.indx2 = iz;
    
    inter_node.left_loc = iy;
    
    inter_node.gamma = inter_node.coord.y_value - yi[iy];
    
    inter_node.jump.u = 0;
    inter_node.jump.betau_xi = 0;
    inter_node.jump.u_eta = 0;
    inter_node.jump.u_tau = 0;
    inter_node.jump.u_dir = 0;
    inter_node.jump.err = 0;
    
    inter_node.diffcoef.inside = beta.Inside(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    inter_node.diffcoef.outside = beta.Outside(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    
    inter_node.local.normal = ex.normal(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    
    Pmatrix_Setup_y(inter_node);
    Auxiliary_eta_y(inter_node);
    Auxiliary_tau_y(inter_node);
}

/*************************************************************************************************
 Initialize irregular interface node on yz plane in x-direction without MIB weights and error
 
 INPUT
 ix         : index of closest left grid node on x-axis
 iy         : index of grid node on y-axis
 iz         : index of grid node on z-axis
 ex         : cartesian surface object
 inter_node : intersection node need to be initialized
 *************************************************************************************************/
void Intersections::Getdata_irr_x(Int_I ix, Int_I iy, Int_I iz, Surface_Cartesian& ex, Intersection_Data& inter_node, Beta& beta)
{
    inter_node.dir = 'x';
    
    inter_node.err.errr.resize(accuracy/2);
    inter_node.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node.err.errr[i] = 0;
        inter_node.err.errl[i] = 0;
    }
    
    inter_node.coord.x_value = ex.Gamma_x(xi[ix],xi[ix+1],yi[iy],zi[iz]);
    inter_node.coord.y_value = yi[iy];
    inter_node.coord.z_value = zi[iz];
    
    inter_node.line.indx1 = iy;
    inter_node.line.indx2 = iz;
    
    inter_node.left_loc = ix;
    
    inter_node.gamma = inter_node.coord.x_value - xi[ix];
    
    inter_node.jump.u = 0;
    inter_node.jump.betau_xi = 0;
    inter_node.jump.u_eta = 0;
    inter_node.jump.u_tau = 0;
    inter_node.jump.u_dir = 0;
    inter_node.jump.err = 0;
    
    inter_node.diffcoef.inside = beta.Inside(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    inter_node.diffcoef.outside = beta.Outside(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    
    inter_node.local.normal = ex.normal(inter_node.coord.x_value, inter_node.coord.y_value, inter_node.coord.z_value);
    
    Pmatrix_Setup_x(inter_node);
    Auxiliary_eta_x(inter_node);
    Auxiliary_tau_x(inter_node);
}

/*************************************************************************************************
 Initialize corner interface nodes on xy plane in z-direction without MIB weights and error
 
 INPUT
 ix               : index of grid node on x-axis
 iy               : index of grid node on y-axis
 iz               : index of closest left grid node on z-axis
 ex               : cartesian surface object
 inter_node_left  : intersection node need to be initialized
 inter_node_right : right intersection node need to be initialized
 *************************************************************************************************/
void Intersections::Getdata_cor_z(Int_I ix, Int_I iy, Int_I iz, Surface_Cartesian& ex, Intersection_Data& inter_node_left, Intersection_Data& inter_node_right, Beta& beta)
{
    //Left corner interface
    inter_node_left.dir = 'z';
    
    inter_node_left.err.errr.resize(accuracy/2);
    inter_node_left.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node_left.err.errr[i] = 0;
        inter_node_left.err.errl[i] = 0;
    }
    
    inter_node_left.coord.x_value = xi[ix];
    inter_node_left.coord.y_value = yi[iy];
    inter_node_left.coord.z_value = ex.Gamma_z(xi[ix],yi[iy],zi[iz],zi[iz+1]);
    
    inter_node_left.line.indx1 = ix;
    inter_node_left.line.indx2 = iy;
    
    inter_node_left.left_loc = iz;
    
    inter_node_left.gamma = inter_node_left.coord.z_value - zi[iz];
    
    inter_node_left.jump.u = 0;
    inter_node_left.jump.betau_xi = 0;
    inter_node_left.jump.u_eta = 0;
    inter_node_left.jump.u_tau = 0;
    inter_node_left.jump.u_dir = 0;
    inter_node_left.jump.err = 0;
    
    inter_node_left.diffcoef.inside = beta.Inside(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    inter_node_left.diffcoef.outside = beta.Outside(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    
    inter_node_left.local.normal = ex.normal(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    
    Pmatrix_Setup_z(inter_node_left);
    Auxiliary_eta_z(inter_node_left);
    Auxiliary_tau_z(inter_node_left);
    
    //Right corner interface
    inter_node_right.dir = 'z';
    
    inter_node_right.err.errr.resize(accuracy/2);
    inter_node_right.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node_right.err.errr[i] = 0;
        inter_node_right.err.errl[i] = 0;
    }
    
    inter_node_right.coord.x_value = xi[ix];
    inter_node_right.coord.y_value = yi[iy];
    inter_node_right.coord.z_value = ex.Gamma_z(xi[ix],yi[iy],zi[iz+1],zi[iz+2]);
    
    inter_node_right.line.indx1 = ix;
    inter_node_right.line.indx2 = iy;
    
    inter_node_right.left_loc = iz+1;
    
    inter_node_right.gamma = inter_node_right.coord.z_value - zi[iz+1];
    
    inter_node_right.jump.u = 0;
    inter_node_right.jump.betau_xi = 0;
    inter_node_right.jump.u_eta = 0;
    inter_node_right.jump.u_tau = 0;
    inter_node_right.jump.u_dir = 0;
    inter_node_right.jump.err = 0;
    
    inter_node_right.diffcoef.inside = beta.Inside(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    inter_node_right.diffcoef.outside = beta.Outside(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    
    inter_node_right.local.normal = ex.normal(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    
    Pmatrix_Setup_z(inter_node_right);
    Auxiliary_eta_z(inter_node_right);
    Auxiliary_tau_z(inter_node_right);
}

/*************************************************************************************************
 Initialize corner interface nodes on xz plane in y-direction without MIB weights and error
 
 INPUT
 ix               : index of grid node on x-axis
 iy               : index of closest left grid node on y-axis
 iz               : index of grid node on z-axis
 ex               : cartesian surface object
 inter_node_left  : intersection node need to be initialized
 inter_node_right : right intersection node need to be initialized
 *************************************************************************************************/
void Intersections::Getdata_cor_y(Int_I ix, Int_I iy, Int_I iz, Surface_Cartesian& ex, Intersection_Data& inter_node_left, Intersection_Data& inter_node_right, Beta& beta)
{
    //Left corner interface
    inter_node_left.dir = 'y';
    
    inter_node_left.err.errr.resize(accuracy/2);
    inter_node_left.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node_left.err.errr[i] = 0;
        inter_node_left.err.errl[i] = 0;
    }
    
    inter_node_left.coord.x_value = xi[ix];
    inter_node_left.coord.y_value = ex.Gamma_y(xi[ix],yi[iy],yi[iy+1],zi[iz]);
    inter_node_left.coord.z_value = zi[iz];
    
    inter_node_left.line.indx1 = ix;
    inter_node_left.line.indx2 = iz;
    
    inter_node_left.left_loc = iy;
    
    inter_node_left.gamma = inter_node_left.coord.y_value - yi[iy];
    
    inter_node_left.jump.u = 0;
    inter_node_left.jump.betau_xi = 0;
    inter_node_left.jump.u_eta = 0;
    inter_node_left.jump.u_tau = 0;
    inter_node_left.jump.u_dir = 0;
    inter_node_left.jump.err = 0;
    
    inter_node_left.diffcoef.inside = beta.Inside(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    inter_node_left.diffcoef.outside = beta.Outside(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    
    inter_node_left.local.normal = ex.normal(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    
    Pmatrix_Setup_y(inter_node_left);
    Auxiliary_eta_y(inter_node_left);
    Auxiliary_tau_y(inter_node_left);
    
    //Right corner interface
    inter_node_right.dir = 'y';
    
    inter_node_right.err.errr.resize(accuracy/2);
    inter_node_right.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node_right.err.errr[i] = 0;
        inter_node_right.err.errl[i] = 0;
    }
    
    inter_node_right.coord.x_value = xi[ix];
    inter_node_right.coord.y_value = ex.Gamma_y(xi[ix],yi[iy+1],yi[iy+2],zi[iz]);
    inter_node_right.coord.z_value = zi[iz];
    
    inter_node_right.line.indx1 = ix;
    inter_node_right.line.indx2 = iz;
    
    inter_node_right.left_loc = iy+1;
    
    inter_node_right.gamma = inter_node_right.coord.y_value - yi[iy+1];
    
    inter_node_right.jump.u = 0;
    inter_node_right.jump.betau_xi = 0;
    inter_node_right.jump.u_eta = 0;
    inter_node_right.jump.u_tau = 0;
    inter_node_right.jump.u_dir = 0;
    inter_node_right.jump.err = 0;
    
    inter_node_right.diffcoef.inside = beta.Inside(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    inter_node_right.diffcoef.outside = beta.Outside(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    
    inter_node_right.local.normal = ex.normal(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    
    Pmatrix_Setup_y(inter_node_right);
    Auxiliary_eta_y(inter_node_right);
    Auxiliary_tau_y(inter_node_right);
}

/*************************************************************************************************
 Initialize corner interface nodes on yz plane in x-direction without MIB weights and error
 
 INPUT
 ix               : index of closest left grid node on x-axis
 iy               : index of grid node on y-axis
 iz               : index of grid node on z-axis
 ex               : cartesian surface object
 inter_node_left  : intersection node need to be initialized
 inter_node_right : right intersection node need to be initialized
 *************************************************************************************************/
void Intersections::Getdata_cor_x(Int_I ix, Int_I iy, Int_I iz, Surface_Cartesian& ex, Intersection_Data& inter_node_left, Intersection_Data& inter_node_right, Beta& beta)
{
    //Left corner interface
    inter_node_left.dir = 'x';
    
    inter_node_left.err.errr.resize(accuracy/2);
    inter_node_left.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node_left.err.errr[i] = 0;
        inter_node_left.err.errl[i] = 0;
    }
    
    inter_node_left.coord.x_value = ex.Gamma_x(xi[ix],xi[ix+1],yi[iy],zi[iz]);
    inter_node_left.coord.y_value = yi[iy];
    inter_node_left.coord.z_value = zi[iz];
    
    inter_node_left.line.indx1 = iy;
    inter_node_left.line.indx2 = iz;
    
    inter_node_left.left_loc = ix;
    
    inter_node_left.gamma = inter_node_left.coord.x_value - xi[ix];
    
    inter_node_left.jump.u = 0;
    inter_node_left.jump.betau_xi = 0;
    inter_node_left.jump.u_eta = 0;
    inter_node_left.jump.u_tau = 0;
    inter_node_left.jump.u_dir = 0;
    inter_node_left.jump.err = 0;
    
    inter_node_left.diffcoef.inside = beta.Inside(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    inter_node_left.diffcoef.outside = beta.Outside(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    
    inter_node_left.local.normal = ex.normal(inter_node_left.coord.x_value, inter_node_left.coord.y_value, inter_node_left.coord.z_value);
    
    Pmatrix_Setup_x(inter_node_left);
    Auxiliary_eta_x(inter_node_left);
    Auxiliary_tau_x(inter_node_left);
    
    //Right corner interface
    inter_node_right.dir = 'x';
    
    inter_node_right.err.errr.resize(accuracy/2);
    inter_node_right.err.errl.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node_right.err.errr[i] = 0;
        inter_node_right.err.errl[i] = 0;
    }
    
    inter_node_right.coord.x_value = ex.Gamma_x(xi[ix+1],xi[ix+2],yi[iy],zi[iz]);
    inter_node_right.coord.y_value = yi[iy];
    inter_node_right.coord.z_value = zi[iz];
    
    inter_node_right.line.indx1 = iy;
    inter_node_right.line.indx2 = iz;
    
    inter_node_right.left_loc = ix+1;
    
    inter_node_right.gamma = inter_node_right.coord.x_value - xi[ix+1];
    
    inter_node_right.jump.u = 0;
    inter_node_right.jump.betau_xi = 0;
    inter_node_right.jump.u_eta = 0;
    inter_node_right.jump.u_tau = 0;
    inter_node_right.jump.u_dir = 0;
    inter_node_right.jump.err = 0;
    
    inter_node_right.diffcoef.inside = beta.Inside(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    inter_node_right.diffcoef.outside = beta.Outside(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    
    inter_node_right.local.normal = ex.normal(inter_node_right.coord.x_value, inter_node_right.coord.y_value, inter_node_right.coord.z_value);
    
    Pmatrix_Setup_x(inter_node_right);
    Auxiliary_eta_x(inter_node_right);
    Auxiliary_tau_x(inter_node_right);
}

/*************************************************************************************************
 Initialize intersection nodes(both irregular and corner) calculated fictitious points' weights
 *************************************************************************************************/
void Intersections::Setup_MIB(Beta& beta)
{
    Intersection_Data inter_node;
    int ip;
    
    //Z-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iy = 1; iy < ny-1; iy++)
        {
            ip = 0;
            while(ip < ifpz[ix][iy].size())
            {
                if(ifpz[ix][iy][ip].ID > 0)
                {
                    Irregular_MIB(ifpz[ix][iy][ip],beta,dz);
                    ip += 1;
                }
                else if(ifpz[ix][iy][ip].ID < 0)
                {
                    Corner_MIB_2nd(ifpz[ix][iy][ip],ifpz[ix][iy][ip+1],beta,dz);
                    ip += 2;
                }
                else
                {
                    cout << "Wrong value for interface node ID in Z-direction" << endl;
                    exit(0);
                }
            }
        }
    }
    
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            ip = 0;
            while(ip < ifpy[ix][iz].size())
            {
                if(ifpy[ix][iz][ip].ID > 0)
                {
                    Irregular_MIB(ifpy[ix][iz][ip],beta,dy);
                    ip += 1;
                }
                else if(ifpy[ix][iz][ip].ID < 0)
                {
                    Corner_MIB_2nd(ifpy[ix][iz][ip],ifpy[ix][iz][ip+1],beta,dy);
                    ip += 2;
                }
                else
                {
                    cout << "Wrong value for interface node ID in Y-direction" << endl;
                    exit(0);
                }
            }
        }
    }
    
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            ip = 0;
            while(ip < ifpx[iy][iz].size())
            {
                if(ifpx[iy][iz][ip].ID > 0)
                {
                    Irregular_MIB(ifpx[iy][iz][ip],beta,dx);
                    ip += 1;
                }
                else if(ifpx[iy][iz][ip].ID < 0)
                {
                    Corner_MIB_2nd(ifpx[iy][iz][ip],ifpx[iy][iz][ip+1],beta,dx);
                    ip += 2;
                }
                else
                {
                    cout << "Wrong value for interface node ID in X-direction" << endl;
                    exit(0);
                }
            }
        }
    }
}

/***************************************************************************
 Recursive MIB weights calculation function
 
 INPUT
 inter_node : intersection node struct
 dv         : grid mesh size at given direction
 ***************************************************************************/
void Intersections::Irregular_MIB(Intersection_Data& inter_node, Beta& beta, Doub_I dv)
{
    VecDoub x;
    int onefp_unknowns,stlength,fpno;
    
    onefp_unknowns = accuracy*2+2;                  //numbers of unknowns for one fictitious points
    stlength = accuracy+1;                          //length of FD stencil
    fpno = accuracy/2;                              //numbers of fictitious points on oneside at current intersection given accuracy
    
    x.resize(2*onefp_unknowns);
    
    inter_node.wei.weil.resize(accuracy/2);
    inter_node.wei.weir.resize(accuracy/2);
    for(int i = 0; i < accuracy/2; i++)
    {
        inter_node.wei.weil[i].resize(onefp_unknowns);
        inter_node.wei.weir[i].resize(onefp_unknowns);
    }
    
    for(int i = 0; i < fpno; i++)
    {
        for(int j = 0; j < 2*onefp_unknowns; j++)
        {
            x[j] = 0;
        }
        
        Irregular_MIB_Recursive(inter_node,beta,dv,onefp_unknowns,stlength,x);
        
        for(int j = 0; j < onefp_unknowns; j++)
        {
            inter_node.wei.weil[i][j] = x[j];
            inter_node.wei.weir[i][j] = x[onefp_unknowns+j];
        }
        
        stlength += 1;
    }
}

/***************************************************************************************
 Core MIB algorithm, solve fictitious points by pairs
 
 INPUT
 inter_node     : intersection node struct
 dv             : grid mesh size at given direction
 onefp_unknowns : number of nodes in oneside of central finite difference
 stlength       : length of central finite difference scheme
 
 OUTPUT
 x : vector of weights for a pair of fictitious points
 ****************************************************************************************/
void Intersections::Irregular_MIB_Recursive(Intersection_Data& inter_node, Beta& beta, Doub_I dv, Int_I onefp_unknowns, Int_I stlength, VecDoub_O& x)
{
    MatrixDoub A,weil,weir;
    VecDoub B;
    int order;
    int dev_order,total_unknowns,oneside_pts;
    double coefl,coefr;                               //diffusion coefficients along interface
    
    total_unknowns = 2*onefp_unknowns;                //total unknowns need to be solved
    oneside_pts = accuracy;                           //numbers of real points used for one fictitous point
    dev_order = 1;                                    //highest order of derivative in two equations
    
    B.resize(total_unknowns);
    A.resize(total_unknowns);
    for(int i = 0; i < total_unknowns; i++)
    {
        A[i].resize(total_unknowns);
    }
    for(int i = 0; i < total_unknowns; i++)
    {
        B[i] = 0;
        for(int j = 0; j < total_unknowns; j++)
        {
            A[i][j] = 0;
        }
    }
    
    weil.resize(dev_order+1);
    weir.resize(dev_order+1);
    for(int i = 0; i < dev_order+1; i++)
    {
        weil[i].resize(stlength);
        weir[i].resize(stlength);
    }
    
    Get_irr_weights(inter_node.gamma,dv,oneside_pts,weil,weir);
    
    //Zero order: U^{-} = U^{+} - [U]
    order = 0;
    
    //All real points weight
    for(int i = 0; i < oneside_pts; i++)
    {
        B[i] += -weir[order][i];                                               //Right FP
        B[oneside_pts+i] += weil[order][stlength-oneside_pts+i];               //Left FP
    }
    //All calculated fictitious points weight
    for(int i = 0; i < stlength-oneside_pts-1; i++)
    {
        for(int j = 0; j < onefp_unknowns; j++)
        {
            B[j] += -weir[order][oneside_pts+i]*inter_node.wei.weir[i][j];
            B[j] += weil[order][stlength-oneside_pts-1-i]*inter_node.wei.weil[i][j];
        }
    }
    //Two fictitious points weight
    for(int i = 0; i < onefp_unknowns; i++)
    {
        A[i][i] = -weil[order][0];                                              //Left FP
        A[i][onefp_unknowns+i] = weir[order][stlength-1];                       //Right FP
    }
    
    if(abs(inter_node.ID)%2 == 0)                   //"-" => "+", U^{-} = U^{+} - [U], F = G - [U]
    {
        B[onefp_unknowns-2] += -1;
    }
    else                                            //"+" => "-", U^{-} = U^{+} - [U], G + [U] = F
    {
        B[onefp_unknowns-2] += 1;
    }
    
    //First order: BETA^{-}*U_{x} ^{-} = BETA^{+}*U_{x}^{+} - [BETA U_{x}]
    order = 1;
    
    if(abs(inter_node.ID)%2 == 0)                   //"-" => "+", BETA^{-}U_{X}^{-} = BETA^{+}U_{X}^{+} - [BETA U_{X}], F = G - [BETA U_{X}]
    {
        coefr = inter_node.diffcoef.inside;
        coefl = inter_node.diffcoef.outside;
    }
    else                                           //"+" => "-", BETA^{-}U_{X}^{-} = BETA^{+}U_{X}^{+} - [BETA U_{X}], G + [BETA U_{X}] = F
    {
        coefr = inter_node.diffcoef.outside;
        coefl = inter_node.diffcoef.inside;
    }
    //All real points weight
    for(int i = 0; i < oneside_pts; i++)
    {
        B[onefp_unknowns+i] += -weir[order][i]*coefr;                                               //Right FP
        B[onefp_unknowns+oneside_pts+i] += weil[order][stlength-oneside_pts+i]*coefl;               //Left FP
    }
    //All calculated fictitious points weight
    for(int i = 0; i < stlength-oneside_pts-1; i++)
    {
        for(int j = 0; j < onefp_unknowns; j++)
        {
            B[onefp_unknowns+j] += -weir[order][oneside_pts+i]*coefr*inter_node.wei.weir[i][j];                     //Right FP
            B[onefp_unknowns+j] += weil[order][stlength-oneside_pts-1-i]*coefl*inter_node.wei.weil[i][j];           //Left FP
        }
    }
    //Two fictitious points weight
    for(int i = 0; i < onefp_unknowns; i++)
    {
        A[onefp_unknowns+i][onefp_unknowns+i] = weir[order][stlength-1]*coefr;          //Right FP
        A[onefp_unknowns+i][i] = -weil[order][0]*coefl;                                 //Left FP
    }
    if(abs(inter_node.ID)%2 == 0)                   //"-" => "+", BETA^{-}U_{X}^{-} = BETA^{+}U_{X}^{+} - [BETA U_{X}], F = G - [BETA U_{X}]
    {
        B[total_unknowns-1] += -1;
    }
    else                                            //"+" => "-", BETA^{-}U_{X}^{-} = BETA^{+}U_{X}^{+} - [BETA U_{X}], G + [BETA U_{X}] = F
    {
        B[total_unknowns-1] += 1;
    }
    
    //Using LU decomposition to solve Ax = B
    LU lu_dcmp(A);
    lu_dcmp.solve(B,x);
    
    if(abs(lu_dcmp.det()) < TOL_SETUP)
    {
        cout << "Warning: when solving IRREGULAR interface, A is a singular matrix locate at: X = "
        << inter_node.coord.x_value << " Y = " << inter_node.coord.y_value
        << " Z = " << inter_node.coord.z_value << endl;
        
        cout << "Determinant of A: " << lu_dcmp.det() << endl;
        
        for(int i = 0; i < total_unknowns; i++)
        {
            double sum = 0;
            
            for(int j = 0; j < total_unknowns; j++)
            {
                sum = sum + A[i][j]*x[j];
            }
            
            if(abs(B[i]-sum) > TOL_SETUP)
            {
                cout << "B[" << i << "] approximation error: " << abs(B[i]-sum) << endl;
                
            }
        }
        cout << endl;
    }
}

/***************************************************************************************************************
 Calculating irregular weights for MIB scheme
 
 INPUT
 gamma : distance between intersection node with its left neighbor grid point
 dv    : grid mesh size
 
 OUTPUT
 weil : finite difference weights for left fictitious point
 weir : finite difference weights for right fictitious point
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 GRAPH EXPLAINATION FOR DIFFERENT SITUATION USED IN MIB
 -------- : MESH LINES
 *    : REAL POINTS
 O    : FP
 x    : INTERFACE
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  IRREGULAR INTERFACE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 -----------*---------------*--------x-------*---------------*-------------
 ---------------------------O--------x-------*---------------*-------------            LEFT FP INTERPOLATION
 -----------*---------------*--------x-------O-----------------------------           RIGHT FP INTERPOLATION
 *****************************************************************************************************************/
void Intersections::Get_irr_weights(Doub_I gamma, Doub_I dv, Int_I oneside_pts, MatrixDoub_O& weil, MatrixDoub_O& weir)
{
    VecDoub vl,vr;
    MatrixDoub temp;
    int ns,order;
    
    ns = (int)weil[0].size();
    if((int)weir[0].size() != ns)
    {
        cout << "Different size for input matrix when getting weights" << endl;
        exit(0);
    }
    
    order = (int)weil.size();
    if((int)weir.size() != order)
    {
        cout << "Different size for input matrix when getting weights" << endl;
        exit(0);
    }
    
    vr.resize(ns);
    vl.resize(ns);
    
    temp.resize(ns);
    for(int i = 0; i < ns; i++)
    {
        temp[i].resize(order);
    }
    
    for(int i = 0; i < ns; i++)
    {
        vr[i] = -(oneside_pts-1)*dv + i*dv;
        vl[i] = (oneside_pts-ns+1)*dv + i*dv;
    }
    
    Weights(gamma,vr,ns,order-1,temp);
    for(int i = 0; i < order; i++)
    {
        for(int j = 0; j < ns; j++)
        {
            weir[i][j] = temp[j][i];
        }
    }
    
    Weights(gamma,vl,ns,order-1,temp);
    for(int i = 0; i < order; i++)
    {
        for(int j = 0; j < ns; j++)
        {
            weil[i][j] = temp[j][i];
        }
    }
}

/******************************************************************************
 Calculate four corner fictitous points' weights in 2nd order precision
 
 INPUT
 inter_node_left  : left interseciton node
 inter_node_right : right intersection node
 dv               : grid mesh size
 ******************************************************************************/
void Intersections::Corner_MIB_2nd(Intersection_Data& inter_node_left, Intersection_Data& inter_node_right, Beta& beta, Doub_I dv)
{
    MatrixDoub A, left_wei_out, left_wei_in, right_wei_out, right_wei_in;
    VecDoub B, x;
    double coef_out, coef_in;
    int order;
    int fictitouspts_no, oneside_unknowns, total_unknowns, dev_order, stlength;
    
    fictitouspts_no = 4;                                 //total numbers of fictitious points
    oneside_unknowns = 5 + 2 * 2;                        //numbers of unknowns for one fictitious points
    total_unknowns = 4 * oneside_unknowns;               //total unknowns need to be solved
    dev_order = 1;                                       //highest order of derivative in two equations
    stlength = 4;                                        //length of FD stencil
    
    B.resize(total_unknowns);
    x.resize(total_unknowns);
    A.resize(total_unknowns);
    for(int i = 0; i < total_unknowns; i++)
    {
        A[i].resize(total_unknowns);
    }
    for(int i = 0; i < total_unknowns; i++)
    {
        B[i] = 0;
        for(int j = 0; j < total_unknowns; j++)
        {
            A[i][j] = 0;
        }
    }
    
    left_wei_out.resize(dev_order+1);
    left_wei_in.resize(dev_order+1);
    right_wei_out.resize(dev_order+1);
    right_wei_in.resize(dev_order+1);
    for(int i = 0; i < dev_order+1; i++)
    {
        left_wei_out[i].resize(stlength);
        left_wei_in[i].resize(stlength);
        right_wei_out[i].resize(stlength);
        right_wei_in[i].resize(stlength);
    }
    
    Get_cor_weights(inter_node_left.gamma, inter_node_right.gamma, dv, left_wei_out, left_wei_in, right_wei_out, right_wei_in);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Left Interface ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Zero order: U^{-} = U^{+} - [U]
    order = 0;
    B[0] = -left_wei_out[order][0];         //GP at ix-1
    B[1] = -left_wei_out[order][1];         //GP at ix
    B[2] =  left_wei_in[order][1];          //GP at ix+1
    B[3] = -left_wei_out[order][3];         //GP at ix+2
    for(int i = 0; i < oneside_unknowns; i++)
    {
        A[i][i] = -left_wei_in[order][0];                       //1st FP, at ix
        A[i][oneside_unknowns+i] = left_wei_out[order][2];      //2nd FP, at ix+1
        A[i][oneside_unknowns*2+i] = -left_wei_in[order][2];    //3rd FP, at ix+2
        A[i][oneside_unknowns*3+i] = -left_wei_in[order][3];    //4th FP, at ix-1 or ix+3
    }
    if(abs(inter_node_left.ID)%2 == 0)       //"-" => "+" => "-"
    {
        B[5] = -1;
    }
    else                                //"+" => "-" => "+"
    {
        B[5] = 1;
    }
    
    //First order: BETA^{-}*U_{x} ^{-} = BETA^{+}*U_{x}^{+} - [BETA U_{x}]
    order = 1;
    
    if(abs(inter_node_left.ID)%2 == 0)              //"-" => "+" => "-"
    {
        coef_in = inter_node_left.diffcoef.outside;
        coef_out = inter_node_left.diffcoef.inside;
    }
    else                                           //"+" => "-" => "+"
    {
        coef_in = inter_node_left.diffcoef.inside;
        coef_out = inter_node_left.diffcoef.outside;
    }
    
    B[oneside_unknowns+0] = -left_wei_out[order][0]*coef_out;             //GP at ix-1
    B[oneside_unknowns+1] = -left_wei_out[order][1]*coef_out;             //GP at ix
    B[oneside_unknowns+2] =  left_wei_in[order][1]*coef_in;               //GP at ix+1
    B[oneside_unknowns+3] = -left_wei_out[order][3]*coef_out;             //GP at ix+2
    for(int i = 0; i < oneside_unknowns; i++)
    {
        A[oneside_unknowns+i][i] = -left_wei_in[order][0]*coef_in;                       //1st FP, at ix
        A[oneside_unknowns+i][oneside_unknowns+i] = left_wei_out[order][2]*coef_out;     //2nd FP, at ix+1
        A[oneside_unknowns+i][oneside_unknowns*2+i] = -left_wei_in[order][2]*coef_in;    //3rd FP, at ix+2
        //4th FP, at ix-1 or ix+3
        if((dv-inter_node_left.gamma) > inter_node_right.gamma)                             //4th FP, at ix-1
        {
            A[oneside_unknowns+i][oneside_unknowns*3+i] = -left_wei_in[order][3]*coef_in;
        }
        else                                                                                 //4th FP, at ix+3
        {
            A[oneside_unknowns+i][oneside_unknowns*3+i] = -left_wei_in[order][3]*coef_in;
        }
    }
    if(abs(inter_node_left.ID)%2 == 0)       //"-" => "+" => "-"
    {
        B[oneside_unknowns+6] = -1;
    }
    else                                     //"+" => "-" => "+"
    {
        B[oneside_unknowns+6] = 1;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Right Interface ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Zero order: U^{-} = U^{+} - [U]
    order = 0;
    B[oneside_unknowns*2+1] = -right_wei_out[order][0];            //GP at ix
    B[oneside_unknowns*2+2] =  right_wei_in[order][1];             //GP at ix+1
    B[oneside_unknowns*2+3] = -right_wei_out[order][2];            //GP at ix+2
    B[oneside_unknowns*2+4] = -right_wei_out[order][3];            //GP at ix+3
    for(int i = 0; i < oneside_unknowns; i++)
    {
        A[oneside_unknowns*2+i][i] = -right_wei_in[order][0];                       //1st FP, at ix
        A[oneside_unknowns*2+i][oneside_unknowns+i] = right_wei_out[order][1];      //2nd FP, at ix+1
        A[oneside_unknowns*2+i][oneside_unknowns*2+i] = -right_wei_in[order][2];    //3rd FP, at ix+2
        A[oneside_unknowns*2+i][oneside_unknowns*3+i] = -right_wei_in[order][3];    //4th FP, at ix-1 or ix+3
    }
    if(abs(inter_node_left.ID)%2 == 0)       //"-" => "+" => "-"
    {
        B[oneside_unknowns*2+7] = -1;
    }
    else                                    //"+" => "-" => "+"
    {
        B[oneside_unknowns*2+7] = 1;
    }
    
    //First order: BETA^{-}*U_{x} ^{-} = BETA^{+}*U_{x}^{+} - [BETA U_{x}]
    order = 1;
    if(abs(inter_node_left.ID)%2 == 0)              //"-" => "+" => "-"
    {
        coef_in = inter_node_right.diffcoef.outside;
        coef_out = inter_node_right.diffcoef.inside;
    }
    else                                            //"+" => "-" => "+"
    {
        coef_in = inter_node_right.diffcoef.inside;
        coef_out = inter_node_right.diffcoef.outside;
    }
    
    B[oneside_unknowns*3+1] = -right_wei_out[order][0]*coef_out;            //GP at ix
    B[oneside_unknowns*3+2] =  right_wei_in[order][1]*coef_in;              //GP at ix+1
    B[oneside_unknowns*3+3] = -right_wei_out[order][2]*coef_out;            //GP at ix+2
    B[oneside_unknowns*3+4] = -right_wei_out[order][3]*coef_out;            //GP at ix+3
    for(int i = 0; i < oneside_unknowns; i++)
    {
        A[oneside_unknowns*3+i][i] = -right_wei_in[order][0]*coef_in;                       //1st FP, at ix
        A[oneside_unknowns*3+i][oneside_unknowns+i] = right_wei_out[order][1]*coef_out;     //2nd FP, at ix+1
        A[oneside_unknowns*3+i][oneside_unknowns*2+i] = -right_wei_in[order][2]*coef_in;    //3rd FP, at ix+2
        //4th FP, at ix-1 or ix+3
        if((dv-inter_node_left.gamma) > inter_node_right.gamma)                                //4th FP, at ix-1
        {
            A[oneside_unknowns*3+i][oneside_unknowns*3+i] = -right_wei_in[order][3]*coef_in;
        }
        else                                                                                    //4th FP, at ix+3
        {
            A[oneside_unknowns*3+i][oneside_unknowns*3+i] = -right_wei_in[order][3]*coef_in;
        }
    }
    if(abs(inter_node_left.ID)%2 == 0)       //"-" => "+" => "-"
    {
        B[oneside_unknowns*3+8] = -1;
    }
    else                                    //"+" => "-" => "+"
    {
        B[oneside_unknowns*3+8] = 1;
    }
    
    //Using LU decomposition to solve Ax = B
    LU lu_dcmp(A);
    lu_dcmp.solve(B,x);
    
    if(abs(lu_dcmp.det()) < TOL_SETUP)
    {
        cout << "Warning: when solving CORNER interface, A is a singular matrix locate at: X = "
        << inter_node_left.coord.x_value << " Y = " << inter_node_left.coord.y_value
        << " Z = " << inter_node_left.coord.z_value << endl;
        
        cout << "Determinant of A: " << lu_dcmp.det() << endl;
        
        for(int i = 0; i < total_unknowns; i++)
        {
            double sum = 0;
            
            for(int j = 0; j < total_unknowns; j++)
            {
                sum = sum + A[i][j]*x[j];
            }
            
            if(abs(B[i]-sum) > TOL_SETUP)
            {
                cout << "B[" << i << "] approximation error: " << abs(B[i]-sum) << endl;
                
            }
        }
        cout << endl;
    }
    
    inter_node_left.wei.weil.resize(1);
    inter_node_left.wei.weir.resize(1);
    inter_node_right.wei.weil.resize(1);
    inter_node_right.wei.weir.resize(1);
    for(int i = 0; i < 1; i++)
    {
        inter_node_left.wei.weil[i].resize(oneside_unknowns);
        inter_node_left.wei.weir[i].resize(oneside_unknowns);
        inter_node_right.wei.weil[i].resize(oneside_unknowns);
        inter_node_right.wei.weir[i].resize(oneside_unknowns);
    }
    
    for(int i = 0; i < 1; i++)
    {
        for(int j = 0; j < oneside_unknowns; j++)
        {
            inter_node_left.wei.weil[i][j] = x[j];
            inter_node_left.wei.weir[i][j] = x[oneside_unknowns+j];
            inter_node_right.wei.weil[i][j] = x[oneside_unknowns+j];
            inter_node_right.wei.weir[i][j] = x[oneside_unknowns*2+j];
        }
    }
}

/**********************************************************************************************************************
 Calculating irregular weights for MIB scheme(with retreated type)
 
 INPUT
 gamma_left  : distance between left intersection node with its left neighbor grid point
 gamma_right : distance between right intersection node with its left neighbor grid point
 dv          : grid mesh size
 
 OUTPUT
 left_wei_out  : weights used for outside interpolation on left interface
 left_wei_in   : weights used for inside interpolation on left interface
 right_wei_out : weights used for outside interpolation on right interface
 right_wei_in  : weights used for inside interpolation on right interface
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 GRAPH EXPLAINATION FOR DIFFERENT SITUATION USED IN MIB
 -------- : MESH LINES
 *    : REAL POINTS
 O    : FP
 x    : INTERFACE
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  CORNER LEFT INTERFACE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     (itype_left,itype_right)
 ---------*---------------*-------x--------*--------x-------*---------------*----------
 OUTSIDE INTERPOLATION, left_v_out
 ---------*---------------*-------x--------O--------x-------*--------------------------
 INSIDE INTERPOLATION, left_v_in
 ---------O---------------O---x------------*--------x-------O--------------------------             ( -2 , 0 )
 -------------------------O-------x--------*------------x---O---------------O----------             (  0 , 2 )
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  CORNER RIGHT INTERFACE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     (itype_left,itype_right)
 ---------*---------------*-------x--------*--------x-------*---------------*----------
 OUTSIDE INTERPOLATION, right_v_out
 -------------------------*-------x--------O--------x-------*---------------*----------
 INSIDE INTERPOLATION, right_v_in
 ---------O---------------O---x------------*--------x-------O--------------------------             ( -2 , 0 )
 -------------------------O-------x--------*------------x---O---------------O----------             (  0 , 2 )
 ********************************************************************************************************************/
void Intersections::Get_cor_weights(Doub_I gamma_left, Doub_I gamma_right, Doub_I dv, MatrixDoub_O& left_wei_out, MatrixDoub_O& left_wei_in, MatrixDoub_O& right_wei_out, MatrixDoub_O& right_wei_in)
{
    VecDoub left_v_out,left_v_in,right_v_out,right_v_in;
    MatrixDoub temp;
    int ns,order;
 
    ns = (int)left_wei_out[0].size();
    if((int)left_wei_in[0].size() != ns || (int)right_wei_out[0].size() != ns || (int)right_wei_in[0].size() != ns)
    {
        cout << "Different size for input matrix when getting weights" << endl;
        exit(0);
    }
    
    order = (int)left_wei_out.size();
    if((int)left_wei_in.size() != order || (int)right_wei_out.size() != order || (int)right_wei_in.size() != order)
    {
        cout << "Different size for input matrix when getting weights" << endl;
        exit(0);
    }
    
    left_v_out.resize(ns);
    left_v_in.resize(ns);
    right_v_out.resize(ns);
    right_v_in.resize(ns);
    
    temp.resize(ns);
    for(int i = 0; i < ns; i++)
    {
        temp[i].resize(order);
    }
    
    for(int i = 0; i < ns; i++)
    {
        left_v_out[i] = -dv + i*dv;
        right_v_out[i] = i*dv;
    }
    
    for(int i = 0; i < ns; i++) //(0,2)
    {
        left_v_in[i] = i*dv;
    }
    if((dv-gamma_left) > gamma_right)
    {
        left_v_in[ns-1] = -dv;  //(-2,0), center point close to right interface, 4th FP changes
    }
    
    for(int i = 0; i < ns; i++)
    {
        right_v_in[i] = left_v_in[i];
    }
    
    Weights(gamma_left,left_v_out,ns,order-1,temp);
    for(int i = 0; i < order; i++)
    {
        for(int j = 0; j < ns; j++)
        {
            left_wei_out[i][j] = temp[j][i];
        }
    }
    
    Weights(gamma_left,left_v_in,ns,order-1,temp);
    for(int i = 0; i < order; i++)
    {
        for(int j = 0; j < ns; j++)
        {
            left_wei_in[i][j] = temp[j][i];
        }
    }
    
    Weights(dv+gamma_right,right_v_out,ns,order-1,temp);
    for(int i = 0; i < order; i++)
    {
        for(int j = 0; j < ns; j++)
        {
            right_wei_out[i][j] = temp[j][i];
        }
    }
    
    Weights(dv+gamma_right,right_v_in,ns,order-1,temp);
    for(int i = 0; i < order; i++)
    {
        for(int j = 0; j < ns; j++)
        {
            right_wei_in[i][j] = temp[j][i];
        }
    }
}

/******************************************************************************
 Calculate error at current time w.r.t given equation object
 
 INPUT
 eq : equation object
 
 OUTPUT
 fictitous points errors for all intersection points
 ******************************************************************************/
void Intersections::Refresh_Fp(Equation& eq)
{
    Intersection_Data inter_node, inter_node_left, inter_node_right;
    CubicDoub u,invu;
    int ix, iy, iz;
    int n, ip, indx, fpno, oneside_pts;
    double suml, sumr, left_suml, left_sumr, right_sumr;
    double x, y, z, left_x, right_x, left_y, right_y, left_z, right_z;
    
    u.resize(nx);
    invu.resize(nx);
    for(int i = 0; i < nx; i++)
    {
        u[i].resize(ny);
        invu[i].resize(ny);
        for(int j = 0; j < ny; j++)
        {
            u[i][j].resize(nz);
            invu[i][j].resize(nz);
        }
    }
    
    for(ix = 0; ix < nx; ix++)
    {
        for(iy = 0; iy < ny; iy++)
        {
            for(iz = 0; iz < nz; iz++)
            {
                x = xi[ix];
                y = yi[iy];
                z = zi[iz];
                
                indx = To1d(ix,iy,iz);
                if(indicator[indx] == 1)
                {
                    u[ix][iy][iz] = eq.Outer_u(x,y,z);
                    invu[ix][iy][iz] = eq.Inner_u(x,y,z);
                }
                else if(indicator[indx] == -1)
                {
                    u[ix][iy][iz] = eq.Inner_u(x,y,z);
                    invu[ix][iy][iz] = eq.Outer_u(x,y,z);
                }
                else
                {
                    cout << "Wrong value for mesh setup" << endl;
                    exit(0);
                }
            }
        }
    }
    
    //Z-direction
    for(ix = 1; ix < nx-1; ix++)
    {
        for(iy = 1; iy < ny-1; iy++)
        {
            ip = 0;
            while(ip < ifpz[ix][iy].size())
            {
                if(ifpz[ix][iy][ip].ID > 0)
                {
                    inter_node = ifpz[ix][iy][ip];
                    
                    fpno = (int)inter_node.err.errr.size();
                    if(fpno != (int)inter_node.err.errl.size())
                    {
                        cout << "Bad size for irregular MIB error" << endl;
                        exit(0);
                    }
                    
                    for(int i = 0; i < fpno; i++)
                    {
                        n = (int)inter_node.wei.weil[i].size();
                        if(n == (int)inter_node.wei.weir[i].size())
                        {
                            n = n-2;
                        }
                        else
                        {
                            cout << "Bad size for irregular MIB weights" << endl;
                            exit(0);
                        }
                        
                        oneside_pts = n/2;
                        
                        x = inter_node.coord.x_value;
                        y = inter_node.coord.y_value;
                        z = inter_node.coord.z_value;
                        
                        iz = inter_node.left_loc;
                        
                        suml = 0;
                        sumr = 0;
                        
                        suml += inter_node.wei.weil[i][n] * eq.Jump_u(x,y,z);
                        suml += inter_node.wei.weil[i][n+1] * eq.Jump_betau_z(x,y,z);
                        sumr += inter_node.wei.weir[i][n] * eq.Jump_u(x,y,z);
                        sumr += inter_node.wei.weir[i][n+1] * eq.Jump_betau_z(x,y,z);
                        
                        for(int j = 0; j < n; j++)
                        {
                            suml += inter_node.wei.weil[i][j] * u[ix][iy][iz-oneside_pts+1+j];
                            sumr += inter_node.wei.weir[i][j] * u[ix][iy][iz-oneside_pts+1+j];
                        }
                        
                        ifpz[ix][iy][ip].err.errl[i] = abs((invu[ix][iy][iz-i] - suml)/invu[ix][iy][iz-i]);
                        ifpz[ix][iy][ip].err.errr[i] = abs((invu[ix][iy][iz+1+i] - sumr)/invu[ix][iy][iz+1+i]);
                    }
                    
                    ip += 1;
                    
                }
                else //Corner interface
                {
                    inter_node_left = ifpz[ix][iy][ip];
                    inter_node_right = ifpz[ix][iy][ip+1];
                    
                    n = (int)inter_node_left.wei.weil[0].size();
                    if(n == (int)inter_node_left.wei.weir[0].size() &&
                       n == (int)inter_node_right.wei.weil[0].size() && n == (int)inter_node_right.wei.weir[0].size())
                    {
                        n = n - 2*2;
                    }
                    else
                    {
                        cout << "Bad size for corner MIB weights" << endl;
                    }
                    
                    left_x = inter_node_left.coord.x_value;
                    left_y = inter_node_left.coord.y_value;
                    left_z = inter_node_left.coord.z_value;
                    
                    right_x = inter_node_right.coord.x_value;
                    right_y = inter_node_right.coord.y_value;
                    right_z = inter_node_right.coord.z_value;
                    
                    iz = inter_node_left.left_loc;
                    
                    left_suml = 0;
                    left_sumr = 0;
                    right_sumr = 0;
                    
                    left_suml += inter_node_left.wei.weil[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    left_suml += inter_node_left.wei.weil[0][n+1] * eq.Jump_betau_z(left_x,left_y,left_z);
                    left_suml += inter_node_left.wei.weil[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    left_suml += inter_node_left.wei.weil[0][n+3] * eq.Jump_betau_z(right_x,right_y,right_z);
                    
                    left_sumr += inter_node_left.wei.weir[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    left_sumr += inter_node_left.wei.weir[0][n+1] * eq.Jump_betau_z(left_x,left_y,left_z);
                    left_sumr += inter_node_left.wei.weir[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    left_sumr += inter_node_left.wei.weir[0][n+3] * eq.Jump_betau_z(right_x,right_y,right_z);
                    
                    right_sumr += inter_node_right.wei.weir[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    right_sumr += inter_node_right.wei.weir[0][n+1] * eq.Jump_betau_z(left_x,left_y,left_z);
                    right_sumr += inter_node_right.wei.weir[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    right_sumr += inter_node_right.wei.weir[0][n+3] * eq.Jump_betau_z(right_x,right_y,right_z);
                    
                    for(int i = 0; i < n; i++)
                    {
                        left_suml += inter_node_left.wei.weil[0][i] * u[ix][iy][iz-1+i];
                        left_sumr += inter_node_left.wei.weir[0][i] * u[ix][iy][iz-1+i];
                        right_sumr += inter_node_right.wei.weir[0][i] * u[ix][iy][iz-1+i];
                    }
                    
                    //Left corner interface
                    ifpz[ix][iy][ip].err.errl[0] = abs((invu[ix][iy][iz]-left_suml)/invu[ix][iy][iz]);
                    ifpz[ix][iy][ip].err.errr[0] = abs((invu[ix][iy][iz+1]-left_sumr)/invu[ix][iy][iz+1]);
                    
                    //Right corner interface
                    ifpz[ix][iy][ip+1].err.errl[0] = ifpz[ix][iy][ip].err.errr[0];
                    ifpz[ix][iy][ip+1].err.errr[0] = abs((invu[ix][iy][iz+2]-right_sumr)/invu[ix][iy][iz+2]);
                    
                    ip += 2;
                }
            }
        }
    }
    
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            ip = 0;
            while(ip < ifpy[ix][iz].size())
            {
                if(ifpy[ix][iz][ip].ID > 0)
                {
                    inter_node = ifpy[ix][iz][ip];
                    
                    fpno = (int)inter_node.err.errr.size();
                    if(fpno != (int)inter_node.err.errl.size())
                    {
                        cout << "Bad size for irregular MIB error" << endl;
                        exit(0);
                    }
                    
                    for(int i = 0; i < fpno; i++)
                    {
                        n = (int)inter_node.wei.weil[i].size();
                        if(n == (int)inter_node.wei.weir[i].size())
                        {
                            n = n-2;
                        }
                        else
                        {
                            cout << "Bad size for irregular MIB weights" << endl;
                            exit(0);
                        }
                        
                        oneside_pts = n/2;
                        
                        x = inter_node.coord.x_value;
                        y = inter_node.coord.y_value;
                        z = inter_node.coord.z_value;
                        
                        iy = inter_node.left_loc;
                        
                        suml = 0;
                        sumr = 0;
                        
                        suml += inter_node.wei.weil[i][n] * eq.Jump_u(x,y,z);
                        suml += inter_node.wei.weil[i][n+1] * eq.Jump_betau_y(x,y,z);
                        sumr += inter_node.wei.weir[i][n] * eq.Jump_u(x,y,z);
                        sumr += inter_node.wei.weir[i][n+1] * eq.Jump_betau_y(x,y,z);
                        
                        for(int j = 0; j < n; j++)
                        {
                            suml += inter_node.wei.weil[i][j] * u[ix][iy-oneside_pts+1+j][iz];
                            sumr += inter_node.wei.weir[i][j] * u[ix][iy-oneside_pts+1+j][iz];;
                        }
                        
                        ifpy[ix][iz][ip].err.errl[i] = abs((invu[ix][iy-i][iz] - suml)/invu[ix][iy-i][iz]);
                        ifpy[ix][iz][ip].err.errr[i] = abs((invu[ix][iy+1+i][iz] - sumr)/invu[ix][iy+1+i][iz]);
                    }
                    
                    ip += 1;
                    
                }
                else //Corner interface
                {
                    inter_node_left = ifpy[ix][iz][ip];
                    inter_node_right = ifpy[ix][iz][ip+1];
                    
                    n = (int)inter_node_left.wei.weil[0].size();
                    if(n == (int)inter_node_left.wei.weir[0].size() &&
                       n == (int)inter_node_right.wei.weil[0].size() && n == (int)inter_node_right.wei.weir[0].size())
                    {
                        n = n - 2*2;
                    }
                    else
                    {
                        cout << "Bad size for corner MIB weights" << endl;
                    }
                    
                    left_x = inter_node_left.coord.x_value;
                    left_y = inter_node_left.coord.y_value;
                    left_z = inter_node_left.coord.z_value;
                    
                    right_x = inter_node_right.coord.x_value;
                    right_y = inter_node_right.coord.y_value;
                    right_z = inter_node_right.coord.z_value;
                    
                    iy = inter_node_left.left_loc;
                    
                    left_suml = 0;
                    left_sumr = 0;
                    right_sumr = 0;
                    
                    left_suml += inter_node_left.wei.weil[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    left_suml += inter_node_left.wei.weil[0][n+1] * eq.Jump_betau_y(left_x,left_y,left_z);
                    left_suml += inter_node_left.wei.weil[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    left_suml += inter_node_left.wei.weil[0][n+3] * eq.Jump_betau_y(right_x,right_y,right_z);
                    
                    left_sumr += inter_node_left.wei.weir[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    left_sumr += inter_node_left.wei.weir[0][n+1] * eq.Jump_betau_y(left_x,left_y,left_z);
                    left_sumr += inter_node_left.wei.weir[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    left_sumr += inter_node_left.wei.weir[0][n+3] * eq.Jump_betau_y(right_x,right_y,right_z);
                    
                    right_sumr += inter_node_right.wei.weir[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    right_sumr += inter_node_right.wei.weir[0][n+1] * eq.Jump_betau_y(left_x,left_y,left_z);
                    right_sumr += inter_node_right.wei.weir[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    right_sumr += inter_node_right.wei.weir[0][n+3] * eq.Jump_betau_y(right_x,right_y,right_z);
                    
                    for(int i = 0; i < n; i++)
                    {
                        left_suml += inter_node_left.wei.weil[0][i] * u[ix][iy-1+i][iz];
                        left_sumr += inter_node_left.wei.weir[0][i] * u[ix][iy-1+i][iz];
                        right_sumr += inter_node_right.wei.weir[0][i] * u[ix][iy-1+i][iz];
                    }
                    
                    //Left corner interface
                    ifpy[ix][iz][ip].err.errl[0] = abs((invu[ix][iy][iz]-left_suml)/invu[ix][iy][iz]);
                    ifpy[ix][iz][ip].err.errr[0] = abs((invu[ix][iy+1][iz]-left_sumr)/invu[ix][iy+1][iz]);
                    
                    //Right corner interface
                    ifpy[ix][iz][ip+1].err.errl[0] = ifpy[ix][iz][ip].err.errr[0];
                    ifpy[ix][iz][ip+1].err.errr[0] = abs((invu[ix][iy+2][iz]-right_sumr)/invu[ix][iy+2][iz]);
                    
                    ip += 2;
                }
            }
        }
    }
    
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            ip = 0;
            while(ip < ifpx[iy][iz].size())
            {
                if(ifpx[iy][iz][ip].ID > 0)
                {
                    inter_node = ifpx[iy][iz][ip];
                    
                    fpno = (int)inter_node.err.errr.size();
                    if(fpno != (int)inter_node.err.errl.size())
                    {
                        cout << "Bad size for irregular MIB error" << endl;
                        exit(0);
                    }
                    
                    for(int i = 0; i < fpno; i++)
                    {
                        n = (int)inter_node.wei.weil[i].size();
                        if(n == (int)inter_node.wei.weir[i].size())
                        {
                            n = n-2;
                        }
                        else
                        {
                            cout << "Bad size for irregular MIB weights" << endl;
                            exit(0);
                        }
                        
                        oneside_pts = n/2;
                        
                        x = inter_node.coord.x_value;
                        y = inter_node.coord.y_value;
                        z = inter_node.coord.z_value;
                        
                        ix = inter_node.left_loc;
                        
                        suml = 0;
                        sumr = 0;
                        
                        suml += inter_node.wei.weil[i][n]*eq.Jump_u(x,y,z);
                        suml += inter_node.wei.weil[i][n+1]*eq.Jump_betau_x(x,y,z);
                        sumr += inter_node.wei.weir[i][n]*eq.Jump_u(x,y,z);
                        sumr += inter_node.wei.weir[i][n+1]*eq.Jump_betau_x(x,y,z);
                        
                        for(int j = 0; j < n; j++)
                        {
                            suml += inter_node.wei.weil[i][j] * u[ix-oneside_pts+1+j][iy][iz];
                            sumr += inter_node.wei.weir[i][j] * u[ix-oneside_pts+1+j][iy][iz];
                        }
                        
                        ifpx[iy][iz][ip].err.errl[i] = abs((invu[ix-i][iy][iz] - suml)/invu[ix-i][iy][iz]);
                        ifpx[iy][iz][ip].err.errr[i] = abs((invu[ix+1+i][iy][iz] - sumr)/invu[ix+1+i][iy][iz]);
                    }
                    
                    ip += 1;
                }
                else //Corner interface
                {
                    inter_node_left = ifpx[iy][iz][ip];
                    inter_node_right = ifpx[iy][iz][ip+1];
                    
                    n = (int)inter_node_left.wei.weil[0].size();
                    if(n == (int)inter_node_left.wei.weir[0].size() &&
                       n == (int)inter_node_right.wei.weil[0].size() && n == (int)inter_node_right.wei.weir[0].size())
                    {
                        n = n - 2*2;
                    }
                    else
                    {
                        cout << "Bad size for corner MIB weights" << endl;
                    }
                    
                    left_x = inter_node_left.coord.x_value;
                    left_y = inter_node_left.coord.y_value;
                    left_z = inter_node_left.coord.z_value;
                    
                    right_x = inter_node_right.coord.x_value;
                    right_y = inter_node_right.coord.y_value;
                    right_z = inter_node_right.coord.z_value;
                    
                    ix = inter_node_left.left_loc;
                    
                    left_suml = 0;
                    left_sumr = 0;
                    right_sumr = 0;
                    
                    left_suml += inter_node_left.wei.weil[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    left_suml += inter_node_left.wei.weil[0][n+1] * eq.Jump_betau_x(left_x,left_y,left_z);
                    left_suml += inter_node_left.wei.weil[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    left_suml += inter_node_left.wei.weil[0][n+3] * eq.Jump_betau_x(right_x,right_y,right_z);
                    
                    left_sumr += inter_node_left.wei.weir[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    left_sumr += inter_node_left.wei.weir[0][n+1] * eq.Jump_betau_x(left_x,left_y,left_z);
                    left_sumr += inter_node_left.wei.weir[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    left_sumr += inter_node_left.wei.weir[0][n+3] * eq.Jump_betau_x(right_x,right_y,right_z);
                    
                    right_sumr += inter_node_right.wei.weir[0][n] * eq.Jump_u(left_x,left_y,left_z);
                    right_sumr += inter_node_right.wei.weir[0][n+1] * eq.Jump_betau_x(left_x,left_y,left_z);
                    right_sumr += inter_node_right.wei.weir[0][n+2] * eq.Jump_u(right_x,right_y,right_z);
                    right_sumr += inter_node_right.wei.weir[0][n+3] * eq.Jump_betau_x(right_x,right_y,right_z);
                    
                    for(int i = 0; i < n; i++)
                    {
                        left_suml += inter_node_left.wei.weil[0][i] * u[ix-1+i][iy][iz];
                        left_sumr += inter_node_left.wei.weir[0][i] * u[ix-1+i][iy][iz];
                        right_sumr += inter_node_right.wei.weir[0][i] * u[ix-1+i][iy][iz];
                    }
                    
                    //Left corner interface
                    ifpx[iy][iz][ip].err.errl[0] = abs((invu[ix][iy][iz]-left_suml)/invu[ix][iy][iz]);
                    ifpx[iy][iz][ip].err.errr[0] = abs((invu[ix+1][iy][iz]-left_sumr)/invu[ix+1][iy][iz]);
                    
                    //Right corner interface
                    ifpx[iy][iz][ip+1].err.errl[0] = ifpx[iy][iz][ip].err.errr[0];
                    ifpx[iy][iz][ip+1].err.errr[0] = abs((invu[ix+2][iy][iz]-right_sumr)/invu[ix+2][iy][iz]);
                    
                    ip += 2;
                }
            }
        }
    }
}

/******************************************************************************
 Calculate approximated jumps at current time w.r.t given equation object
 
 INPUT
 eq : equation object
 uh : solution at current time
 
 OUTPUT
 inter_node.jump : accurate jumps at current time
 ******************************************************************************/
void Intersections::Refresh_Jump(Equation& eq, CubicDoub& uh, Beta& beta)
{
    double jumpbeta_xi, jumpbeta_tau, jumpbeta_eta;
    double coord_x, coord_y, coord_z;
    double eta, tau;
    MatrixDoub p;
    
    p.resize(3);
    for(int i = 0; i < 3; i++)
    {
        p[i].resize(3);
    }
    
    //Z-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iy = 1; iy < ny-1; iy++)
        {
            for(int ip = 0; ip < ifpz[ix][iy].size(); ip++)
            {
                coord_x = ifpz[ix][iy][ip].coord.x_value;
                coord_y = ifpz[ix][iy][ip].coord.y_value;
                coord_z = ifpz[ix][iy][ip].coord.z_value;
                
                ifpz[ix][iy][ip].jump.u = eq.Jump_u(coord_x,coord_y,coord_z);
                ifpz[ix][iy][ip].jump.betau_xi = eq.Jump_betau_xi(coord_x,coord_y,coord_z,
                                                                  ifpz[ix][iy][ip].p[0][0],ifpz[ix][iy][ip].p[0][1],ifpz[ix][iy][ip].p[0][2]);
                ifpz[ix][iy][ip].jump.u_eta = eq.Jump_u_eta(coord_x,coord_y,coord_z,
                                                            ifpz[ix][iy][ip].p[1][0],ifpz[ix][iy][ip].p[1][1],ifpz[ix][iy][ip].p[1][2]);
                ifpz[ix][iy][ip].jump.u_tau = eq.Jump_u_tau(coord_x,coord_y,coord_z,
                                                            ifpz[ix][iy][ip].p[2][0],ifpz[ix][iy][ip].p[2][1],ifpz[ix][iy][ip].p[2][2]);
                
                jumpbeta_xi = ifpz[ix][iy][ip].jump.betau_xi;
                jumpbeta_eta = Eta_z(ifpz[ix][iy][ip],uh,beta);
                jumpbeta_tau = Tau_z(ifpz[ix][iy][ip],uh,beta);
                ifpz[ix][iy][ip].jump.betau_eta = jumpbeta_eta;
                ifpz[ix][iy][ip].jump.betau_tau = jumpbeta_tau;
                
                
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        p[i][j] = ifpz[ix][iy][ip].p[i][j];
                    }
                }
                
                eta = p[1][0]*eq.Jump_betau_x(coord_x,coord_y,coord_z)+p[1][1]*eq.Jump_betau_y(coord_x,coord_y,coord_z)+p[1][2]*eq.Jump_betau_z(coord_x,coord_y,coord_z);
                ifpz[ix][iy][ip].jump.eta_err = abs(jumpbeta_eta-eta)/eta;
                tau = p[2][0]*eq.Jump_betau_x(coord_x,coord_y,coord_z)+p[2][1]*eq.Jump_betau_y(coord_x,coord_y,coord_z)+p[2][2]*eq.Jump_betau_z(coord_x,coord_y,coord_z);
                ifpz[ix][iy][ip].jump.tau_err = abs(jumpbeta_tau-tau)/tau;
                
                LU lu_inv(p);
                lu_inv.inverse(p);
                
                ifpz[ix][iy][ip].jump.u_dir = p[2][0]*jumpbeta_xi+p[2][1]*jumpbeta_eta+p[2][2]*jumpbeta_tau;
                
                ifpz[ix][iy][ip].jump.err = abs(ifpz[ix][iy][ip].jump.u_dir-eq.Jump_betau_z(coord_x,coord_y,coord_z))/eq.Jump_betau_z(coord_x,coord_y,coord_z);
            }
        }
    }
    
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpy[ix][iz].size(); ip++)
            {
                coord_x = ifpy[ix][iz][ip].coord.x_value;
                coord_y = ifpy[ix][iz][ip].coord.y_value;
                coord_z = ifpy[ix][iz][ip].coord.z_value;
                
                ifpy[ix][iz][ip].jump.u = eq.Jump_u(coord_x,coord_y,coord_z);
                ifpy[ix][iz][ip].jump.betau_xi = eq.Jump_betau_xi(coord_x,coord_y,coord_z,ifpy[ix][iz][ip].p[0][0],ifpy[ix][iz][ip].p[0][1],
                                                                  ifpy[ix][iz][ip].p[0][2]);
                ifpy[ix][iz][ip].jump.u_eta = eq.Jump_u_eta(coord_x,coord_y,coord_z,ifpy[ix][iz][ip].p[1][0],ifpy[ix][iz][ip].p[1][1],
                                                            ifpy[ix][iz][ip].p[1][2]);
                ifpy[ix][iz][ip].jump.u_tau = eq.Jump_u_tau(coord_x,coord_y,coord_z,ifpy[ix][iz][ip].p[2][0],ifpy[ix][iz][ip].p[2][1],
                                                            ifpy[ix][iz][ip].p[2][2]);
                
                jumpbeta_xi = ifpy[ix][iz][ip].jump.betau_xi;
                jumpbeta_eta = Eta_y(ifpy[ix][iz][ip],uh,beta);
                jumpbeta_tau = Tau_y(ifpy[ix][iz][ip],uh,beta);
                ifpy[ix][iz][ip].jump.betau_eta = jumpbeta_eta;
                ifpy[ix][iz][ip].jump.betau_tau = jumpbeta_tau;
                
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        p[i][j] = ifpy[ix][iz][ip].p[i][j];
                    }
                }
                
                eta = p[1][0]*eq.Jump_betau_x(coord_x,coord_y,coord_z)+p[1][1]*eq.Jump_betau_y(coord_x,coord_y,coord_z)+p[1][2]*eq.Jump_betau_z(coord_x,coord_y,coord_z);
                ifpy[ix][iz][ip].jump.eta_err = abs(jumpbeta_eta-eta)/eta;
                tau = p[2][0]*eq.Jump_betau_x(coord_x,coord_y,coord_z)+p[2][1]*eq.Jump_betau_y(coord_x,coord_y,coord_z)+p[2][2]*eq.Jump_betau_z(coord_x,coord_y,coord_z);
                ifpy[ix][iz][ip].jump.tau_err = abs(jumpbeta_tau-tau)/tau;
                
                LU lu_inv(p);
                lu_inv.inverse(p);
                
                ifpy[ix][iz][ip].jump.u_dir = p[1][0]*jumpbeta_xi+p[1][1]*jumpbeta_eta+p[1][2]*jumpbeta_tau;
                
                ifpy[ix][iz][ip].jump.err = abs(ifpy[ix][iz][ip].jump.u_dir-eq.Jump_betau_y(coord_x,coord_y,coord_z))/eq.Jump_betau_y(coord_x,coord_y,coord_z);
            }
        }
    }
    
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpx[iy][iz].size(); ip++)
            {
                coord_x = ifpx[iy][iz][ip].coord.x_value;
                coord_y = ifpx[iy][iz][ip].coord.y_value;
                coord_z = ifpx[iy][iz][ip].coord.z_value;
                
                ifpx[iy][iz][ip].jump.u = eq.Jump_u(coord_x,coord_y,coord_z);
                ifpx[iy][iz][ip].jump.betau_xi = eq.Jump_betau_xi(coord_x,coord_y,coord_z,ifpx[iy][iz][ip].p[0][0],ifpx[iy][iz][ip].p[0][1],ifpx[iy][iz][ip].p[0][2]);
                ifpx[iy][iz][ip].jump.u_eta = eq.Jump_u_eta(coord_x,coord_y,coord_z,ifpx[iy][iz][ip].p[1][0],ifpx[iy][iz][ip].p[1][1],
                                                            ifpx[iy][iz][ip].p[1][2]);
                ifpx[iy][iz][ip].jump.u_tau = eq.Jump_u_tau(coord_x,coord_y,coord_z,ifpx[iy][iz][ip].p[2][0],ifpx[iy][iz][ip].p[2][1],
                                                            ifpx[iy][iz][ip].p[2][2]);
                
                jumpbeta_xi = ifpx[iy][iz][ip].jump.betau_xi;
                jumpbeta_eta = Eta_x(ifpx[iy][iz][ip],uh,beta);
                jumpbeta_tau = Tau_x(ifpx[iy][iz][ip],uh,beta);
                ifpx[iy][iz][ip].jump.betau_eta = jumpbeta_eta;
                ifpx[iy][iz][ip].jump.betau_tau = jumpbeta_tau;
                
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        p[i][j] = ifpx[iy][iz][ip].p[i][j];
                    }
                }
                
                eta = p[1][0]*eq.Jump_betau_x(coord_x,coord_y,coord_z)+p[1][1]*eq.Jump_betau_y(coord_x,coord_y,coord_z)+p[1][2]*eq.Jump_betau_z(coord_x,coord_y,coord_z);
                ifpx[iy][iz][ip].jump.eta_err = abs(jumpbeta_eta-eta)/eta;
                tau = p[2][0]*eq.Jump_betau_x(coord_x,coord_y,coord_z)+p[2][1]*eq.Jump_betau_y(coord_x,coord_y,coord_z)+p[2][2]*eq.Jump_betau_z(coord_x,coord_y,coord_z);
                ifpx[iy][iz][ip].jump.tau_err = abs(jumpbeta_tau-tau)/tau;
                
                LU lu_inv(p);
                lu_inv.inverse(p);
                
                ifpx[iy][iz][ip].jump.u_dir = p[0][0]*jumpbeta_xi+p[0][1]*jumpbeta_eta+p[0][2]*jumpbeta_tau;
                
                ifpx[iy][iz][ip].jump.err = abs(ifpx[iy][iz][ip].jump.u_dir-eq.Jump_betau_x(coord_x,coord_y,coord_z))/eq.Jump_betau_x(coord_x,coord_y,coord_z);
            }
        }
    }
}

/************************************************************************************************************************
 Normalized local coordinate transformation matrix initialization at given intersection node on yz plane in x-direction
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.p : local coordinate transformation matrix for given intersection node
 ************************************************************************************************************************/
void Intersections::Pmatrix_Setup_x(Intersection_Data& inter_node)
{
    VecDoub a;
    
    for(int i = 0; i < 3; i++)
    {
        a.push_back(inter_node.local.normal[i]);
    }
    
    inter_node.p.resize(3);
    for(int i = 0; i < 3; i++)
    {
        inter_node.p[i].resize(3);
    }
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.p[0][i] = inter_node.local.normal[i];
    }
    
    inter_node.p[1][0] = a[2]/sqrt(a[0]*a[0]+a[2]*a[2]);
    inter_node.p[1][1] = 0;
    inter_node.p[1][2] = -a[0]/sqrt(a[0]*a[0]+a[2]*a[2]);
    
    inter_node.p[2][0] = a[1]/sqrt(a[0]*a[0]+a[1]*a[1]);
    inter_node.p[2][1] = -a[0]/sqrt(a[0]*a[0]+a[1]*a[1]);
    inter_node.p[2][2] = 0;
}

/************************************************************************************************************************
 Normalized local coordinate transformation matrix initialization at given intersection node on xz plane in y-direction
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.p : local coordinate transformation matrix for given intersection node
 ************************************************************************************************************************/
void Intersections::Pmatrix_Setup_y(Intersection_Data& inter_node)
{
    VecDoub a;
    
    for(int i = 0; i < 3; i++)
    {
        a.push_back(inter_node.local.normal[i]);
    }
    
    inter_node.p.resize(3);
    for(int i = 0; i < 3; i++)
    {
        inter_node.p[i].resize(3);
    }
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.p[0][i] = a[i];
    }
    
    inter_node.p[1][0] = 0;
    inter_node.p[1][1] = a[2]/sqrt(a[1]*a[1]+a[2]*a[2]);
    inter_node.p[1][2] = -a[1]/sqrt(a[1]*a[1]+a[2]*a[2]);
    
    inter_node.p[2][0] = a[1]/sqrt(a[0]*a[0]+a[1]*a[1]);
    inter_node.p[2][1] = -a[0]/sqrt(a[0]*a[0]+a[1]*a[1]);
    inter_node.p[2][2] = 0;
}

/************************************************************************************************************************
 Normalized local coordinate transformation matrix initialization at given intersection node on xy plane in z-direction
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.p : local coordinate transformation matrix for given intersection node
 ************************************************************************************************************************/
void Intersections::Pmatrix_Setup_z(Intersection_Data& inter_node)
{
    VecDoub a;
    
    for(int i = 0; i < 3; i++)
    {
        a.push_back(inter_node.local.normal[i]);
    }
    
    inter_node.p.resize(3);
    for(int i = 0; i < 3; i++)
    {
        inter_node.p[i].resize(3);
    }
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.p[0][i] = a[i];
    }
    
    inter_node.p[1][0] = 0;
    inter_node.p[1][1] = a[2]/sqrt(a[1]*a[1]+a[2]*a[2]);
    inter_node.p[1][2] = -a[1]/sqrt(a[1]*a[1]+a[2]*a[2]);
    
    inter_node.p[2][0] = a[2]/sqrt(a[0]*a[0]+a[2]*a[2]);
    inter_node.p[2][1] = 0;
    inter_node.p[2][2] = -a[0]/sqrt(a[0]*a[0]+a[2]*a[2]);
}

/******************************************************************************************************************
 Finding best three nodes to approximate auxilary point on given direction
 
 INPUT
 axis : direction(x, y or z)
 ix   : left node x direction indx
 iy   : left node y direction indx
 iz   : left node z direction indx
 
 OUTPUT
 outside : indices of three nodes found on outside subdomain
 inside  : indices of three nodes found on inside subdomain
 ******************************************************************************************************************/
void Intersections::Search_indx(Char_I axis, Int_I ix, Int_I iy, Int_I iz, VecInt_O& outside, VecInt_O& inside)
{
    VecInt indx;
    int temp, dis_out, dis_in;
    
    indx.resize(3);
    
    //Initialize all indices of outside approximation on-grid nodes and inside approximation on-grid nodes with -1
    for(int i = 0; i < 3; i++)
    {
        outside[i] = -1;
        inside[i] = -1;
    }
    
    if(axis == 'x')
    {
        dis_out = nx*3;
        dis_in = nx*3;
        
        //i is the indx for tree search in given direction
        for(int i = 0; i < nx-2; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                indx[j] = To1d(i+j,iy,iz);
            }
            if(indicator[indx[0]]>0 && indicator[indx[1]]>0 && indicator[indx[2]]>0)
            {
                //relative distance
                temp = abs(i-ix+(i+1)-ix+(i+2)-ix);
                //if closer, replace indx
                if(temp < dis_out)
                {
                    dis_out = temp;
                    for(int j = 0; j < 3; j++)
                    {
                        outside[j] = i+j;
                    }
                }
            }
            if(indicator[indx[0]]<0 && indicator[indx[1]]<0 && indicator[indx[2]]<0)
            {
                temp = abs(i-ix+(i+1)-ix+(i+2)-ix);
                if(temp < dis_in)
                {
                    dis_in = temp;
                    for(int j = 0; j < 3; j++)
                    {
                        inside[j] = i+j;
                    }
                }
            }
        }
    }
    else if(axis == 'y')
    {
        dis_out = ny*3;
        dis_in = ny*3;
        
        //i is the indx for tree search in given direction
        for(int i = 0; i < nx-2; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                indx[j] = To1d(ix,i+j,iz);
            }
            if(indicator[indx[0]]>0 && indicator[indx[1]]>0 && indicator[indx[2]]>0)
            {
                //relative distance
                temp = abs(i-iy+(i+1)-iy+(i+2)-iy);
                //if closer, replace indx
                if(temp < dis_out)
                {
                    dis_out = temp;
                    for(int j = 0; j < 3; j++)
                    {
                        outside[j] = i+j;
                    }
                }
            }
            if(indicator[indx[0]]<0 && indicator[indx[1]]<0 && indicator[indx[2]]<0)
            {
                temp = abs(i-iy+(i+1)-iy+(i+2)-iy);
                if(temp < dis_in)
                {
                    dis_in = temp;
                    for(int j = 0; j < 3; j++)
                    {
                        inside[j] = i+j;
                    }
                }
            }
        }
    }
    else if(axis == 'z')
    {
        dis_out = nz*3;
        dis_in = nz*3;
        
        //i is the indx for tree search in given direction
        for(int i = 0; i < nz-2; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                indx[j] = To1d(ix,iy,i+j);
            }
            if(indicator[indx[0]]>0 && indicator[indx[1]]>0 && indicator[indx[2]]>0)
            {
                //relative distance
                temp = abs(i-iz+(i+1)-iz+(i+2)-iz);
                //if closer, replace indx
                if(temp < dis_out)
                {
                    dis_out = temp;
                    for(int j = 0; j < 3; j++)
                    {
                        outside[j] = i+j;
                    }
                }
            }
            if(indicator[indx[0]]<0 && indicator[indx[1]]<0 && indicator[indx[2]]<0)
            {
                temp = abs(i-iz+(i+1)-iz+(i+2)-iz);
                if(temp < dis_in)
                {
                    dis_in = temp;
                    for(int j = 0; j < 3; j++)
                    {
                        inside[j] = i+j;
                    }
                }
            }
        }
    }
    else
    {
        cout << "No axis is found" << endl;
        exit(0);
    }
}

/******************************************************************************************************************
 Initialize informations of auxiliary nodes on the Plane Y=IY with intersecting line of Y=IY Plane and Tangent Plane
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.auxlaxis       : the axis used for approximate auxiliary points (both upper and lower)
 inter_node.auxl           : coordinate location of auxiliary points
 inter_node.uin_auxlnodes  : indices of upper inside auxiliary nodes
 inter_node.uout_auxlnodes : indices of upper outside auxiliary nodes
 inter_node.lin_auxlnodes  : indices of lower inside auxiliary nodes
 inter_node.lout_auxlnodes : indices of lower outside auxiliary nodes
 inter_node.region         : approximation region (Omega^{+} or Omega^{-})
 ******************************************************************************************************************/
void Intersections::Auxiliary_eta_x(Intersection_Data& inter_node)
{
    bool flag_outside, flag_inside;
    int ix,iz,iy;                                                         //index of left node of intersection
    int upper, lower;                                                      //left index of upper and lower auxiliary points
    double outside_distance, inside_distance;                              //distance for comparing outside and inside
    double x0, z0;                                                         //x and z coordinate of intersections
    double a[3];                                                           //normalized normal direction for current intersection
    VecInt upper_outside, upper_inside, lower_outside, lower_inside;       //indices of each three nodes for auxiliary points on outside and inside
    
    //int x_upper, x_lower;
    int z_upper, z_lower;
    //double x_upper_auxl, x_lower_auxl, x_distance,;
    double z_upper_auxl, z_lower_auxl, z_distance;
    
    inter_node.eta.auxl.resize(2);
    //inter_node.eta.ul_indx.resize(2);
    inter_node.eta.auxlaxis.resize(2);
    
    upper_outside.resize(3);
    upper_inside.resize(3);
    lower_outside.resize(3);
    lower_inside.resize(3);
    
    inter_node.eta.uin_auxlnodes.resize(3);
    inter_node.eta.uout_auxlnodes.resize(3);
    inter_node.eta.lin_auxlnodes.resize(3);
    inter_node.eta.lout_auxlnodes.resize(3);
    
    ix = inter_node.left_loc;
    iy = inter_node.line.indx1;
    iz = inter_node.line.indx2;
    
    x0 = inter_node.coord.x_value;
    z0 = inter_node.coord.z_value;
    
    for(int i = 0; i < 3; i++)
    {
        a[i] = inter_node.local.normal[i];
    }
    
    /*
     for(int i = 0; i < nx; i++)
     {
     if((x0>xi[i]) && (x0<xi[i+1]))
     {
     x_upper = i+1;
     x_lower = i;
     }
     }
     */
    
    z_upper = iz+1;
    z_lower = iz-1;
    
    /*
     //Z coordinate of auxilary points on upper line X=x_upper and lower line X=x_lower
     x_upper_auxl = (a[0]*(x0-xi[x_upper])+a[2]*z0)/a[2];
     x_lower_auxl = (a[0]*(x0-xi[x_lower])+a[2]*z0)/a[2];
     //Distance of two auxilary points on two X lines
     x_distance = sqrt(dx*dx+abs(x_upper_auxl-x_lower_auxl));
     */
    
    //X coordinate of auxilary points on upper line Z=IZ+1 and lower line Z=IZ-1
    z_upper_auxl = (a[0]*x0+a[2]*(z0-zi[z_upper]))/a[0];
    z_lower_auxl = (a[0]*x0+a[2]*(z0-zi[z_lower]))/a[0];
    //Distance of two auxilary points on two Z lines
    //z_distance = sqrt(4*dz*dz+abs(z_upper_auxl-z_lower_auxl));
    
    //Choose direction which has a smaller distance
    //if(z_distance < x_distance)
    //{
    inter_node.eta.auxl[0] = z_upper_auxl;
    inter_node.eta.auxl[1] = z_lower_auxl;
    
    inter_node.eta.auxlaxis[0] = 'z';
    inter_node.eta.auxlaxis[1] = 'z';
    
    //inter_node.eta.ul_indx[0] = z_upper;
    //inter_node.eta.ul_indx[1] = z_lower;
    
    //Check if the auxilary points are in domain
    if((inter_node.eta.auxl[0]<xi[0]) || (inter_node.eta.auxl[0]>xi[nx-1]) || (inter_node.eta.auxl[1]<xi[0]) || (inter_node.eta.auxl[1]>xi[nx-1]))
    {
        cout << "Auxilary nodes in eta direction are out of boundary on x direction of Plane Y"<< endl;
        exit(0);
    }
    
    //Find index of left node for upper and lower auxiliary points on x-axis
    for(int i = 0; i < nx; i++)
    {
        if((xi[i] < inter_node.eta.auxl[0]) && (inter_node.eta.auxl[0] < xi[i+1]))
        {
            upper = i;
        }
        if((xi[i] < inter_node.eta.auxl[1]) && (inter_node.eta.auxl[1] < xi[i+1]))
        {
            lower = i;
        }
    }
    
    Search_indx('x',upper,iy,z_upper,upper_outside,upper_inside);
    Search_indx('x',lower,iy,z_lower,lower_outside,lower_inside);
    // }
    /*
     else
     {
     inter_node.eta.auxl[0] = x_upper_auxl;
     inter_node.eta.auxl[1] = x_lower_auxl;
     
     inter_node.eta.auxlaxis[0] = 'x';
     inter_node.eta.auxlaxis[1] = 'x';
     
     inter_node.eta.ul_indx[0] = x_upper;
     inter_node.eta.ul_indx[1] = x_lower;
     
     //Check if the auxilary points are in domain
     if((inter_node.eta.auxl[0]<zi[0]) || (inter_node.eta.auxl[0]>zi[nz-1]) || (inter_node.eta.auxl[1]<zi[0]) || (inter_node.eta.auxl[1]>zi[nz-1]))
     {
     cout << "Auxilary nodes in eta direction are out of boundary on z direction of Plane Y"<< endl;
     exit(0);
     }
     
     //Find index of left node for upper and lower auxiliary points on z-axis
     for(int i = 0; i < nx; i++)
     {
     if((zi[i] < inter_node.eta.auxl[0]) && (inter_node.eta.auxl[0] < zi[i+1]))
     {
     upper = i;
     }
     if((zi[i] < inter_node.eta.auxl[1]) && (inter_node.eta.auxl[1] < zi[i+1]))
     {
     lower = i;
     }
     }
     
     Search_indx('z',x_upper,iy,upper,upper_outside,upper_inside);
     Search_indx('z',x_lower,iy,lower,lower_outside,lower_inside);
     }
     */
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.eta.uin_auxlnodes[i] = upper_inside[i];
        inter_node.eta.uout_auxlnodes[i] = upper_outside[i];
        inter_node.eta.lin_auxlnodes[i] = lower_inside[i];
        inter_node.eta.lout_auxlnodes[i] = lower_outside[i];
    }
    
    //Initialize lower&upper distance with 0
    inside_distance = 0;
    outside_distance = 0;
    
    //Flags for four different cases when using outside or inside nodes for auxiliary points
    flag_inside = (upper_inside[0]>=0) && (lower_inside[0]>=0);
    flag_outside = (upper_outside[0]>=0) && (lower_outside[0]>=0);
    
    //Case 1, cannot find three nodes on same side for both upper and lower auxiliary points, approximation fails
    if(!flag_inside && !flag_outside)
    {
        cout << "Eta approximation fails" << endl;
        exit(0);
    }
    //Case 2, find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    else if(!flag_inside && flag_outside)
    {
        inter_node.eta.region = 'o';
    }
    //Case 3, find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(flag_inside && !flag_outside)
    {
        inter_node.eta.region = 'i';
    }
    //Case 4, both inside and outside domain available
    else
    {
        //Compare which six nodes is closer to two auxiliary points
        if(inter_node.eta.auxlaxis[0] == 'z')
        {
            for(int i = 0; i < 3; i++)
            {
                outside_distance += abs(xi[upper_outside[i]]-inter_node.eta.auxl[0])+abs(xi[lower_outside[i]]-inter_node.eta.auxl[1]);
                inside_distance += abs(xi[upper_inside[i]]-inter_node.eta.auxl[0])+abs(xi[lower_inside[i]]-inter_node.eta.auxl[1]);
            }
        }
        else if(inter_node.eta.auxlaxis[0] == 'x')
        {
            for(int i = 0; i < 3; i++)
            {
                outside_distance += abs(zi[upper_outside[i]]-inter_node.eta.auxl[0])+abs(zi[lower_outside[i]]-inter_node.eta.auxl[1]);
                inside_distance += abs(zi[upper_inside[i]]-inter_node.eta.auxl[0])+abs(zi[lower_inside[i]]-inter_node.eta.auxl[1]);
            }
        }
        
        //Find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
        if(outside_distance < inside_distance)
        {
            inter_node.eta.region = 'o';
        }
        //Find six nodes available on inside Omega^{-} domain, use all inside nodes for auxiliary points approximation
        else
        {
            inter_node.eta.region = 'i';
        }
    }
}

/******************************************************************************************************************
 Initialize informations of auxiliary nodes on the Plane Z=IZ with intersecting line of Z=IZ Plane and Tangent Plane
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.auxlaxis       : the axis used for approximate auxiliary points (both upper and lower)
 inter_node.auxl           : coordinate location of auxiliary points
 inter_node.uin_auxlnodes  : indices of upper inside auxiliary nodes
 inter_node.uout_auxlnodes : indices of upper outside auxiliary nodes
 inter_node.lin_auxlnodes  : indices of lower inside auxiliary nodes
 inter_node.lout_auxlnodes : indices of lower outside auxiliary nodes
 inter_node.region         : approximation region (Omega^{+} or Omega^{-})
 ******************************************************************************************************************/
void Intersections::Auxiliary_tau_x(Intersection_Data& inter_node)
{
    bool flag_outside, flag_inside;
    int ix,iz,iy;                                               //index of left node of intersection
    int ix_upper, ix_lower;                                     //left index of upper and lower auxiliary points
    double outside_distance, inside_distance;                   //distance for comparing outside and inside
    double x0, y0;                                              //x and z coordinate of intersections; u_{tau}^{+} or u_{tau}^{-}
    double a[3];                                                //normalized normal direction for current intersection
    VecInt upper_outside, upper_inside, lower_outside, lower_inside;       //indices of each three nodes for auxiliary points on outside and inside
    
    upper_outside.resize(3);
    upper_inside.resize(3);
    lower_outside.resize(3);
    lower_inside.resize(3);
    
    inter_node.tau.uin_auxlnodes.resize(3);
    inter_node.tau.uout_auxlnodes.resize(3);
    inter_node.tau.lin_auxlnodes.resize(3);
    inter_node.tau.lout_auxlnodes.resize(3);
    
    ix = inter_node.left_loc;
    iy = inter_node.line.indx1;
    iz = inter_node.line.indx2;
    
    x0 = inter_node.coord.x_value;
    y0 = inter_node.coord.y_value;
    
    for(int i = 0; i < 3; i++)
    {
        a[i] = inter_node.local.normal[i];
    }
    
    inter_node.tau.auxlaxis.resize(2);
    inter_node.tau.auxlaxis[0] = 'y';
    inter_node.tau.auxlaxis[1] = 'y';
    
    //X coordinate of auxilary points on upper line Y=IY+1 and lower line Y=IY-1
    inter_node.tau.auxl.resize(2);
    inter_node.tau.auxl[0] = (a[0]*x0+a[1]*(y0-yi[iy+1]))/a[0];
    inter_node.tau.auxl[1] = (a[0]*x0+a[1]*(y0-yi[iy-1]))/a[0];
    
    //Check if the auxilary points are in domain
    if((inter_node.tau.auxl[0]<xi[0]) || (inter_node.tau.auxl[0]>xi[nx-1]) || (inter_node.tau.auxl[1]<xi[0]) || (inter_node.tau.auxl[1]>xi[nx-1]))
    {
        cout << "Auxilary nodes in tau direction are out of boundary on x direction of Plane Z"<< endl;
        exit(0);
    }
    
    //Find index of left node for upper and lower auxiliary points
    for(int i = 0; i < nx; i++)
    {
        if((xi[i] < inter_node.tau.auxl[0]) && (inter_node.tau.auxl[0] < xi[i+1]))
        {
            ix_upper = i;
        }
        if((xi[i] < inter_node.tau.auxl[1]) && (inter_node.tau.auxl[1]< xi[i+1]))
        {
            ix_lower = i;
        }
    }
    
    Search_indx('x',ix_upper,iy+1,iz,upper_outside,upper_inside);
    Search_indx('x',ix_lower,iy-1,iz,lower_outside,lower_inside);
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.tau.uin_auxlnodes[i] = upper_inside[i];
        inter_node.tau.uout_auxlnodes[i] = upper_outside[i];
        inter_node.tau.lin_auxlnodes[i] = lower_inside[i];
        inter_node.tau.lout_auxlnodes[i] = lower_outside[i];
    }
    
    //Initialize lower&upper distance with 0
    inside_distance = 0;
    outside_distance = 0;
    
    //Flags for four different cases when using outside or inside nodes for auxiliary points
    flag_inside = (upper_inside[0]>=0) && (lower_inside[0]>=0);
    flag_outside = (upper_outside[0]>=0) && (lower_outside[0]>=0);
    
    //Case 1, cannot find three nodes on same side for both upper and lower auxiliary points, approximation fails
    if(!flag_inside && !flag_outside)
    {
        cout << "Tau approximation fails" << endl;
        exit(0);
    }
    //Case 2, find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    else if(!flag_inside && flag_outside)
    {
        inter_node.tau.region = 'o';
    }
    //Case 3, find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(flag_inside && !flag_outside)
    {
        inter_node.tau.region = 'i';
    }
    //Case 4, both inside and outside domain available
    else
    {
        //Compare which six nodes is closer to two auxiliary points
        for(int i = 0; i < 3; i++)
        {
            outside_distance += abs(xi[upper_outside[i]]-inter_node.tau.auxl[0])+abs(xi[lower_outside[i]]-inter_node.tau.auxl[1]);
            inside_distance += abs(xi[upper_inside[i]]-inter_node.tau.auxl[0])+abs(xi[lower_inside[i]]-inter_node.tau.auxl[1]);
        }
        
        //Find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
        if(outside_distance < inside_distance)
        {
            inter_node.tau.region = 'o';
        }
        //Find six nodes available on inside Omega^{-} domain, use all inside nodes for auxiliary points approximation
        else
        {
            inter_node.tau.region = 'i';
        }
    }
}

/******************************************************************************************************************
 Initialize informations of auxiliary nodes on the Plane X=IX with intersecting line of X=IX Plane and Tangent Plane
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.auxlaxis       : the axis used for approximate auxiliary points (both upper and lower)
 inter_node.auxl           : coordinate location of auxiliary points
 inter_node.uin_auxlnodes  : indices of upper inside auxiliary nodes
 inter_node.uout_auxlnodes : indices of upper outside auxiliary nodes
 inter_node.lin_auxlnodes  : indices of lower inside auxiliary nodes
 inter_node.lout_auxlnodes : indices of lower outside auxiliary nodes
 inter_node.region         : approximation region (Omega^{+} or Omega^{-})
 ******************************************************************************************************************/
void Intersections::Auxiliary_eta_y(Intersection_Data& inter_node)
{
    bool flag_outside, flag_inside;
    int ix,iz,iy;                                               //index of left node of intersection
    int iy_upper, iy_lower;                                     //left index of upper and lower auxiliary points
    double outside_distance, inside_distance;                   //distance for comparing outside and inside
    double y0, z0;                                              //y and z coordinate of intersections; u_{eta}^{+} or u_{eta}^{-}
    double a[3];                                                //normalized normal direction for current intersection
    VecInt upper_outside, upper_inside, lower_outside, lower_inside;       //indices of each three nodes for auxiliary points on outside and inside
    
    upper_outside.resize(3);
    upper_inside.resize(3);
    lower_outside.resize(3);
    lower_inside.resize(3);
    
    inter_node.eta.uin_auxlnodes.resize(3);
    inter_node.eta.uout_auxlnodes.resize(3);
    inter_node.eta.lin_auxlnodes.resize(3);
    inter_node.eta.lout_auxlnodes.resize(3);
    
    ix = inter_node.line.indx1;
    iy = inter_node.left_loc;
    iz = inter_node.line.indx2;
    
    y0 = inter_node.coord.y_value;
    z0 = inter_node.coord.z_value;
    
    for(int i = 0; i < 3; i++)
    {
        a[i] = inter_node.local.normal[i];
    }
    
    inter_node.eta.auxlaxis.resize(2);
    inter_node.eta.auxlaxis[0] = 'z';
    inter_node.eta.auxlaxis[1] = 'z';
    
    //Y coordinate of auxilary points on upper line Z=IZ+1 and lower line Z=IZ-1
    inter_node.eta.auxl.resize(2);
    inter_node.eta.auxl[0] = (a[1]*y0+a[2]*(z0-zi[iz+1]))/a[1];
    inter_node.eta.auxl[1] = (a[1]*y0+a[2]*(z0-zi[iz-1]))/a[1];
    
    //Check if the auxilary points are in domain
    if((inter_node.eta.auxl[0]<yi[0]) || (inter_node.eta.auxl[0]>yi[ny-1]) || (inter_node.eta.auxl[1]<yi[0]) || (inter_node.eta.auxl[1]>yi[ny-1]))
    {
        cout << "Auxilary nodes in eta direction are out of boundary on y direction of Plane X"<< endl;
        exit(0);
    }
    
    //Find index of left node for upper and lower auxiliary points
    for(int i = 0; i < ny; i++)
    {
        if((yi[i] < inter_node.eta.auxl[0]) && (inter_node.eta.auxl[0] < yi[i+1]))
        {
            iy_upper = i;
        }
        if((yi[i] < inter_node.eta.auxl[1]) && (inter_node.eta.auxl[1] < yi[i+1]))
        {
            iy_lower = i;
        }
    }
    
    Search_indx('y',ix,iy_upper,iz+1,upper_outside,upper_inside);
    Search_indx('y',ix,iy_lower,iz-1,lower_outside,lower_inside);
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.eta.uin_auxlnodes[i] = upper_inside[i];
        inter_node.eta.uout_auxlnodes[i] = upper_outside[i];
        inter_node.eta.lin_auxlnodes[i] = lower_inside[i];
        inter_node.eta.lout_auxlnodes[i] = lower_outside[i];
    }
    
    //Initialize lower&upper distance with 0
    inside_distance = 0;
    outside_distance = 0;
    
    //Flags for four different cases when using outside or inside nodes for auxiliary points
    flag_inside = (upper_inside[0]>=0) && (lower_inside[0]>=0);
    flag_outside = (upper_outside[0]>=0) && (lower_outside[0]>=0);
    
    //Case 1, cannot find three nodes on same side for both upper and lower auxiliary points, approximation fails
    if(!flag_inside && !flag_outside)
    {
        cout << "Eta approximation fails" << endl;
        exit(0);
    }
    //Case 2, find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    else if(!flag_inside && flag_outside)
    {
        inter_node.eta.region = 'o';
    }
    //Case 3, find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(flag_inside && !flag_outside)
    {
        inter_node.eta.region = 'i';
    }
    //Case 4, both inside and outside domain available
    else
    {
        //Compare which six nodes is closer to two auxiliary points
        for(int i = 0; i < 3; i++)
        {
            outside_distance += abs(yi[upper_outside[i]]-inter_node.eta.auxl[0])+abs(yi[lower_outside[i]]-inter_node.eta.auxl[1]);
            inside_distance += abs(yi[upper_inside[i]]-inter_node.eta.auxl[0])+abs(yi[lower_inside[i]]-inter_node.eta.auxl[1]);
        }
        
        //Find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
        if(outside_distance < inside_distance)
        {
            inter_node.eta.region = 'o';
        }
        //Find six nodes available on inside Omega^{-} domain, use all inside nodes for auxiliary points approximation
        else
        {
            inter_node.eta.region = 'i';
        }
    }
}

/******************************************************************************************************************
 Initialize informations of auxiliary nodes on the Plane Z=IZ with intersecting line of Z=IZ Plane and Tangent Plane
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.auxlaxis       : the axis used for approximate auxiliary points (both upper and lower)
 inter_node.auxl           : coordinate location of auxiliary points
 inter_node.uin_auxlnodes  : indices of upper inside auxiliary nodes
 inter_node.uout_auxlnodes : indices of upper outside auxiliary nodes
 inter_node.lin_auxlnodes  : indices of lower inside auxiliary nodes
 inter_node.lout_auxlnodes : indices of lower outside auxiliary nodes
 inter_node.region         : approximation region (Omega^{+} or Omega^{-})
 ******************************************************************************************************************/
void Intersections::Auxiliary_tau_y(Intersection_Data& inter_node)
{
    bool flag_outside, flag_inside;
    int ix,iz,iy;                                               //index of left node of intersection
    int iy_upper, iy_lower;                                     //left index of upper and lower auxiliary points
    double outside_distance, inside_distance;                   //distance for comparing outside and inside
    double y0, x0;                                              //y and x coordinate of intersections; u_{tau}^{+} or u_{tau}^{-}
    double a[3];                                                //normalized normal direction for current intersection
    VecInt upper_outside, upper_inside, lower_outside, lower_inside;       //indices of each three nodes for auxiliary points on outside and inside
    
    upper_outside.resize(3);
    upper_inside.resize(3);
    lower_outside.resize(3);
    lower_inside.resize(3);
    
    inter_node.tau.uin_auxlnodes.resize(3);
    inter_node.tau.uout_auxlnodes.resize(3);
    inter_node.tau.lin_auxlnodes.resize(3);
    inter_node.tau.lout_auxlnodes.resize(3);
    
    ix = inter_node.line.indx1;
    iy = inter_node.left_loc;
    iz = inter_node.line.indx2;
    
    y0 = inter_node.coord.y_value;
    x0 = inter_node.coord.x_value;
    
    for(int i = 0; i < 3; i++)
    {
        a[i] = inter_node.local.normal[i];
    }
    
    inter_node.tau.auxlaxis.resize(2);
    inter_node.tau.auxlaxis[0] = 'x';
    inter_node.tau.auxlaxis[1] = 'x';
    
    //Y coordinate of auxilary points on upper line X=IX+1 and lower line X=IX-1
    inter_node.tau.auxl.resize(2);
    inter_node.tau.auxl[0] = (a[1]*y0+a[0]*(x0-xi[ix+1]))/a[1];
    inter_node.tau.auxl[1] = (a[1]*y0+a[0]*(x0-xi[ix-1]))/a[1];
    
    //Check if the auxilary points are in domain
    if((inter_node.tau.auxl[0]<yi[0]) || (inter_node.tau.auxl[0]>yi[ny-1]) || (inter_node.tau.auxl[1]<yi[0]) || (inter_node.tau.auxl[1]>yi[ny-1]))
    {
        cout << "Auxilary nodes in tau direction are out of boundary on y direction of Plane Z"<< endl;
        exit(0);
    }
    
    //Find index of left node for upper and lower auxiliary points
    for(int i = 0; i < ny; i++)
    {
        if((yi[i] < inter_node.tau.auxl[0]) && (inter_node.tau.auxl[0] < yi[i+1]))
        {
            iy_upper = i;
        }
        if((yi[i] < inter_node.tau.auxl[1]) && (inter_node.tau.auxl[1] < yi[i+1]))
        {
            iy_lower = i;
        }
    }
    
    Search_indx('y',ix+1,iy_upper,iz,upper_outside,upper_inside);
    Search_indx('y',ix-1,iy_lower,iz,lower_outside,lower_inside);
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.tau.uin_auxlnodes[i] = upper_inside[i];
        inter_node.tau.uout_auxlnodes[i] = upper_outside[i];
        inter_node.tau.lin_auxlnodes[i] = lower_inside[i];
        inter_node.tau.lout_auxlnodes[i] = lower_outside[i];
    }
    
    //Initialize lower&upper distance with 0
    inside_distance = 0;
    outside_distance = 0;
    
    //Flags for four different cases when using outside or inside nodes for auxiliary points
    flag_inside = (upper_inside[0]>=0) && (lower_inside[0]>=0);
    flag_outside = (upper_outside[0]>=0) && (lower_outside[0]>=0);
    
    //Case 1, cannot find three nodes on same side for both upper and lower auxiliary points, approximation fails
    if(!flag_inside && !flag_outside)
    {
        cout << "Tau approximation fails" << endl;
        exit(0);
    }
    //Case 2, find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    else if(!flag_inside && flag_outside)
    {
        inter_node.tau.region = 'o';
    }
    //Case 3, find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(flag_inside && !flag_outside)
    {
        inter_node.tau.region = 'i';
    }
    //Case 4, both inside and outside domain available
    else
    {
        //Compare which six nodes is closer to two auxiliary points
        for(int i = 0; i < 3; i++)
        {
            outside_distance += abs(yi[upper_outside[i]]-inter_node.tau.auxl[0])+abs(yi[lower_outside[i]]-inter_node.tau.auxl[1]);
            inside_distance += abs(yi[upper_inside[i]]-inter_node.tau.auxl[0])+abs(yi[lower_inside[i]]-inter_node.tau.auxl[1]);
        }
        
        //Find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
        if(outside_distance < inside_distance)
        {
            inter_node.tau.region = 'o';
        }
        //Find six nodes available on inside Omega^{-} domain, use all inside nodes for auxiliary points approximation
        else
        {
            inter_node.tau.region = 'i';
        }
    }
}

/******************************************************************************************************************
 Initialize informations of auxiliary nodes on the Plane X=IX with intersecting line of X=IX Plane and Tangent Plane
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.auxlaxis       : the axis used for approximate auxiliary points (both upper and lower)
 inter_node.auxl           : coordinate location of auxiliary points
 inter_node.uin_auxlnodes  : indices of upper inside auxiliary nodes
 inter_node.uout_auxlnodes : indices of upper outside auxiliary nodes
 inter_node.lin_auxlnodes  : indices of lower inside auxiliary nodes
 inter_node.lout_auxlnodes : indices of lower outside auxiliary nodes
 inter_node.region         : approximation region (Omega^{+} or Omega^{-})
 ******************************************************************************************************************/
void Intersections::Auxiliary_eta_z(Intersection_Data& inter_node)
{
    bool flag_outside, flag_inside;
    int ix,iz,iy;                                               //index of left node of intersection
    int iz_upper, iz_lower;                                     //left index of upper and lower auxiliary points
    double outside_distance, inside_distance;                   //distance for comparing outside and inside
    double y0, z0;                                              //y and z coordinate of intersections; u_{eta}^{+} or u_{eta}^{-}
    double a[3];                                                //normalized normal direction for current intersection
    VecInt upper_outside, upper_inside, lower_outside, lower_inside;       //indices of each three nodes for auxiliary points on outside and inside
    
    upper_outside.resize(3);
    upper_inside.resize(3);
    lower_outside.resize(3);
    lower_inside.resize(3);
    
    inter_node.eta.uin_auxlnodes.resize(3);
    inter_node.eta.uout_auxlnodes.resize(3);
    inter_node.eta.lin_auxlnodes.resize(3);
    inter_node.eta.lout_auxlnodes.resize(3);
    
    ix = inter_node.line.indx1;
    iy = inter_node.line.indx2;
    iz = inter_node.left_loc;
    
    y0 = inter_node.coord.y_value;
    z0 = inter_node.coord.z_value;
    
    for(int i = 0; i < 3; i++)
    {
        a[i] = inter_node.local.normal[i];
    }
    
    inter_node.eta.auxlaxis.resize(2);
    inter_node.eta.auxlaxis[0] = 'y';
    inter_node.eta.auxlaxis[1] = 'y';
    
    //Z coordinate of auxilary points on upper line Y=IY+1 and lower line Y=IY-1
    inter_node.eta.auxl.resize(2);
    inter_node.eta.auxl[0] = (a[2]*z0+a[1]*(y0-yi[iy+1]))/a[2];
    inter_node.eta.auxl[1] = (a[2]*z0+a[1]*(y0-yi[iy-1]))/a[2];
    
    //Check if the auxilary points are in domain
    if((inter_node.eta.auxl[0]<zi[0]) || (inter_node.eta.auxl[0]>zi[nz-1]) || (inter_node.eta.auxl[1]<zi[0]) || (inter_node.eta.auxl[1]>zi[nz-1]))
    {
        cout << "Auxilary nodes in eta direction are out of boundary on z direction of Plane X"<< endl;
        exit(0);
    }
    
    //Find index of left node for upper and lower auxiliary points
    for(int i = 0; i < nz; i++)
    {
        if((zi[i] < inter_node.eta.auxl[0]) && (inter_node.eta.auxl[0] < zi[i+1]))
        {
            iz_upper = i;
        }
        if((zi[i] < inter_node.eta.auxl[1]) && (inter_node.eta.auxl[1] < zi[i+1]))
        {
            iz_lower = i;
        }
    }
    
    Search_indx('z',ix,iy+1,iz_upper,upper_outside,upper_inside);
    Search_indx('z',ix,iy-1,iz_lower,lower_outside,lower_inside);
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.eta.uin_auxlnodes[i] = upper_inside[i];
        inter_node.eta.uout_auxlnodes[i] = upper_outside[i];
        inter_node.eta.lin_auxlnodes[i] = lower_inside[i];
        inter_node.eta.lout_auxlnodes[i] = lower_outside[i];
    }
    
    //Initialize lower&upper distance with 0
    inside_distance = 0;
    outside_distance = 0;
    
    //Flags for four different cases when using outside or inside nodes for auxiliary points
    flag_inside = (upper_inside[0]>=0) && (lower_inside[0]>=0);
    flag_outside = (upper_outside[0]>=0) && (lower_outside[0]>=0);
    
    //Case 1, cannot find three nodes on same side for both upper and lower auxiliary points, approximation fails
    if(!flag_inside && !flag_outside)
    {
        cout << "Eta approximation fails" << endl;
        exit(0);
    }
    //Case 2, find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    else if(!flag_inside && flag_outside)
    {
        inter_node.eta.region = 'o';
    }
    //Case 3, find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(flag_inside && !flag_outside)
    {
        inter_node.eta.region = 'i';
    }
    //Case 4, both inside and outside domain available
    else
    {
        //Compare which six nodes is closer to two auxiliary points
        for(int i = 0; i < 3; i++)
        {
            outside_distance += abs(zi[upper_outside[i]]-inter_node.eta.auxl[0])+abs(zi[lower_outside[i]]-inter_node.eta.auxl[1]);
            inside_distance += abs(zi[upper_inside[i]]-inter_node.eta.auxl[0])+abs(zi[lower_inside[i]]-inter_node.eta.auxl[1]);
        }
        
        //Find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
        if(outside_distance < inside_distance)
        {
            inter_node.eta.region = 'o';
        }
        //Find six nodes available on inside Omega^{-} domain, use all inside nodes for auxiliary points approximation
        else
        {
            inter_node.eta.region = 'i';
        }
    }
}

/******************************************************************************************************************
 Initialize informations of auxiliary nodes on the Plane Y=IY with intersecting line of Y=IY Plane and Tangent Plane
 
 INPUT
 inter_node : given intersection node
 
 OUTPUT
 inter_node.auxlaxis       : the axis used for approximate auxiliary points (both upper and lower)
 inter_node.auxl           : coordinate location of auxiliary points
 inter_node.uin_auxlnodes  : indices of upper inside auxiliary nodes
 inter_node.uout_auxlnodes : indices of upper outside auxiliary nodes
 inter_node.lin_auxlnodes  : indices of lower inside auxiliary nodes
 inter_node.lout_auxlnodes : indices of lower outside auxiliary nodes
 inter_node.region         : approximation region (Omega^{+} or Omega^{-})
 ******************************************************************************************************************/
void Intersections::Auxiliary_tau_z(Intersection_Data& inter_node)
{
    bool flag_outside, flag_inside;
    int ix,iz,iy;                                               //index of left node of intersection
    int iz_upper, iz_lower;                                     //left index of upper and lower auxiliary points
    double outside_distance, inside_distance;                   //distance for comparing outside and inside
    double x0, z0;                                              //x and z coordinate of intersections; u_{tau}^{+} or u_{tau}^{-}
    double a[3];                                                //normalized normal direction for current intersection
    VecInt upper_outside, upper_inside, lower_outside, lower_inside;       //indices of each three nodes for auxiliary points on outside and inside
    
    upper_outside.resize(3);
    upper_inside.resize(3);
    lower_outside.resize(3);
    lower_inside.resize(3);
    
    inter_node.tau.uin_auxlnodes.resize(3);
    inter_node.tau.uout_auxlnodes.resize(3);
    inter_node.tau.lin_auxlnodes.resize(3);
    inter_node.tau.lout_auxlnodes.resize(3);
    
    ix = inter_node.line.indx1;
    iy = inter_node.line.indx2;
    iz = inter_node.left_loc;
    
    x0 = inter_node.coord.x_value;
    z0 = inter_node.coord.z_value;
    
    for(int i = 0; i < 3; i++)
    {
        a[i] = inter_node.local.normal[i];
    }
    
    inter_node.tau.auxlaxis.resize(2);
    inter_node.tau.auxlaxis[0] = 'x';
    inter_node.tau.auxlaxis[1] = 'x';
    
    //Z coordinate of auxilary points on upper line X=IX+1 and lower line X=IX-1
    inter_node.tau.auxl.resize(2);
    inter_node.tau.auxl[0] = (a[2]*z0+a[0]*(x0-xi[ix+1]))/a[2];
    inter_node.tau.auxl[1] = (a[2]*z0+a[0]*(x0-xi[ix-1]))/a[2];
    
    //Check if the auxilary points are in domain
    if((inter_node.tau.auxl[0]<zi[0]) || (inter_node.tau.auxl[0]>zi[nz-1]) || (inter_node.tau.auxl[1]<zi[0]) || (inter_node.tau.auxl[1]>zi[nz-1]))
    {
        cout << "Auxilary nodes in tau direction are out of boundary on z direction of Plane Y"<< endl;
        exit(0);
    }
    
    //Find index of left node for upper and lower auxiliary points
    for(int i = 0; i < nz; i++)
    {
        if((zi[i] < inter_node.tau.auxl[0]) && (inter_node.tau.auxl[0] < zi[i+1]))
        {
            iz_upper = i;
        }
        if((zi[i] < inter_node.tau.auxl[1]) && (inter_node.tau.auxl[1] < zi[i+1]))
        {
            iz_lower = i;
        }
    }
    
    Search_indx('z',ix+1,iy,iz_upper,upper_outside,upper_inside);
    Search_indx('z',ix-1,iy,iz_lower,lower_outside,lower_inside);
    
    for(int i = 0; i < 3; i++)
    {
        inter_node.tau.uin_auxlnodes[i] = upper_inside[i];
        inter_node.tau.uout_auxlnodes[i] = upper_outside[i];
        inter_node.tau.lin_auxlnodes[i] = lower_inside[i];
        inter_node.tau.lout_auxlnodes[i] = lower_outside[i];
    }
    
    //Initialize lower&upper distance with 0
    inside_distance = 0;
    outside_distance = 0;
    
    //Flags for four different cases when using outside or inside nodes for auxiliary points
    flag_inside = (upper_inside[0]>=0) && (lower_inside[0]>=0);
    flag_outside = (upper_outside[0]>=0) && (lower_outside[0]>=0);
    
    //Case 1, cannot find three nodes on same side for both upper and lower auxiliary points, approximation fails
    if(!flag_inside && !flag_outside)
    {
        cout << "Tau approximation fails" << endl;
        exit(0);
    }
    //Case 2, find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    else if(!flag_inside && flag_outside)
    {
        inter_node.tau.region = 'o';
    }
    //Case 3, find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(flag_inside && !flag_outside)
    {
        inter_node.tau.region = 'i';
    }
    //Case 4, both inside and outside domain available
    else
    {
        //Compare which six nodes is closer to two auxiliary points
        for(int i = 0; i < 3; i++)
        {
            outside_distance += abs(zi[upper_outside[i]]-inter_node.tau.auxl[0])+abs(zi[lower_outside[i]]-inter_node.tau.auxl[1]);
            inside_distance += abs(zi[upper_inside[i]]-inter_node.tau.auxl[0])+abs(zi[lower_inside[i]]-inter_node.tau.auxl[1]);
        }
        
        //Find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
        if(outside_distance < inside_distance)
        {
            inter_node.tau.region = 'o';
        }
        //Find six nodes available on inside Omega^{-} domain, use all inside nodes for auxiliary points approximation
        else
        {
            inter_node.tau.region = 'i';
        }
    }
}

/*******************************************************************************************
 Calculating of jump beta Eta for given intersection on x-direction based on current solution
 
 INPUT
 inter_node : given intersection node
 uh         : current solution
 beta       : beta object
 
 OUTPUT
 jumpbeta_eta : jump beta eta for given intersection based on current solution
 *******************************************************************************************/
double Intersections::Eta_x(Intersection_Data& inter_node, CubicDoub& uh, Beta& beta)
{
    int ix, iy, iz;                                                 //index of left node of intersection
    double upper_distance, lower_distance;                          //distance for eta approximation; distance for comparing outside and inside
    double x_upper, x_lower, upper_value, lower_value;              //x coordinate and value of two auxiliary points
    double x0, u_eta;                                               //u_{eta}^{+} or u_{eta}^{-}
    VecDoub upper, lower, eta;                                      //absolute location for approximating two auxiliary points and u_{eta}
    MatrixDoub wei_upper, wei_lower, wei_eta;                       //weights for two auxiliary points and weights for u_{eta}
    double jumpbeta_eta;
    
    
    ix = inter_node.left_loc;
    iy = inter_node.line.indx1;
    iz = inter_node.line.indx2;
    
    x0 = inter_node.coord.x_value;
    
    x_upper = inter_node.eta.auxl[0];
    x_lower = inter_node.eta.auxl[1];
    
    //Initialize two auxiliary points with 0
    upper_value = 0;
    lower_value = 0;
    
    upper_distance = sqrt(abs(x_upper-x0)*abs(x_upper-x0)+dz*dz);
    lower_distance = sqrt(abs(x_lower-x0)*abs(x_lower-x0)+dz*dz);
    
    /*
     upper_distance = sqrt(abs(x_upper-x0)*abs(x_upper-x0)*inter_node.p[1][0]*inter_node.p[1][0]+
     dz*dz*inter_node.p[1][2]*inter_node.p[1][2]);
     lower_distance = sqrt(abs(x_lower-x0)*abs(x_lower-x0)*inter_node.p[1][0]*inter_node.p[1][0]+
     dz*dz*inter_node.p[1][2]*inter_node.p[1][2]);
     */
    
    //Get Eta derivative along z-direction, make sure it's the same direction of given vector
    if(inter_node.p[1][2] > 0)
    {
        eta.push_back(-lower_distance);
        eta.push_back(0.0);
        eta.push_back(upper_distance);
    }
    else if(inter_node.p[1][2] < 0)
    {
        eta.push_back(lower_distance);
        eta.push_back(0.0);
        eta.push_back(-upper_distance);
    }
    else
    {
        cout << "x0 is on x-axis" << endl;
        exit(0);
    }
    
    wei_eta.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_eta[i].resize(2);
    }
    Weights(0,eta,3,1,wei_eta);
    
    //Weights for upper and lower auxiliary points approximation
    wei_upper.resize(3);
    wei_lower.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_upper[i].resize(1);
        wei_lower[i].resize(1);
    }
    
    //find six nodes available on outside Omega^{+} domain only, use all outside nodes for auxiliary points approximation
    if((REG == 'o') || (inter_node.eta.region == 'o'))
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(xi[inter_node.eta.uout_auxlnodes[i]]);
            lower.push_back(xi[inter_node.eta.lout_auxlnodes[i]]);
        }
        
        Weights(x_upper,upper,3,0,wei_upper);
        Weights(x_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[inter_node.eta.uout_auxlnodes[i]][iy][iz+1]*wei_upper[i][0];
            lower_value += uh[inter_node.eta.lout_auxlnodes[i]][iy][iz-1]*wei_lower[i][0];
        }
        
        u_eta = wei_eta[0][1]*lower_value+wei_eta[2][1]*upper_value;
        
        jumpbeta_eta = Eta_Plus(inter_node,beta,u_eta);
    }
    //find six nodes available on inside Omega^{-} domain only, use all inside nodes for auxiliary points approximation
    else if(inter_node.eta.region == 'i')
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(xi[inter_node.eta.uin_auxlnodes[i]]);
            lower.push_back(xi[inter_node.eta.lin_auxlnodes[i]]);
        }
        
        Weights(x_upper,upper,3,0,wei_upper);
        Weights(x_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[inter_node.eta.uin_auxlnodes[i]][iy][iz+1]*wei_upper[i][0];
            lower_value += uh[inter_node.eta.lin_auxlnodes[i]][iy][iz-1]*wei_lower[i][0];
        }
        
        u_eta = wei_eta[0][1]*lower_value+wei_eta[2][1]*upper_value;
        
        jumpbeta_eta = Eta_Minus(inter_node,beta,u_eta);
    }
    else
    {
        cout << "Approximation region isn't found, should only be 'o' for outside or 'i' for inside" << endl;
        exit(0);
    }
    
    return jumpbeta_eta;
}

/*******************************************************************************************
 Calculating of jump beta Tau for given intersection on x-direction based on current solution
 
 INPUT
 inter_node : given intersection node
 uh         : current solution
 beta       : beta object
 
 OUTPUT
 jumpbeta_tau : jump beta tau for given intersection based on current solution
 *******************************************************************************************/
double Intersections::Tau_x(Intersection_Data& inter_node, CubicDoub& uh, Beta& beta)
{
    int ix, iy, iz;                                                 //index of left node of intersection
    double upper_distance, lower_distance;                          //distance for eta approximation; distance for comparing outside and inside
    double x_upper, x_lower, upper_value, lower_value;              //x coordinate and value of two auxiliary points
    double x0, u_tau;                                               //u_{tau}^{+} or u_{tau}^{-}
    VecDoub upper, lower, tau;                                      //absolute location for approximating two auxiliary points and u_{tau}
    MatrixDoub wei_upper, wei_lower, wei_tau;                       //weights for two auxiliary points and weights for u_{tau}
    double jumpbeta_tau;
    
    ix = inter_node.left_loc;
    iy = inter_node.line.indx1;
    iz = inter_node.line.indx2;
    
    x0 = inter_node.coord.x_value;
    
    x_upper = inter_node.tau.auxl[0];
    x_lower = inter_node.tau.auxl[1];
    
    //Initialize two auxiliary points with 0
    upper_value = 0;
    lower_value = 0;
    
    upper_distance = sqrt(abs(x_upper-x0)*abs(x_upper-x0)+dy*dy);
    lower_distance = sqrt(abs(x_lower-x0)*abs(x_lower-x0)+dy*dy);
    
    /*
     upper_distance = sqrt(abs(x_upper-x0)*abs(x_upper-x0)*inter_node.p[2][0]*inter_node.p[2][0]+
     dy*dy*inter_node.p[2][1]*inter_node.p[2][1]);
     lower_distance = sqrt(abs(x_lower-x0)*abs(x_lower-x0)*inter_node.p[2][0]*inter_node.p[2][0]+
     dy*dy*inter_node.p[2][1]*inter_node.p[2][1]);
     */
    
    //Get Tau derivative along y-direction, make sure it's the same direction of given vector
    if(inter_node.p[2][1] > 0)
    {
        tau.push_back(-lower_distance);
        tau.push_back(0.0);
        tau.push_back(upper_distance);
    }
    else if(inter_node.p[2][1] < 0)
    {
        tau.push_back(lower_distance);
        tau.push_back(0.0);
        tau.push_back(-upper_distance);
    }
    else
    {
        cout << "x0 is on x-axis" << endl;
        exit(0);
    }
    
    wei_tau.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_tau[i].resize(2);
    }
    Weights(0,tau,3,1,wei_tau);
    
    //Weights for upper and lower auxiliary points approximation
    wei_upper.resize(3);
    wei_lower.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_upper[i].resize(1);
        wei_lower[i].resize(1);
    }
    
    if((REG == 'o') || (inter_node.tau.region == 'o'))
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(xi[inter_node.tau.uout_auxlnodes[i]]);
            lower.push_back(xi[inter_node.tau.lout_auxlnodes[i]]);
        }
        
        Weights(x_upper,upper,3,0,wei_upper);
        Weights(x_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[inter_node.tau.uout_auxlnodes[i]][iy+1][iz]*wei_upper[i][0];
            lower_value += uh[inter_node.tau.lout_auxlnodes[i]][iy-1][iz]*wei_lower[i][0];
        }
        
        u_tau = wei_tau[0][1]*lower_value+wei_tau[2][1]*upper_value;
        
        jumpbeta_tau = Tau_Plus(inter_node,beta,u_tau);
    }
    else if(inter_node.tau.region == 'i')
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(xi[inter_node.tau.uin_auxlnodes[i]]);
            lower.push_back(xi[inter_node.tau.lin_auxlnodes[i]]);
        }
        
        Weights(x_upper,upper,3,0,wei_upper);
        Weights(x_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[inter_node.tau.uin_auxlnodes[i]][iy+1][iz]*wei_upper[i][0];
            lower_value += uh[inter_node.tau.lin_auxlnodes[i]][iy-1][iz]*wei_lower[i][0];
        }
        
        u_tau = wei_tau[0][1]*lower_value+wei_tau[2][1]*upper_value;
        
        jumpbeta_tau = Tau_Minus(inter_node,beta,u_tau);
    }
    else
    {
        cout << "Approximation region isn't found, should only be 'o' for outside or 'i' for inside" << endl;
        exit(0);
    }
    
    return jumpbeta_tau;
}

/*******************************************************************************************
 Calculating of jump beta Eta for given intersection on y-direction based on current solution
 
 INPUT
 inter_node : given intersection node
 uh         : current solution
 beta       : beta object
 
 OUTPUT
 jumpbeta_eta : jump beta eta for given intersection based on current solution
 *******************************************************************************************/
double Intersections::Eta_y(Intersection_Data& inter_node, CubicDoub& uh, Beta& beta)
{
    int ix,iz,iy;                                               //index of left node of intersection
    double upper_distance, lower_distance;                      //distance for eta approximation; distance for comparing outside and inside
    double y_upper, y_lower, upper_value, lower_value;          //y coordinate and value of two auxiliary points
    double y0, u_eta;                                           //y coordinate of intersections; u_{eta}^{+} or u_{eta}^{-}
    VecDoub upper, lower, eta;                                  //absolute location for approximating two auxiliary points and u_{eta}
    MatrixDoub wei_upper, wei_lower, wei_eta;                   //weights for two auxiliary points and weights for u_{eta}
    double jumpbeta_eta;
    
    ix = inter_node.line.indx1;
    iy = inter_node.left_loc;
    iz = inter_node.line.indx2;
    
    y0 = inter_node.coord.y_value;
    
    y_upper = inter_node.eta.auxl[0];
    y_lower = inter_node.eta.auxl[1];
    
    //Initialize two auxiliary points with 0
    upper_value = 0;
    lower_value = 0;
    
    upper_distance = sqrt(abs(y_upper-y0)*abs(y_upper-y0)+dz*dz);
    lower_distance = sqrt(abs(y_lower-y0)*abs(y_lower-y0)+dz*dz);
    
    /*
     upper_distance = sqrt(abs(y_upper-y0)*abs(y_upper-y0)*inter_node.p[1][1]*inter_node.p[1][1]+
     dz*dz*inter_node.p[1][2]*inter_node.p[1][2]);
     lower_distance = sqrt(abs(y_lower-y0)*abs(y_lower-y0)*inter_node.p[1][1]*inter_node.p[1][1]+
     dz*dz*inter_node.p[1][2]*inter_node.p[1][2]);
     */
    
    //Get Eta derivative along z-direction, make sure it's the same direction of given vector
    if(inter_node.p[1][2] > 0)
    {
        eta.push_back(-lower_distance);
        eta.push_back(0.0);
        eta.push_back(upper_distance);
    }
    else if(inter_node.p[1][2] < 0)
    {
        eta.push_back(lower_distance);
        eta.push_back(0.0);
        eta.push_back(-upper_distance);
    }
    else
    {
        cout << "y0 is on y-axis" << endl;
        exit(0);
    }
    
    wei_eta.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_eta[i].resize(2);
    }
    Weights(0,eta,3,1,wei_eta);
    
    //Weights for upper and lower auxiliary points approximation
    wei_upper.resize(3);
    wei_lower.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_upper[i].resize(1);
        wei_lower[i].resize(1);
    }
    
    if((REG == 'o') || (inter_node.eta.region == 'o'))
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(yi[inter_node.eta.uout_auxlnodes[i]]);
            lower.push_back(yi[inter_node.eta.lout_auxlnodes[i]]);
        }
        
        Weights(y_upper,upper,3,0,wei_upper);
        Weights(y_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix][inter_node.eta.uout_auxlnodes[i]][iz+1]*wei_upper[i][0];
            lower_value += uh[ix][inter_node.eta.lout_auxlnodes[i]][iz-1]*wei_lower[i][0];
        }
        
        u_eta = wei_eta[0][1]*lower_value+wei_eta[2][1]*upper_value;
        
        jumpbeta_eta = Eta_Plus(inter_node,beta,u_eta);
    }
    else if(inter_node.eta.region == 'i')
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(yi[inter_node.eta.uin_auxlnodes[i]]);
            lower.push_back(yi[inter_node.eta.lin_auxlnodes[i]]);
        }
        
        Weights(y_upper,upper,3,0,wei_upper);
        Weights(y_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix][inter_node.eta.uin_auxlnodes[i]][iz+1]*wei_upper[i][0];
            lower_value += uh[ix][inter_node.eta.lin_auxlnodes[i]][iz-1]*wei_lower[i][0];
        }
        
        u_eta = wei_eta[0][1]*lower_value+wei_eta[2][1]*upper_value;
        
        jumpbeta_eta = Eta_Minus(inter_node,beta,u_eta);
    }
    else
    {
        cout << "Approximation region isn't found, should only be 'o' for outside or 'i' for inside" << endl;
        exit(0);
    }
    
    return jumpbeta_eta;
}

/*******************************************************************************************
 Calculating of jump beta Tau for given intersection on y-direction based on current solution
 
 INPUT
 inter_node : given intersection node
 uh         : current solution
 beta       : beta object
 
 OUTPUT
 jumpbeta_tau : jump beta tau for given intersection based on current solution
 *******************************************************************************************/
double Intersections::Tau_y(Intersection_Data& inter_node, CubicDoub& uh, Beta& beta)
{
    int ix,iz,iy;                                               //index of left node of intersection
    double upper_distance, lower_distance;                      //distance for eta approximation; distance for comparing outside and inside
    double y_upper, y_lower, upper_value, lower_value;          //y coordinate and value of two auxiliary points
    double y0, u_tau;                                           //y coordinate of intersections; u_{tau}^{+} or u_{tau}^{-}
    VecDoub upper, lower, tau;                                  //absolute location for approximating two auxiliary points and u_{tau}
    MatrixDoub wei_upper, wei_lower, wei_tau;                   //weights for two auxiliary points and weights for u_{tau}
    double jumpbeta_tau;
    
    ix = inter_node.line.indx1;
    iy = inter_node.left_loc;
    iz = inter_node.line.indx2;
    
    y0 = inter_node.coord.y_value;
    
    y_upper = inter_node.tau.auxl[0];
    y_lower = inter_node.tau.auxl[1];
    
    //Initialize two auxiliary points with 0
    upper_value = 0;
    lower_value = 0;
    
    upper_distance = sqrt(abs(y_upper-y0)*abs(y_upper-y0)+dx*dx);
    lower_distance = sqrt(abs(y_lower-y0)*abs(y_lower-y0)+dx*dx);
    
    /*
     upper_distance = sqrt(abs(y_upper-y0)*abs(y_upper-y0)*inter_node.p[2][1]*inter_node.p[2][1]+
     dx*dx*inter_node.p[2][0]*inter_node.p[2][0]);
     lower_distance = sqrt(abs(y_lower-y0)*abs(y_lower-y0)*inter_node.p[2][1]*inter_node.p[2][1]+
     dx*dx*inter_node.p[2][0]*inter_node.p[2][0]);
     */
    
    //Get Tau derivative along x-direction, make sure it's the same direction of given vector
    if(inter_node.p[2][0] > 0)
    {
        tau.push_back(-lower_distance);
        tau.push_back(0.0);
        tau.push_back(upper_distance);
    }
    else if(inter_node.p[2][0] < 0)
    {
        tau.push_back(lower_distance);
        tau.push_back(0.0);
        tau.push_back(-upper_distance);
    }
    else
    {
        cout << "y0 is on y-axis" << endl;
        exit(0);
    }
    
    wei_tau.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_tau[i].resize(2);
    }
    Weights(0,tau,3,1,wei_tau);
    
    //Weights for upper and lower auxiliary points approximation
    wei_upper.resize(3);
    wei_lower.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_upper[i].resize(1);
        wei_lower[i].resize(1);
    }
    
    
    if((REG == 'o') || (inter_node.tau.region == 'o'))
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(yi[inter_node.tau.uout_auxlnodes[i]]);
            lower.push_back(yi[inter_node.tau.lout_auxlnodes[i]]);
        }
        
        Weights(y_upper,upper,3,0,wei_upper);
        Weights(y_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix+1][inter_node.tau.uout_auxlnodes[i]][iz]*wei_upper[i][0];
            lower_value += uh[ix-1][inter_node.tau.lout_auxlnodes[i]][iz]*wei_lower[i][0];
        }
        
        u_tau = wei_tau[0][1]*lower_value+wei_tau[2][1]*upper_value;
        
        jumpbeta_tau = Tau_Plus(inter_node,beta,u_tau);
    }
    else if(inter_node.tau.region == 'i')
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(yi[inter_node.tau.uin_auxlnodes[i]]);
            lower.push_back(yi[inter_node.tau.lin_auxlnodes[i]]);
        }
        
        Weights(y_upper,upper,3,0,wei_upper);
        Weights(y_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix+1][inter_node.tau.uin_auxlnodes[i]][iz]*wei_upper[i][0];
            lower_value += uh[ix-1][inter_node.tau.lin_auxlnodes[i]][iz]*wei_lower[i][0];
        }
        
        u_tau = wei_tau[0][1]*lower_value+wei_tau[2][1]*upper_value;
        
        jumpbeta_tau = Tau_Minus(inter_node,beta,u_tau);
    }
    else
    {
        cout << "Approximation region isn't found, should only be 'o' for outside or 'i' for inside" << endl;
        exit(0);
    }
    
    return jumpbeta_tau;
}

/*******************************************************************************************
 Calculating of jump beta Eta for given intersection on z-direction based on current solution
 
 INPUT
 inter_node : given intersection node
 uh         : current solution
 beta       : beta object
 
 OUTPUT
 jumpbeta_eta : jump beta eta for given intersection based on current solution
 *******************************************************************************************/
double Intersections::Eta_z(Intersection_Data& inter_node, CubicDoub& uh, Beta& beta)
{
    int ix,iz,iy;                                               //index of left node of intersection
    double upper_distance, lower_distance;                      //distance for eta approximation; distance for comparing outside and inside
    double z_upper, z_lower, upper_value, lower_value;          //z coordinate and value of two auxiliary points
    double z0, u_eta;                                           //z coordinate of intersections; u_{eta}^{+} or u_{eta}^{-}
    VecDoub upper, lower, eta;                                  //absolute location for approximating two auxiliary points and u_{eta}
    MatrixDoub wei_upper, wei_lower, wei_eta;                   //weights for two auxiliary points and weights for u_{eta}
    double jumpbeta_eta;
    
    ix = inter_node.line.indx1;
    iy = inter_node.line.indx2;
    iz = inter_node.left_loc;
    
    z0 = inter_node.coord.z_value;
    
    z_upper = inter_node.eta.auxl[0];
    z_lower = inter_node.eta.auxl[1];
    
    //Initialize two auxiliary points with 0
    upper_value = 0;
    lower_value = 0;
    
    upper_distance = sqrt(abs(z_upper-z0)*abs(z_upper-z0)+dy*dy);
    lower_distance = sqrt(abs(z_lower-z0)*abs(z_lower-z0)+dy*dy);
    
    /*
     upper_distance = sqrt(abs(z_upper-z0)*abs(z_upper-z0)*inter_node.p[1][2]*inter_node.p[1][2]+
     dy*dy*inter_node.p[1][1]*inter_node.p[1][1]);
     lower_distance = sqrt(abs(z_lower-z0)*abs(z_lower-z0)*inter_node.p[1][2]*inter_node.p[1][2]+
     dy*dy*inter_node.p[1][1]*inter_node.p[1][1]);
     */
    
    //Get Eta derivative along y-direction, make sure it's the same direction of given vector
    if(inter_node.p[1][1] > 0)
    {
        eta.push_back(-lower_distance);
        eta.push_back(0.0);
        eta.push_back(upper_distance);
    }
    else if(inter_node.p[1][1] < 0)
    {
        eta.push_back(lower_distance);
        eta.push_back(0.0);
        eta.push_back(-upper_distance);
    }
    else
    {
        cout << "z0 is on z-axis" << endl;
        exit(0);
    }
    
    wei_eta.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_eta[i].resize(2);
    }
    Weights(0,eta,3,1,wei_eta);
    
    //Weights for upper and lower auxiliary points approximation
    wei_upper.resize(3);
    wei_lower.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_upper[i].resize(1);
        wei_lower[i].resize(1);
    }
    
    
    if((REG == 'o') || (inter_node.eta.region == 'o'))
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(zi[inter_node.eta.uout_auxlnodes[i]]);
            lower.push_back(zi[inter_node.eta.lout_auxlnodes[i]]);
        }
        
        Weights(z_upper,upper,3,0,wei_upper);
        Weights(z_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix][iy+1][inter_node.eta.uout_auxlnodes[i]]*wei_upper[i][0];
            lower_value += uh[ix][iy-1][inter_node.eta.lout_auxlnodes[i]]*wei_lower[i][0];
        }
        
        u_eta = wei_eta[0][1]*lower_value+wei_eta[2][1]*upper_value;
        
        jumpbeta_eta = Eta_Plus(inter_node,beta,u_eta);
    }
    else if(inter_node.eta.region == 'i')
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(zi[inter_node.eta.uin_auxlnodes[i]]);
            lower.push_back(zi[inter_node.eta.lin_auxlnodes[i]]);
        }
        
        Weights(z_upper,upper,3,0,wei_upper);
        Weights(z_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix][iy+1][inter_node.eta.uin_auxlnodes[i]]*wei_upper[i][0];
            lower_value += uh[ix][iy-1][inter_node.eta.lin_auxlnodes[i]]*wei_lower[i][0];
        }
        
        u_eta = wei_eta[0][1]*lower_value+wei_eta[2][1]*upper_value;
        
        jumpbeta_eta = Eta_Minus(inter_node,beta,u_eta);
    }
    else
    {
        cout << "Approximation region isn't found, should only be 'o' for outside or 'i' for inside" << endl;
        exit(0);
    }
    
    return jumpbeta_eta;
}

/*******************************************************************************************
 Calculating of jump beta Tau for given intersection on z-direction based on current solution
 
 INPUT
 inter_node : given intersection node
 uh         : current solution
 beta       : beta object
 
 OUTPUT
 jumpbeta_tau : jump beta tau for given intersection based on current solution
 *******************************************************************************************/
double Intersections::Tau_z(Intersection_Data& inter_node, CubicDoub& uh, Beta& beta)
{
    int ix,iz,iy;                                               //index of left node of intersection
    double upper_distance, lower_distance;                      //distance for eta approximation; distance for comparing outside and inside
    double z_upper, z_lower, upper_value, lower_value;          //z coordinate and value of two auxiliary points
    double z0, u_tau;                                           //z coordinate of intersections; u_{tau}^{+} or u_{tau}^{-}
    VecDoub upper, lower, tau;                                  //absolute location for approximating two auxiliary points and u_{tau}
    MatrixDoub wei_upper, wei_lower, wei_tau;                   //weights for two auxiliary points and weights for u_{tau}
    
    VecDoub x_vec, z_vec;
    MatrixDoub wei_x, wei_z;
    double x_dir, z_dir;
    double jumpbeta_tau;
    
    ix = inter_node.line.indx1;
    iy = inter_node.line.indx2;
    iz = inter_node.left_loc;
    
    z0 = inter_node.coord.z_value;
    
    z_upper = inter_node.tau.auxl[0];
    z_lower = inter_node.tau.auxl[1];
    
    //Initialize two auxiliary points with 0
    upper_value = 0;
    lower_value = 0;
    
    upper_distance = sqrt(abs(z_upper-z0)*abs(z_upper-z0)+dx*dx);
    lower_distance = sqrt(abs(z_lower-z0)*abs(z_lower-z0)+dx*dx);
    
    /*
     upper_distance = sqrt(abs(z_upper-z0)*abs(z_upper-z0)*inter_node.p[2][2]*inter_node.p[2][2]+
     dx*dx*inter_node.p[2][0]*inter_node.p[2][0]);
     lower_distance = sqrt(abs(z_lower-z0)*abs(z_lower-z0)*inter_node.p[2][2]*inter_node.p[2][2]+
     dx*dx*inter_node.p[2][0]*inter_node.p[2][0]);
     */
    
    //Get Tau derivative along x-direction, make sure it's the same direction of given vector
    if(inter_node.p[2][0] > 0)
    {
        tau.push_back(-lower_distance);
        tau.push_back(0.0);
        tau.push_back(upper_distance);
    }
    else if(inter_node.p[2][0] < 0)
    {
        tau.push_back(lower_distance);
        tau.push_back(0.0);
        tau.push_back(-upper_distance);
    }
    else
    {
        cout << "z0 is on z-axis" << endl;
        exit(0);
    }
    
    wei_tau.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_tau[i].resize(2);
    }
    Weights(0,tau,3,1,wei_tau);
    
    //Weights for upper and lower auxiliary points approximation
    wei_upper.resize(3);
    wei_lower.resize(3);
    for(int i = 0; i < 3; i++)
    {
        wei_upper[i].resize(1);
        wei_lower[i].resize(1);
    }
    
    if((REG == 'o') || (inter_node.tau.region == 'o'))
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(zi[inter_node.tau.uout_auxlnodes[i]]);
            lower.push_back(zi[inter_node.tau.lout_auxlnodes[i]]);
        }
        
        Weights(z_upper,upper,3,0,wei_upper);
        Weights(z_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix+1][iy][inter_node.tau.uout_auxlnodes[i]]*wei_upper[i][0];
            lower_value += uh[ix-1][iy][inter_node.tau.lout_auxlnodes[i]]*wei_lower[i][0];
        }
        
        u_tau = wei_tau[0][1]*lower_value+wei_tau[2][1]*upper_value;
        
        jumpbeta_tau = Tau_Plus(inter_node,beta,u_tau);
    }
    else if(inter_node.tau.region == 'i')
    {
        for(int i = 0; i < 3; i++)
        {
            upper.push_back(zi[inter_node.tau.uin_auxlnodes[i]]);
            lower.push_back(zi[inter_node.tau.lin_auxlnodes[i]]);
        }
        
        Weights(z_upper,upper,3,0,wei_upper);
        Weights(z_lower,lower,3,0,wei_lower);
        
        for(int i = 0; i < 3; i++)
        {
            upper_value += uh[ix+1][iy][inter_node.tau.uin_auxlnodes[i]]*wei_upper[i][0];
            lower_value += uh[ix-1][iy][inter_node.tau.lin_auxlnodes[i]]*wei_lower[i][0];
        }
        
        u_tau = wei_tau[0][1]*lower_value+wei_tau[2][1]*upper_value;
        
        jumpbeta_tau = Tau_Minus(inter_node,beta,u_tau);
    }
    else
    {
        cout << "Approximation region isn't found, should only be 'o' for outside or 'i' for inside" << endl;
        exit(0);
    }
    
    return jumpbeta_tau;
}

/*****************************************************************
 Calculating of jump beta Eta with given Eta^{+}
 
 INPUT
 inter_node : given intersection node
 beta       : beta object
 eta_plus   : jump eta^{+}
 
 OUTPUT
 eta_plus : jump beta eta^{+}
 *****************************************************************/
double Intersections::Eta_Plus(Intersection_Data& inter_node, Beta& beta, Doub_I eta_plus)
{
    double temp;
    
    temp = beta.Inside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)*inter_node.jump.u_eta+
          (beta.Outside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)-
           beta.Inside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value))*eta_plus;
    
    return temp;
}

/*****************************************************************
 Calculating of jump beta Eta with given Eta^{-}
 
 INPUT
 inter_node : given intersection node
 beta       : beta object
 eta_plus   : jump eta^{-}
 
 OUTPUT
 eta_plus : jump beta eta^{-}
 *****************************************************************/
double Intersections::Eta_Minus(Intersection_Data& inter_node, Beta& beta, Doub_I eta_minus)
{
    double temp;
    
    temp = beta.Outside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)*inter_node.jump.u_eta+
          (beta.Outside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)-
           beta.Inside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value))*eta_minus;
    
    return temp;
}

/*****************************************************************
 Calculating of jump beta Tau with given Tau^{+}
 
 INPUT
 inter_node : given intersection node
 beta       : beta object
 tau_plus   : jump tau^{+}
 
 OUTPUT
 tau_plus : jump beta tau^{+}
 *****************************************************************/
double Intersections::Tau_Plus(Intersection_Data& inter_node, Beta& beta, Doub_I tau_plus)
{
    double temp;
    
    temp = beta.Inside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)*inter_node.jump.u_tau+
          (beta.Outside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)-
           beta.Inside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value))*tau_plus;

    return temp;
}

/*****************************************************************
 Calculating of jump beta Tau with giveN Tau^{-}
 
 INPUT
 inter_node : given intersection node
 beta       : beta object
 tau_plus   : jump tau^{-}
 
 OUTPUT
 tau_plus : jump beta tau^{-}
 *****************************************************************/
double Intersections::Tau_Minus(Intersection_Data& inter_node, Beta& beta, Doub_I tau_minus)
{
    double temp;
    
    temp = beta.Outside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)*inter_node.jump.u_tau+
          (beta.Outside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value)-
           beta.Inside(inter_node.coord.x_value,inter_node.coord.y_value,inter_node.coord.z_value))*tau_minus;
    
    return temp;
}

/*****************************************************************
 Show all intersections' informations in x,y,z direction
 *****************************************************************/
void Intersections::Display()
{
    cout << endl << "----------------------------------- Intersections Info. Display -------------------------------" << endl;
    
    //Z-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iy = 1; iy < ny-1; iy++)
        {
            for(int ip = 0; ip < ifpz[ix][iy].size(); ip++)
            {
                if(ifpz[ix][iy][ip].ID != 0)
                {
                    Intersection_display(ifpz[ix][iy][ip],1);
                }
            }
        }
    }
    
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpy[ix][iz].size(); ip++)
            {
                if(ifpy[ix][iz][ip].ID != 0)
                {
                    Intersection_display(ifpy[ix][iz][ip],2);
                }
            }
        }
    }
    
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpx[iy][iz].size(); ip++)
            {
                if(ifpx[iy][iz][ip].ID != 0)
                {
                    Intersection_display(ifpx[iy][iz][ip],3);
                }
            }
        }
    }
}

/****************************************************************************
 Check fictitious points error, show in L2 norm and L^{infinity} norm
 ****************************************************************************/
void Intersections::Error_Jump(ofstream& out_file)
{
    double temp, sum, lmax, l2;
    int no;
    
    temp = 0;
    sum = 0;
    no = 0;
    
    //Z-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iy = 1; iy < ny-1; iy++)
        {
            for(int ip = 0; ip < ifpz[ix][iy].size(); ip++)
            {
                if(ifpz[ix][iy][ip].ID != 0)
                {
                    no += 1;
                    
                    sum += ifpz[ix][iy][ip].jump.err * ifpz[ix][iy][ip].jump.err;
                    
                    if(abs(ifpz[ix][iy][ip].jump.err) > temp)
                    {
                        temp = abs(ifpz[ix][iy][ip].jump.err);
                    }
                }
            }
        }
    }
    
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpy[ix][iz].size(); ip++)
            {
                if(ifpy[ix][iz][ip].ID != 0)
                {
                    no += 1;
                    
                    sum += ifpy[ix][iz][ip].jump.err * ifpy[ix][iz][ip].jump.err;
                    
                    if(abs(ifpy[ix][iz][ip].jump.err) > temp)
                    {
                        temp = abs(ifpy[ix][iz][ip].jump.err);
                    }
                }
            }
        }
    }
    
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpx[iy][iz].size(); ip++)
            {
                if(ifpx[iy][iz][ip].ID != 0)
                {
                    no += 1;
                    
                    sum += ifpx[iy][iz][ip].jump.err * ifpx[iy][iz][ip].jump.err;
                    
                    if(abs(ifpx[iy][iz][ip].jump.err) > temp)
                    {
                        temp = abs(ifpx[iy][iz][ip].jump.err);
                    }
                }
            }
        }
    }
    
    lmax = temp;
    l2 = sqrt(sum/no);
    
    out_file << "------------------------ Error of Jumps ------------------------" << endl;
    out_file << setprecision(3) << scientific;
    out_file << "Lmax = " << lmax << endl;
    out_file << "L2   = " << l2 << endl;
    out_file << fixed << endl;
}

/****************************************************************************
 Check fictitious points error, show in L2 norm and L^{infinity} norm
 ****************************************************************************/
void Intersections::Error_Fp(ofstream& out_file)
{
    double temp, sum, lmax, l2;
    int no, fpno;                  //fpno : numbers of fictitious point pairs
    
    temp = 0;
    sum = 0;
    no = 0;
    
    //Z-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iy = 1; iy < ny-1; iy++)
        {
            for(int ip = 0; ip < ifpz[ix][iy].size(); ip++)
            {
                if(ifpz[ix][iy][ip].ID != 0)
                {
                    fpno = (int)ifpz[ix][iy][ip].err.errr.size();
                    if(fpno != (int)ifpz[ix][iy][ip].err.errl.size())
                    {
                        cout << "Bad size for irregular MIB error" << endl;
                        exit(0);
                    }
                    
                    for(int i = 0; i < fpno; i++)
                    {
                        no += 2;
                        
                        sum += ifpz[ix][iy][ip].err.errl[i] * ifpz[ix][iy][ip].err.errl[i];
                        sum += ifpz[ix][iy][ip].err.errr[i] * ifpz[ix][iy][ip].err.errr[i];
                        
                        if(abs(ifpz[ix][iy][ip].err.errl[i]) > temp)
                        {
                            temp = abs(ifpz[ix][iy][ip].err.errl[i]);
                        }
                        if(abs(ifpz[ix][iy][ip].err.errr[i]) > temp)
                        {
                            temp = abs(ifpz[ix][iy][ip].err.errr[i]);
                        }
                    }
                }
            }
        }
    }
    
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpy[ix][iz].size(); ip++)
            {
                if(ifpy[ix][iz][ip].ID != 0)
                {
                    fpno = (int)ifpy[ix][iz][ip].err.errr.size();
                    if(fpno != (int)ifpy[ix][iz][ip].err.errl.size())
                    {
                        cout << "Bad size for irregular MIB error" << endl;
                        exit(0);
                    }
                    
                    for(int i = 0; i < fpno; i++)
                    {
                        no += 2;
                        
                        sum += ifpy[ix][iz][ip].err.errl[i] * ifpy[ix][iz][ip].err.errl[i];
                        sum += ifpy[ix][iz][ip].err.errr[i] * ifpy[ix][iz][ip].err.errr[i];
                        
                        if(abs(ifpy[ix][iz][ip].err.errl[i]) > temp)
                        {
                            temp = abs(ifpy[ix][iz][ip].err.errl[i]);
                        }
                        if(abs(ifpy[ix][iz][ip].err.errr[i]) > temp)
                        {
                            temp = abs(ifpy[ix][iz][ip].err.errr[i]);
                        }
                    }
                }
            }
        }
    }
    
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpx[iy][iz].size(); ip++)
            {
                if(ifpx[iy][iz][ip].ID != 0)
                {
                    fpno = (int)ifpx[iy][iz][ip].err.errr.size();
                    if(fpno != (int)ifpx[iy][iz][ip].err.errl.size())
                    {
                        cout << "Bad size for irregular MIB error" << endl;
                        exit(0);
                    }
                    
                    for(int i = 0; i < fpno; i++)
                    {
                        no += 2;
                        
                        sum += ifpx[iy][iz][ip].err.errl[i] * ifpx[iy][iz][ip].err.errl[i];
                        sum += ifpx[iy][iz][ip].err.errr[i] * ifpx[iy][iz][ip].err.errr[i];
                        
                        if(abs(ifpx[iy][iz][ip].err.errl[i]) > temp)
                        {
                            temp = abs(ifpx[iy][iz][ip].err.errl[i]);
                        }
                        if(abs(ifpx[iy][iz][ip].err.errr[i]) > temp)
                        {
                            temp = abs(ifpx[iy][iz][ip].err.errr[i]);
                        }
                    }
                }
            }
        }
    }
    
    lmax = temp;
    l2 = sqrt(sum/no);
    
    out_file << "------------------------ Error of FPs ------------------------" << endl;
    out_file << setprecision(3) << scientific;
    out_file << "Lmax = " << lmax << endl;
    out_file << "L2   = " << l2 << endl;
    out_file << fixed << endl;
}

/****************************************************
 Check size in z-direction
 ****************************************************/
void Intersections::Check_size_xy()
{
    cout << "On X-Y Plane, by Z-direction" << endl;
    cout << "x_size: " << ifpz.size() << " y_size: " << ifpz[0].size() << endl;
    for(int ix = 0; ix < ifpz.size(); ix++)
    {
        for(int iy = 0; iy < ifpz[ix].size(); iy++)
        {
            cout << "ix: " << ix << " iy: " << iy << " inode: " << ifpz[ix][iy].size() << endl;
        }
    }
}

/****************************************************
 Check size in y-direction
 ****************************************************/
void Intersections::Check_size_xz()
{
    cout << "On X-Z Plane, by Y-direction" << endl;
    cout << "x_size: " << ifpy.size() << " z_size: " << ifpy[0].size() << endl;
    for(int ix = 0; ix < ifpy.size(); ix++)
    {
        for(int iz = 0; iz < ifpy[ix].size(); iz++)
        {
            cout << "ix: " << ix << " iz: " << iz << " inode: " << ifpy[ix][iz].size() << endl;
        }
    }
}

/****************************************************
 Check size in x-direction
 ****************************************************/
void Intersections::Check_size_yz()
{
    cout << "On Y-Z Plane, by X-direction" << endl;
    cout << "y_size: " << ifpx.size() << " z_size: " << ifpx[0].size() << endl;
    for(int iy = 0; iy < ifpx.size(); iy++)
    {
        for(int iz = 0; iz < ifpx[iy].size(); iz++)
        {
            cout << "iy: " << iy << " iz: " << iz << " inode: " << ifpx[iy][iz].size() << endl;
        }
    }
}

/*************************************************************************
 Check coordinate approximation w.r.t given surface in x,y,z direction
 *************************************************************************/
void Intersections::Check_coord(Surface_Cartesian& ex)
{
    double temp;
    
    cout << "Z-Direction" << endl << endl;
    //Z-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iy = 1; iy < ny-1; iy++)
        {
            for(int ip = 0; ip < ifpz[ix][iy].size(); ip++)
            {
                temp = ex.check(ifpz[ix][iy][ip].coord.x_value,ifpz[ix][iy][ip].coord.y_value,ifpz[ix][iy][ip].coord.z_value);
                cout << temp << endl;
            }
        }
    }
    cout << endl;
    
    cout << "Y-Direction" << endl << endl;
    //Y-direction
    for(int ix = 1; ix < nx-1; ix++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpy[ix][iz].size(); ip++)
            {
                temp = ex.check(ifpy[ix][iz][ip].coord.x_value,ifpy[ix][iz][ip].coord.y_value,ifpy[ix][iz][ip].coord.z_value);
                cout << temp << endl;
            }
        }
    }
    cout << endl;
    
    cout << "X-Direction" << endl << endl;
    //X-direction
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpx[iy][iz].size(); ip++)
            {
                temp = ex.check(ifpx[iy][iz][ip].coord.x_value,ifpx[iy][iz][ip].coord.y_value,ifpx[iy][iz][ip].coord.z_value);
                cout << temp << endl;
            }
        }
    }
}

/*********************************************************************************************
 Show the smallest distance between intersection point and its left point in x-direction
 *********************************************************************************************/
void Intersections::Smallest_gamma_x()
{
    double temp;
    int xindx, yindx, zindx;
    
    temp = dx;
    
    for(int iy = 1; iy < ny-1; iy++)
    {
        for(int iz = 1; iz < nz-1; iz++)
        {
            for(int ip = 0; ip < ifpx[iy][iz].size(); ip++)
            {
                if(ifpx[iy][iz][ip].gamma < temp)
                {
                    xindx = ifpx[iy][iz][ip].left_loc;
                    yindx = ifpx[iy][iz][ip].line.indx1;
                    zindx = ifpx[iy][iz][ip].line.indx2;
                    temp = ifpx[iy][iz][ip].gamma;
                }
            }
        }
    }
    
    cout << "X index: " << xindx << endl;
    cout << "Y index: " << yindx << endl;
    cout << "Z index: " << zindx << endl;
    cout << "Smallest Gamma: " << temp << endl;
    cout << "Dx: " << dx << endl;
}

/************************************************************************
 Show all informations for one intersection node
 
 INPUT
 inter_node : given intersection node
 ************************************************************************/
void Intersections::Intersection_display(Intersection_Data inter_node, int ind)
{
    double dv;
    
    if(inter_node.ID < 0)
    {
        cout << "It's a corner interface!" << endl;
    }
    
    /*
     cout << "Direction: " << inter_node.dir << endl;
     cout << "X-value: " << inter_node.coord.x_value << " Y-value: " << inter_node.coord.y_value
     << " Z-value: " << inter_node.coord.z_value << endl;
     
     if(inter_node.dir == 'x')
     {
     cout << "Y-direction: " << inter_node.line.indx1 << " Z-direction: " << inter_node.line.indx2 << endl;
     }
     else if(inter_node.dir == 'y')
     {
     cout << "X-direction: " << inter_node.line.indx1 << " Z-direction: " << inter_node.line.indx2 << endl;
     }
     else if(inter_node.dir == 'z')
     {
     cout << "X-direction: " << inter_node.line.indx1 << " Y-direction: " << inter_node.line.indx2 << endl;
     }
     else
     {
     cout << "Error data for node's direction" << endl;
     exit(0);
     }
     
     cout << "Left grid point location: " << inter_node.left_loc << endl;
     
     //cout << "zenith: " << inter_node.angle.zenith << " azimuth: " << inter_node.angle.azimuth << endl;
     cout << "type: " << inter_node.type << endl;
     cout << "ID: " <<  inter_node.ID << endl;
     cout << "gamma: " << inter_node.gamma << endl;
     
     
     cout << "FP weights left: " << endl;
     for(int i = 0; i < inter_node.wei.weil.size(); i++)
     {
     for(int j = 0; j < inter_node.wei.weil[0].size(); j++)
     {
     cout << i << " " << j << ": " << inter_node.wei.weil[i][j] << " ";
     }
     }
     cout << endl;
     
     cout << "FP weights right: " << endl;
     for(int i = 0; i < inter_node.wei.weir.size(); i++)
     {
     for(int j = 0; j < inter_node.wei.weir[0].size(); j++)
     {
     cout << i << " " << j << ": " << inter_node.wei.weir[i][j] << " ";
     }
     }
     cout << endl;
     
     cout << "Error left: " << endl;
     for(int i = 0; i < inter_node.err.errr.size(); i++)
     {
     cout << i << ": " << inter_node.err.errl[i] << endl;
     }
     cout << "Error right: " << endl;
     for(int i = 0; i < inter_node.err.errl.size(); i++)
     {
     cout << i << ": " << inter_node.err.errr[i] << endl;
     }
     cout << endl;
     
     cout << "Jumps approximation relative error: " << inter_node.jump.err << endl;
     cout << endl;
     */
    
    
    
    if(abs(inter_node.jump.err) > 0.5)
    {
        cout << endl;
        cout << inter_node.dir << endl;
        
        cout << "Eta_error" << endl;
        cout << inter_node.jump.eta_err << endl;
        
        cout << "Tau_error" << endl;
        cout << inter_node.jump.tau_err << endl;
    }
    
    
    /*
     if(ind == 1)
     {
     dv = dz;
     }
     else if(ind == 2)
     {
     dv = dy;
     }
     else if(ind == 3)
     {
     dv = dx;
     }
     
     if(inter_node.jump.err > 3)
     {
     cout << endl;
     cout << inter_node.dir << endl;
     
     cout << "Eta Auxi Info." << endl;
     
     cout << "Upper" << endl;
     cout << (int)((inter_node.eta.auxl[0]+1.99)/dv) << endl;
     cout << "In" << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.eta.uin_auxlnodes[i] << " ";
     }
     cout << endl;
     cout << "Out" << endl;
     cout << (int)((inter_node.eta.auxl[0]+1.99)/dv) << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.eta.uout_auxlnodes[i] << " ";
     }
     cout << endl;
     
     cout << "Lower" << endl;
     cout << (int)((inter_node.eta.auxl[1]+1.99)/dv) << endl;
     cout << "In" << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.eta.lin_auxlnodes[i] << " ";
     }
     cout << endl;
     cout << "Out" << endl;
     cout << (int)((inter_node.eta.auxl[1]+1.99)/dv) << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.eta.lout_auxlnodes[i] << " ";
     }
     cout << endl;
     
     cout << "Region: " << inter_node.eta.region << endl << endl;
     
     cout << "Tau Auxi Info." << endl;
     
     cout << "Upper" << endl;
     cout << (int)((inter_node.tau.auxl[0]+1.99)/dv) << endl;
     cout << "In" << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.tau.uin_auxlnodes[i] << " ";
     }
     cout << endl;
     cout << "Out" << endl;
     cout << (int)((inter_node.tau.auxl[0]+1.99)/dv) << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.tau.uout_auxlnodes[i] << " ";
     }
     cout << endl;
     
     cout << "Lower" << endl;
     cout << (int)((inter_node.tau.auxl[1]+1.99)/dv) << endl;
     cout << "In" << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.tau.lin_auxlnodes[i] << " ";
     }
     cout << endl;
     cout << "Out" << endl;
     cout << (int)((inter_node.tau.auxl[1]+1.99)/dv) << endl;
     for(int i = 0; i < 3; i++)
     {
     cout << inter_node.tau.lout_auxlnodes[i] << " ";
     }
     cout << endl;
     
     cout << "Region: " << inter_node.tau.region << endl;
     }
     
     */
    
}

/******************************************************************************
 Calculate weights of Lagrange polynomial interpolation
 
 INPUT
 z: location where approximatetions are to be accurated
 x: x(n) grid points relative locations
 n: dimension of vector x(n), numbers of grid points
 m: highest derivative for which weghts are sought
 
 OUTPUT
 c: c(0:n-1,0:m) weights at grid locations x(n) for derivatives of order 0:m
 ******************************************************************************/
void Intersections::Weights(Doub_I z, VecDoub_I& x, Int_I n, Int_I m, MatrixDoub_O& c)
{
    double c1,c2,c3,c4,c5;
    int mn;
    
    c1 = 1;
    c4 = x[0]-z;
    
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m+1; j++)
        {
            c[i][j] = 0;
        }
    }
    
    c[0][0] = 1;
    
    for(int i = 1; i < n; i++)
    {
        mn = min(i,m);
        c2 = 1;
        c5 = c4;
        c4 = x[i] - z;
        for(int j = 0; j < i; j++)
        {
            c3 = x[i] - x[j];
            c2 = c2 * c3;
            if(j == i-1)
            {
                for(int k = mn; k > 0; k--)
                {
                    c[i][k] = c1 * (k * c[i-1][k-1] - c5 * c[i-1][k]) / c2;
                }
                c[i][0] = -c1 * c5 * c[i-1][0] / c2;
            }
            for(int k = mn; k > 0; k--)
            {
                c[j][k] = (c4 * c[j][k] - k * c[j][k-1]) / c3;
            }
            c[j][0] = c4 * c[j][0] / c3;
        }
        c1 = c2;
    }
}