#ifndef __DATA_H__
#define __DATA_H__

#include "Constant.h"
#include "New_Data_Type.h"

class Data
{
private:
    enum Prt_name
    {
        e_TOL_ITYPE,
        e_xl,e_xr,e_yl,e_yr,e_zl,e_zr,
        e_t_start,e_t_finish,e_t_step,
        e_nx,e_ny,e_nz,
        e_beta_inside,e_beta_outside,
        e_surface,e_equation,
        e_accuracy
    };
    
    double TOL_ITYPE;                       //TOL_ITYPE is the multiplication of dx/dy/dz
    double beta_inside,beta_outside;
    double xl,xr,yl,yr,zl,zr;
    double t_start,t_finish,t_step;
    int nx,ny,nz;
    char surface;
    int equation;
    int accuracy;
    
    Prt_name Translation(const string& in_string);
public:
    Data(const string& file_name);
    
    void Display();
    
    VecDoub Get_Domain() const;
    VecInt Get_Size() const;
    Beta Get_Beta() const;
    VecDoub Get_Time() const;
    char Get_Surface() const;
    int Get_Equation() const;
    int Get_Accuracy() const;
    double Get_Tol() const;
};
#endif
