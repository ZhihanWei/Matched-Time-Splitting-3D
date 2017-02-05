#ifndef __Surface_Cartesian_H__
#define __Surface_Cartesian_H__

#include <string>
#include "Constant.h"

class Surface_Cartesian
{
private:
    virtual double func(const Doub_I, const Doub_I, const Doub_I) const = 0;
    virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const = 0;
    virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const = 0;
    virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const = 0;
    
public:
    virtual double Set_Node(Doub_I, Doub_I, Doub_I) const = 0;
    virtual void display() = 0;
    virtual VecDoub normal(Doub_I, Doub_I, Doub_I) = 0;
    virtual double check(Doub_I, Doub_I, Doub_I) = 0;
    
    virtual double Gamma_x(Doub_I, Doub_I, Doub_I, Doub_I);
    virtual double Gamma_y(Doub_I, Doub_I, Doub_I, Doub_I);
    virtual double Gamma_z(Doub_I, Doub_I, Doub_I, Doub_I);
    
    double Rtsafe_x(const Doub_I, const Doub_I, const Doub_I, const Doub_I);
    double Rtsafe_y(const Doub_I, const Doub_I, const Doub_I, const Doub_I);
    double Rtsafe_z(const Doub_I, const Doub_I, const Doub_I, const Doub_I);
    
    //double Get_zenith(Doub_I, Doub_I, Doub_I);
    //double Get_azimuth(Doub_I, Doub_I);
    
};

#endif
