#ifndef __Eq_1_H__
#define __Eq_1_H__

#include "Constant.h"
#include "New_Data_Type.h"
#include "Equation.h"

class Eq_1: public Equation
{
private:
    virtual double Inner_dux(Doub_I, Doub_I, Doub_I) const;
    virtual double Outer_dux(Doub_I, Doub_I, Doub_I) const;
    virtual double Inner_duy(Doub_I, Doub_I, Doub_I) const;
    virtual double Outer_duy(Doub_I, Doub_I, Doub_I) const;
    virtual double Inner_duz(Doub_I, Doub_I, Doub_I) const;
    virtual double Outer_duz(Doub_I, Doub_I, Doub_I) const;
public:
    Eq_1(Doub_I, Beta);
    virtual double Outer_u(Doub_I, Doub_I, Doub_I) const;                 //Outer source term
    virtual double Inner_u(Doub_I, Doub_I, Doub_I) const;                 //Inner source term
    virtual double Outer_f(Doub_I, Doub_I, Doub_I) const;                 //Outer source term
    virtual double Inner_f(Doub_I, Doub_I, Doub_I) const;                 //Inner source term
};

#endif
