#include "Eq_2.h"
#include <iostream>

/******************************************************************
 Constructor
 
 INPUT
 ti   : current time
 beta : vector of 2 double values represent beta^{-} and beta^{+}
 ******************************************************************/
Eq_2::Eq_2(Doub_I ti, Beta beta)
{
    t = ti;
    beta_minus = beta.in;
    beta_plus = beta.out;
}

/****************************************************
 Analytical solution in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of analytical solution at given node
 *****************************************************/
double Eq_2::Inner_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 10*exp(-x*x)*exp(-y*y)*exp(-z*z)+exp(a+t);
    
    return temp;
}

/****************************************************
 Analytical solution in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of analytical solution at given node
 *****************************************************/
double Eq_2::Outer_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 5*exp(-x*x)*exp(-y*y)*exp(-z*z)+exp(a+t);
    
    return temp;
}

/****************************************************
 Source term in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of source term at given node
 *****************************************************/
double Eq_2::Inner_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -beta_minus*10*(-6+4*(x*x+y*y+z*z))*exp(-x*x)*exp(-y*y)*exp(-z*z)+exp(a+t);
    
    return temp;
}

/****************************************************
 Source term in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of source term at given node
 *****************************************************/
double Eq_2::Outer_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -beta_plus*5*(-6+4*(x*x+y*y+z*z))*exp(-x*x)*exp(-y*y)*exp(-z*z)+exp(a+t);
    
    return temp;
}

/**********************************************************************
 1st derivative in x-direction of analytical solution in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of 1st derivative in x-direction at given node
 ***********************************************************************/
double Eq_2::Inner_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -2*10*x*exp(-x*x)*exp(-y*y)*exp(-z*z);
    
    return temp;
}

/**********************************************************************
 1st derivative in x-direction of analytical solution in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of 1st derivative in x-direction at given node
 ***********************************************************************/
double Eq_2::Outer_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -2*5*x*exp(-x*x)*exp(-y*y)*exp(-z*z);
    
    return temp;
}

/**********************************************************************
 1st derivative in y-direction of analytical solution in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of 1st derivative in y-direction at given node
 ***********************************************************************/
double Eq_2::Inner_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -2*10*y*exp(-x*x)*exp(-y*y)*exp(-z*z);
    
    return temp;
}

/**********************************************************************
 1st derivative in y-direction of analytical solution in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of 1st derivative in y-direction at given node
 ***********************************************************************/
double Eq_2::Outer_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -2*5*y*exp(-x*x)*exp(-y*y)*exp(-z*z);
    
    return temp;
}

/**********************************************************************
 1st derivative in z-direction of analytical solution in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of 1st derivative in z-direction at given node
 ***********************************************************************/
double Eq_2::Inner_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -2*10*z*exp(-x*x)*exp(-y*y)*exp(-z*z);
    
    return temp;
}

/**********************************************************************
 1st derivative in z-direction of analytical solution in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of 1st derivative in z-direction at given node
 ***********************************************************************/
double Eq_2::Outer_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -2*5*z*exp(-x*x)*exp(-y*y)*exp(-z*z);
    
    return temp;
}