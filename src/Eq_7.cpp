#include "Eq_7.h"

/******************************************************************
 Constructor
 
 INPUT
 ti   : current time
 beta : vector of 2 double values represent beta^{-} and beta^{+}
 ******************************************************************/
Eq_7::Eq_7(Doub_I ti, Beta beta)
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
double Eq_7::Inner_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = cos(k*x)*sin(k*y)*cos(k*z)+sin(t);
    
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
double Eq_7::Outer_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = sin(k*x)*cos(k*y)*sin(k*z)+sin(t);
    
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
double Eq_7::Inner_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 3*beta_minus*k*k*cos(k*x)*sin(k*y)*cos(k*z)+cos(t);
    
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
double Eq_7::Outer_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 3*beta_plus*k*k*sin(k*x)*cos(k*y)*sin(k*z)+cos(t);
    
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
double Eq_7::Inner_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -k*sin(k*x)*sin(k*y)*cos(k*z);
    
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
double Eq_7::Outer_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = k*cos(k*x)*cos(k*y)*sin(k*z);
    
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
double Eq_7::Inner_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = k*cos(k*x)*cos(k*y)*cos(k*z);
    
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
double Eq_7::Outer_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -k*sin(k*x)*sin(k*y)*sin(k*z);
    
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
double Eq_7::Inner_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -k*cos(k*x)*sin(k*y)*sin(k*z);
    
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
double Eq_7::Outer_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = k*sin(k*x)*cos(k*y)*cos(k*z);
    
    return temp;
}