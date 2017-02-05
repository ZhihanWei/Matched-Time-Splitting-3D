#include <cmath>
#include "Eq_0.h"

/******************************************************************
                       Constructor
 
 INPUT
 ti   : current time
 beta : vector of 2 double values represent beta^{-} and beta^{+}
 ******************************************************************/
Eq_0::Eq_0(Doub_I ti, Beta beta)
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
double Eq_0::Inner_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 80;
    
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
double Eq_0::Outer_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 80;
    
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
double Eq_0::Inner_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
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
double Eq_0::Outer_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
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
double Eq_0::Inner_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
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
double Eq_0::Outer_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
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
double Eq_0::Inner_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;

    temp = 0;
    
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
double Eq_0::Outer_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
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
double Eq_0::Inner_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
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
double Eq_0::Outer_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    
    return temp;
}