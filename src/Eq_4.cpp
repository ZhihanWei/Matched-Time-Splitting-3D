#include "Eq_4.h"

/******************************************************************
                            Constructor
 
 INPUT
 ti   : current time
 beta : vector of 2 double values represent beta^{-} and beta^{+}
 ******************************************************************/
Eq_4::Eq_4(Doub_I ti, Beta& in_beta):beta(in_beta)
{
    t = ti;
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
double Eq_4::Inner_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = cos(k*x)*sin(k*y)*cos(k*z)*cos(t)-1;
    
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
double Eq_4::Outer_u(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = sin(k*x)*cos(k*y)*sin(k*z)*cos(t)+1;
    
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
double Eq_4::Inner_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -sin(t)*cos(k*x)*sin(k*y)*cos(k*z)+
            3*beta.Inside(x,y,z)*k*k*cos(t)*cos(k*x)*sin(k*y)*cos(k*z)-
            beta.Inside_Dx(x,y,z)*Inner_dux(x,y,z)-
            beta.Inside_Dy(x,y,z)*Inner_duy(x,y,z)-
            beta.Inside_Dz(x,y,z)*Inner_duz(x,y,z);
    
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
double Eq_4::Outer_f(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -sin(t)*sin(k*x)*cos(k*y)*sin(k*z)+
            3*beta.Outside(x,y,z)*k*k*cos(t)*sin(k*x)*cos(k*y)*sin(k*z)-
            beta.Outside_Dx(x,y,z)*Outer_dux(x,y,z)-
            beta.Outside_Dy(x,y,z)*Outer_duy(x,y,z)-
            beta.Outside_Dz(x,y,z)*Outer_duz(x,y,z);
    
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
double Eq_4::Inner_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -k*sin(k*x)*sin(k*y)*cos(k*z)*cos(t);
    
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
double Eq_4::Outer_dux(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = k*cos(k*x)*cos(k*y)*sin(k*z)*cos(t);
    
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
double Eq_4::Inner_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = k*cos(k*x)*cos(k*y)*cos(k*z)*cos(t);
    
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
double Eq_4::Outer_duy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -k*sin(k*x)*sin(k*y)*sin(k*z)*cos(t);
    
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
double Eq_4::Inner_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -k*cos(k*x)*sin(k*y)*sin(k*z)*cos(t);
    
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
double Eq_4::Outer_duz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = k*sin(k*x)*cos(k*y)*cos(k*z)*cos(t);
    
    return temp;
}
/***************************************************************
 Analytical jump [beta_u] in x-direction
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 analutical of jump [beta_u_x]
 ***************************************************************/
double Eq_4::Jump_betau_x(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = beta.Outside(x,y,z)*Outer_dux(x,y,z) - beta.Inside(x,y,z)*Inner_dux(x,y,z);
    
    return temp;
}

/***************************************************************
 Analytical jump [beta_u] in y-direction
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 analutical of jump [beta_u_y]
 ***************************************************************/
double Eq_4::Jump_betau_y(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = beta.Outside(x,y,z)*Outer_duy(x,y,z) - beta.Inside(x,y,z)*Inner_duy(x,y,z);
    
    return temp;
}

/***************************************************************
 Analytical jump [beta_u] in z-direction
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 analutical of jump [beta_u_z]
 ***************************************************************/
double Eq_4::Jump_betau_z(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = beta.Outside(x,y,z)*Outer_duz(x,y,z) - beta.Inside(x,y,z)*Inner_duz(x,y,z);
    
    return temp;
}