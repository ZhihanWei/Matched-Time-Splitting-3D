#include "Beta_3.h"

/******************************************************************
 Constructor
 ******************************************************************/
Beta_3::Beta_3()
{
    
}

/****************************************************
 Variable coefficient in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 *****************************************************/
double Beta_3::Inside(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp((x+y+z)/6);
    
    return temp;
}

/****************************************************
 Variable coefficient in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at the point in Omega^{+}
 *****************************************************/
double Beta_3::Outside(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp(-(x+y+z)/6);
    
    return temp;
}

/**************************************************************
 First derivative variable coefficient x-direction in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 ***************************************************************/
double Beta_3::Inside_Dx(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp((x+y+z)/6)/6;
    
    return temp;
}

/****************************************************
 First derivative variable coefficient x-direction in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at the point in Omega^{+}
 *****************************************************/
double Beta_3::Outside_Dx(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp(-(x+y+z)/6)/6;
    
    return temp;
}

/**************************************************************
 First derivative variable coefficient x-direction in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 ***************************************************************/
double Beta_3::Inside_Dy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp((x+y+z)/6)/6;
    
    return temp;
}

/****************************************************
 First derivative variable coefficient x-direction in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at the point in Omega^{+}
 *****************************************************/
double Beta_3::Outside_Dy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp(-(x+y+z)/6)/6;
    
    return temp;
}

/**************************************************************
 First derivative variable coefficient x-direction in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 ***************************************************************/
double Beta_3::Inside_Dz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp((x+y+z)/6)/6;
    
    return temp;
}

/****************************************************
 First derivative variable coefficient x-direction in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at the point in Omega^{+}
 *****************************************************/
double Beta_3::Outside_Dz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = exp(-(x+y+z)/6)/6;
    
    return temp;
}
