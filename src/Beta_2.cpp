#include "Beta_2.h"

/******************************************************************
 Constructor
 ******************************************************************/
Beta_2::Beta_2()
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
double Beta_2::Inside(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;

    temp = sin(x+y+z)+2;
    
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
double Beta_2::Outside(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -cos(x+y+z)+2;
    
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
double Beta_2::Inside_Dx(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = cos(x+y+z);
    
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
double Beta_2::Outside_Dx(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = sin(x+y+z);
    
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
double Beta_2::Inside_Dy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = cos(x+y+z);
    
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
double Beta_2::Outside_Dy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = sin(x+y+z);
    
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
double Beta_2::Inside_Dz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = cos(x+y+z);
    
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
double Beta_2::Outside_Dz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = sin(x+y+z);
    
    return temp;
}
