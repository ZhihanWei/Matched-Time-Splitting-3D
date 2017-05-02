#include "Beta_1.h"

/******************************************************************
                            Constructor
 ******************************************************************/
Beta_1::Beta_1()
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
double Beta_1::Inside(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = y;
    //temp = 1/(x*x+y*y+z*z+1);
    
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
double Beta_1::Outside(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -y;
    //temp = -1/(x*x+y*y+z*z+1);
    
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
double Beta_1::Inside_Dx(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    //temp= -2*x/(x*x+y*y+z*z+1)/(x*x+y*y+z*z+1);
    
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
double Beta_1::Outside_Dx(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    //temp = 2*x/(x*x+y*y+z*z+1)/(x*x+y*y+z*z+1);
    
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
double Beta_1::Inside_Dy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 1;
    //temp = -2*y/(x*x+y*y+z*z+1)/(x*x+y*y+z*z+1);
    
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
double Beta_1::Outside_Dy(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = -1;
    //temp = 2*y/(x*x+y*y+z*z+1)/(x*x+y*y+z*z+1);
    
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
double Beta_1::Inside_Dz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    //temp = -2*z/(x*x+y*y+z*z+1)/(x*x+y*y+z*z+1);
    
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
double Beta_1::Outside_Dz(Doub_I x, Doub_I y, Doub_I z) const
{
    double temp;
    
    temp = 0;
    //temp = 2*z/(x*x+y*y+z*z+1)/(x*x+y*y+z*z+1);
    
    return temp;
}