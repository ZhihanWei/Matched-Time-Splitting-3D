#include "Beta_1.h"

/****************************************************
          Diffusion coefficient in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 *****************************************************/
double Beta_1::Inside(Doub_I x, Doub_I y, Doub_I z) const
{
    return 1/(x*x+y*y+z*z+1);
}

/****************************************************
         Diffusion coefficient in Omega^{+}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at the point in Omega^{+}
 *****************************************************/
double Beta_1::Outside(Doub_I x, Doub_I y, Doub_I z) const
{
    return -1/(x*x+y*y+z*z+1);
}