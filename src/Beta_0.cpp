#include "Beta_0.h"

/******************************************************************
 Constructor
 ******************************************************************/
Beta_0::Beta_0() {}

/****************************************************
 Variable coefficient in Omega^{-}

 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point

 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 *****************************************************/
double Beta_0::Inside(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 1;

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
double Beta_0::Outside(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 10;

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
double Beta_0::Inside_Dx(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0;

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
double Beta_0::Outside_Dx(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0;

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
double Beta_0::Inside_Dy(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0;

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
double Beta_0::Outside_Dy(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0;

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
double Beta_0::Inside_Dz(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0;

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
double Beta_0::Outside_Dz(Doub_I x, Doub_I y, Doub_I z) const {
  double temp;

  temp = 0;

  return temp;
}
