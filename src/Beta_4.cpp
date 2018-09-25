#include "Beta_4.h"
/******************************************************************
 Constructor
 ******************************************************************/
Beta_4::Beta_4() {}

/****************************************************
 Variable coefficient in Omega^{-}
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 Diffusion coefficient at that point in Omega^{-}
 *****************************************************/
double Beta_4::Inside(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = exp((x * x / (2 * in_sigma_x * in_sigma_x) + y * y / (2 * in_sigma_y * in_sigma_y) + z * z / (2 * in_sigma_z * in_sigma_z)));
	
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
double Beta_4::Outside(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = exp((x * x / (2 * out_sigma_x * out_sigma_x) + y * y / (2 * out_sigma_y * out_sigma_y) + z * z / (2 * out_sigma_z * out_sigma_z)));
	
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
double Beta_4::Inside_Dx(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = (x / in_sigma_x / in_sigma_x) * exp((x * x / (2 * in_sigma_x * in_sigma_x) + y * y / (2 * in_sigma_y * in_sigma_y) + z * z / (2 * in_sigma_z * in_sigma_z)));
	
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
double Beta_4::Outside_Dx(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = (x / out_sigma_x / out_sigma_x) * exp((x * x / (2 * out_sigma_x * out_sigma_x) + y * y / (2 * out_sigma_y * out_sigma_y) + z * z / (2 * out_sigma_z * out_sigma_z)));
	
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
double Beta_4::Inside_Dy(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = (y / in_sigma_y / in_sigma_y) * exp((x * x / (2 * in_sigma_x * in_sigma_x) + y * y / (2 * in_sigma_y * in_sigma_y) + z * z / (2 * in_sigma_z * in_sigma_z)));
	
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
double Beta_4::Outside_Dy(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = (y / out_sigma_y / out_sigma_y) * exp((x * x / (2 * out_sigma_x * out_sigma_x) + y * y / (2 * out_sigma_y * out_sigma_y) + z * z / (2 * out_sigma_z * out_sigma_z)));
	
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
double Beta_4::Inside_Dz(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = (z / in_sigma_z / in_sigma_z) * exp((x * x / (2 * in_sigma_x * in_sigma_x) + y * y / (2 * in_sigma_y * in_sigma_y) + z * z / (2 * in_sigma_z * in_sigma_z)));
	
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
double Beta_4::Outside_Dz(Doub_I x, Doub_I y, Doub_I z) const {
	double temp;
	
	temp = (z / out_sigma_z / out_sigma_z) * exp((x * x / (2 * out_sigma_x * out_sigma_x) + y * y / (2 * out_sigma_y * out_sigma_y) + z * z / (2 * out_sigma_z * out_sigma_z)));
	
	return temp;
}
