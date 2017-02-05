#include <iostream>
#include <cmath>
#include "Constant.h"
//#include "Surface_Cartesian.h"
#include "Surface_Heart.h"

using namespace std;

/****************************************************
                    Constructor
 
 INPUT
 arg : arguments of Heart object
 *****************************************************/
Surface_Heart::Surface_Heart(VecDoub_I arg)
{
    for(int i = 0; i < arg.size(); i++)
    {
        par[i] = arg[i];
    }
}

/**********************************************************************
 Implicit function with little perturbation
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of implicit function with perturbation
 ***********************************************************************/
double Surface_Heart::Set_Node(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double diff;
    
    diff = pow((par[0]*x*x + par[1]*y*y + par[2]*z*z + par[3]),3) + par[4]*x*x*z*z*z + par[5]*y*y*z*z*z + TOL_SETUP;
    
    return diff;
}

/*********************************************************************************
 Private member function to get value of implicit function without perturbation
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of function
 **********************************************************************************/
double Surface_Heart::func(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;
    
    temp = pow((par[0]*x*x + par[1]*y*y + par[2]*z*z + par[3]),3) + par[4]*x*x*z*z*z + par[5]*y*y*z*z*z;
    
    return temp;
}

/**********************************************************************
 Calculate derivative of function in x-direction
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 derivative of function in x-direction
 ***********************************************************************/
double Surface_Heart::dfx(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;
        
    temp = par[0]*6*x*pow((par[0]*x*x + par[1]*y*y + par[2]*z*z + par[3]),2) + par[4]*2*x*z*z*z;
        
    return temp;
}

/**********************************************************************
 Calculate derivative of function in y-direction
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 derivative of function in y-direction
 ***********************************************************************/
double Surface_Heart::dfy(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;

    temp = 6*par[1]*y*pow((par[0]*x*x + par[1]*y*y + par[2]*z*z + par[3]),2) + 2*par[5]*y*z*z*z;
        
    return temp;
}

/**********************************************************************
 Calculate derivative of function in z-direction
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 derivative of function in z-direction
 ***********************************************************************/
double Surface_Heart::dfz(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;

    temp = par[2]*6*z*pow((par[0]*x*x + par[1]*y*y + par[2]*z*z + par[3]),2) + par[4]*3*x*x*z*z + par[5]*3*y*y*z*z;
        
    return temp;
}

/**********************************************************************
 Get normal direction at given point
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 vector of normal direction with 3 components
 ***********************************************************************/
VecDoub Surface_Heart::normal(Doub_I x, Doub_I y, Doub_I z)
{
    VecDoub normal;
    double total,temp;
    
    total = sqrt(pow(par[0]*6*x*pow((par[0]*x*x+par[1]*y*y+par[2]*z*z+par[3]),2)+par[4]*2*x*z*z*z,2)+
                 pow(6*par[1]*y*pow((par[0]*x*x+par[1]*y*y+par[2]*z*z+par[3]),2)+2*par[5]*y*z*z*z,2)+
                 pow(par[2]*6*z*pow((par[0]*x*x+par[1]*y*y+par[2]*z*z+par[3]),2)+par[4]*3*x*x*z*z+par[5]*3*y*y*z*z,2));
    
    temp = (par[0]*6*x*pow((par[0]*x*x+par[1]*y*y+par[2]*z*z+par[3]),2)+par[4]*2*x*z*z*z)/total;
    normal.push_back(temp);
    temp = (6*par[1]*y*pow((par[0]*x*x+par[1]*y*y+par[2]*z*z+par[3]),2)+2*par[5]*y*z*z*z)/total;
    normal.push_back(temp);
    temp = (par[2]*6*z*pow((par[0]*x*x+par[1]*y*y+par[2]*z*z+par[3]),2)+par[4]*3*x*x*z*z+par[5]*3*y*y*z*z)/total;
    normal.push_back(temp);
    
    return normal;
}

/**********************************************************************
 Show the quation of implicit function
 ***********************************************************************/
void Surface_Heart::display()
{
    cout << "The equation of Heart surface is:" << endl;
    cout << "f(x) = [(" << par[0] << ")*x^2+(" << par[1] << ")*y^2+(" << par[2] <<")*z^2+(" << par[3] << ")]^3+(";
    cout << par[4] << ")*x^2*z^3+(" << par[5] << ")*y^2*z^3 = 0" << endl;
}

/*****************************************************************************************
 Public member function to get the value of implicit function without perturbation
 
 INPUT
 x : x coordinate of given point
 y : y coordinate of given point
 z : z coordinate of given point
 
 OUTPUT
 value of implicit function without perturbation
 ******************************************************************************************/
double Surface_Heart::check(Doub_I x, Doub_I y, Doub_I z)
{
    double temp;
    
    temp = pow((par[0]*x*x + par[1]*y*y + par[2]*z*z + par[3]),3) + par[4]*x*x*z*z*z + par[5]*y*y*z*z*z;
    
    return temp;
}