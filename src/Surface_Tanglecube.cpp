#include <iostream>
#include <cmath>
#include "Constant.h"
#include "Surface_Tanglecube.h"

using namespace std;

/****************************************************
                    Constructor
 
 INPUT
 arg : arguments of Tanglecube object
 *****************************************************/
Surface_Tanglecube::Surface_Tanglecube(VecDoub_I arg)
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
double Surface_Tanglecube::Set_Node(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double diff;
    
    diff = par[0]*pow(x,4) + par[1]*pow(x,2) + par[2]*pow(y,4) + par[3]*pow(y,2) + par[4]*pow(z,4) + par[5]*pow(z,2) + par[6] + TOL_SETUP;
    
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
double Surface_Tanglecube::func(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;
    
    temp = par[0]*pow(x,4) + par[1]*pow(x,2) + par[2]*pow(y,4) + par[3]*pow(y,2) + par[4]*pow(z,4) + par[5]*pow(z,2) + par[6];
    
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
double Surface_Tanglecube::dfx(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;
    
    temp = 4*par[0]*pow(x,3) + 2*par[1]*x;
    
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
double Surface_Tanglecube::dfy(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;
    
    temp = 4*par[2]*pow(y,3) + 2*par[3]*y;
    
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
double Surface_Tanglecube::dfz(const Doub_I x, const Doub_I y, const Doub_I z) const
{
    double temp;
    
    temp = 4*par[4]*pow(z,3) + 2*par[5]*z;
    
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
VecDoub Surface_Tanglecube::normal(Doub_I x, Doub_I y, Doub_I z)
{
    VecDoub normal;
    double total,temp;
    
    total = sqrt(pow(4*par[0]*pow(x,3)+2*par[1]*x,2)+pow(4*par[2]*pow(y,3)+2*par[3]*y,2)+pow(4*par[4]*pow(z,3)+2*par[5]*z,2));
    
    temp = (4*par[0]*pow(x,3)+2*par[1]*x)/total;
    normal.push_back(temp);
    temp = (4*par[2]*pow(y,3)+2*par[3]*y)/total;
    normal.push_back(temp);
    temp = (4*par[4]*pow(z,3)+2*par[5]*z)/total;
    normal.push_back(temp);
    
    return normal;
}

/**********************************************************************
 Show the quation of implicit function
 ***********************************************************************/
void Surface_Tanglecube::display()
{
    cout << "The equation of Tanglecube surface is:" << endl;
    cout << "f(x) = (" << par[0] << ")*x^4+(" << par[1] << ")*x^2+("
        << par[2] << ")*y^4+(" << par[3] << ")*y^2+(" << par[4] <<")*z^4+(" << par[5] << ")*z^2+(" << par[6] << ")=0" << endl;
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
double Surface_Tanglecube::check(Doub_I x, Doub_I y, Doub_I z)
{
    double temp;
    
    temp = par[0]*pow(x,4) + par[1]*pow(x,2) + par[2]*pow(y,4) + par[3]*pow(y,2) + par[4]*pow(z,4) + par[5]*pow(z,2) + par[6];
    
    return temp;
}