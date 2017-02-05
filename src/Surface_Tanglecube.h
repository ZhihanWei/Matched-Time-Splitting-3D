#ifndef __Surface_Tanglecube_H__
#define __Surface_Tanglecube_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/************************************************************
 Tanglecube Surface with parameter [1,-5,1,-5,1,-5,10]
 Suggested domain: [-4,4;-4,4;-5,5]
 ************************************************************/
class Surface_Tanglecube: public Surface_Cartesian
{
private:
    double par[7];
    
    virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
    virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
    virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
    virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;
public:
    Surface_Tanglecube(VecDoub);
    virtual double Set_Node(const Doub_I, const Doub_I, const Doub_I) const;
    virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
    virtual void display();
    virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif
