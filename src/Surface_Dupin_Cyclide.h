#ifndef Surface_Dupin_Cyclide_H__
#define Surface_Dupin_Cyclide_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/****************************************************************
 Dupin Cuclide Surface with parameter [1,0.3,0.6]
 Suggested domain: [-4,4;-4,4;-4,4] for grid (20,40,80,160)
 ****************************************************************/
class Surface_Dupin_Cyclide: public Surface_Cartesian
{
private:
    double a, b, c;
    
    virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
    virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
    virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
    virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;
public:
    Surface_Dupin_Cyclide(VecDoub_I);
    virtual double Set_Node(Doub_I, Doub_I, Doub_I) const;
    virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
    virtual void display();
    virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif

