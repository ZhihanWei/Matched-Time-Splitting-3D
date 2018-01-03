#ifndef __Surface_Torus_H__
#define __Surface_Torus_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/************************************************************
 Torus Surface with parameter: [a, b]
 
 Equation : (a - sqrt(x^2 + y^2))^2 + z^2 - b^2 = 0
 
 Suggested domain: [-4,4;-4,4;-4,4] for grid (20,40,80,160)
 ************************************************************/
class Surface_Torus : public Surface_Cartesian {
 private:
  double a, b;

  virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;

 public:
  Surface_Torus(VecDoub_I);
  virtual double Set_Node(Doub_I, Doub_I, Doub_I) const;
  virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
  virtual void display();
  virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif
