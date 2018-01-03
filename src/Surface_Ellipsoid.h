#ifndef __Surface_Ellipsoid_H__
#define __Surface_Ellipsoid_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/************************************************************
 Ellipsoid Surface with parameter: [a, b, c, r]
 
 Equation : x^2 / a^2 + y^2 / b^2 + z^2 / c^2 - r^2 = 0
 
 Suggested domain: [-4,4;-2,2;-2,2] for grid (20,40,80,160)
 ************************************************************/
class Surface_Ellipsoid : public Surface_Cartesian {
 private:
  double radius;
  double a, b, c;

  virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;

 public:
  Surface_Ellipsoid(VecDoub_I);
  virtual double Set_Node(Doub_I, Doub_I, Doub_I) const;
  virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
  virtual void display();
  virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif
