#ifndef __Surface_Molecular_H__
#define __Surface_Molecular_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/************************************************************
 Molecular Surface
 
 Equation : (x^2 + y^2 + z^2 + 0.6)^2 - 3.5 * y^2 - 0.6 = 0
 
 Suggested domain: [-2,2;-3,3;-2,2] for grid (20,40,80,160)
 ************************************************************/
class Surface_Molecular : public Surface_Cartesian {
 private:
  virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;

 public:
  Surface_Molecular();
  virtual double Set_Node(Doub_I, Doub_I, Doub_I) const;
  virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
  virtual void display();
  virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif
