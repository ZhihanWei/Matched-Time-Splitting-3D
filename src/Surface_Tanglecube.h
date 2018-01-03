#ifndef __Surface_Tanglecube_H__
#define __Surface_Tanglecube_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/**************************************************************************
 Tanglecube Surface with parameter [a1, a2, b1, b2, c1, c2, d]
 
 Equation : a1 * x^4 + a2 * x^2 + b1 * y^4 + b2 * y^2 + c1 * z^4 + c2 * z^2 + d = 0
 
 Suggested domain: [-4,4;-4,4;-5,5]
 **************************************************************************/
class Surface_Tanglecube : public Surface_Cartesian {
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
