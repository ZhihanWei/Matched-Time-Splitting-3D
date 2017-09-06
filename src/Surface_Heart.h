#ifndef __Surface_Heart_H__
#define __Surface_Heart_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

/************************************************************
 Heart Surface with parameter [1,9/4,1,-1,-1,-9/80]
 Suggested domain: [-4,4;-2,2;-4,4]
 ************************************************************/
class Surface_Heart : public Surface_Cartesian {
 private:
  double par[6];

  virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;

 public:
  Surface_Heart(VecDoub_I);
  virtual double Set_Node(const Doub_I, const Doub_I, const Doub_I) const;
  virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
  virtual void display();
  virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif
