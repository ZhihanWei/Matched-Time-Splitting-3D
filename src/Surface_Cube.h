#ifndef __Surface_Cube_H__
#define __Surface_Cube_H__

#include <string>
#include "Constant.h"
#include "Surface_Cartesian.h"

class Surface_Cube : public Surface_Cartesian {
 private:
  double xl, xr, yl, yr, zl, zr;

  // Not defined
  virtual double func(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfx(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfy(const Doub_I, const Doub_I, const Doub_I) const;
  virtual double dfz(const Doub_I, const Doub_I, const Doub_I) const;

 public:
  Surface_Cube(VecDoub_I);
  virtual double Set_Node(Doub_I, Doub_I, Doub_I) const;
  virtual VecDoub normal(Doub_I, Doub_I, Doub_I);
  virtual void display();

  // Override the member funcitons defined in Surface_Cartesian
  virtual double Gamma_x(Doub_I, Doub_I, Doub_I, Doub_I);
  virtual double Gamma_y(Doub_I, Doub_I, Doub_I, Doub_I);
  virtual double Gamma_z(Doub_I, Doub_I, Doub_I, Doub_I);

  // Not defined
  virtual double check(Doub_I, Doub_I, Doub_I);
};

#endif
