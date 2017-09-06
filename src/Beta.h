#ifndef __Beta_H__
#define __Beta_H__

#include "Constant.h"

class Beta {
 public:
  virtual double Inside(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outside(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Inside_Dx(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outside_Dx(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Inside_Dy(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outside_Dy(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Inside_Dz(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outside_Dz(Doub_I, Doub_I, Doub_I) const = 0;
};

#endif
