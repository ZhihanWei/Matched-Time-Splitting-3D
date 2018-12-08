#pragma once

#include "diffusion/beta.h"
#include "constant.h"

class Beta_1 : public Beta {
public:
  Beta_1();

  virtual double Inside(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dx(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dx(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dy(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dy(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dz(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dz(Doub_I, Doub_I, Doub_I) const;
};
