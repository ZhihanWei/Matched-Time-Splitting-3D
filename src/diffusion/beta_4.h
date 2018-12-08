#pragma once

#include "diffusion/beta.h"
#include "constant.h"

class Beta_4 : public Beta {
private:
  double in_sigma_x = 3;
  double in_sigma_y = 3;
  double in_sigma_z = 3;
  double out_sigma_x = 4;
  double out_sigma_y = 4;
  double out_sigma_z = 4;

public:
  Beta_4();

  virtual double Inside(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dx(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dx(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dy(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dy(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dz(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dz(Doub_I, Doub_I, Doub_I) const;
};
