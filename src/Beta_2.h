#ifndef __Beta_2_H__
#define __Beta_2_H__

#include "Beta.h"
#include "Constant.h"

class Beta_2 : public Beta {
 public:
  Beta_2();

  virtual double Inside(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dx(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dx(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dy(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dy(Doub_I, Doub_I, Doub_I) const;
  virtual double Inside_Dz(Doub_I, Doub_I, Doub_I) const;
  virtual double Outside_Dz(Doub_I, Doub_I, Doub_I) const;
};

#endif
