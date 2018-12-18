#pragma once

#include "constant.h"
#include "diffusion/beta.h"
#include "equation/equation.h"

class Eq_1 : public Equation {
private:
  Beta &beta;

  virtual double Inner_dux(Doub_I, Doub_I, Doub_I) const;
  virtual double Outer_dux(Doub_I, Doub_I, Doub_I) const;
  virtual double Inner_duy(Doub_I, Doub_I, Doub_I) const;
  virtual double Outer_duy(Doub_I, Doub_I, Doub_I) const;
  virtual double Inner_duz(Doub_I, Doub_I, Doub_I) const;
  virtual double Outer_duz(Doub_I, Doub_I, Doub_I) const;

public:
  Eq_1(Doub_I, Beta &);

  virtual double Outer_u(Doub_I, Doub_I, Doub_I) const;  // Outer source term
  virtual double Inner_u(Doub_I, Doub_I, Doub_I) const;  // Inner source term
  virtual double Outer_f(Doub_I, Doub_I, Doub_I) const;  // Outer source term
  virtual double Inner_f(Doub_I, Doub_I, Doub_I) const;  // Inner source term

  virtual double Jump_betau_x(Doub_I, Doub_I, Doub_I) const;
  virtual double Jump_betau_y(Doub_I, Doub_I, Doub_I) const;
  virtual double Jump_betau_z(Doub_I, Doub_I, Doub_I) const;
};
