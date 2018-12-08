#pragma once

#include "constant.h"

class Equation {
private:
  virtual double Inner_dux(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outer_dux(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Inner_duy(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outer_duy(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Inner_duz(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Outer_duz(Doub_I, Doub_I, Doub_I) const = 0;

protected:
  double t;

public:
  // Outer analytical solution
  virtual double Outer_u(Doub_I, Doub_I, Doub_I) const = 0;
  // Inner analytical solution
  virtual double Inner_u(Doub_I, Doub_I, Doub_I) const = 0;
  // Outer source term
  virtual double Outer_f(Doub_I, Doub_I, Doub_I) const = 0;
  // Inner source term
  virtual double Inner_f(Doub_I, Doub_I, Doub_I) const = 0;

  virtual double Jump_betau_x(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Jump_betau_y(Doub_I, Doub_I, Doub_I) const = 0;
  virtual double Jump_betau_z(Doub_I, Doub_I, Doub_I) const = 0;

  double Jump_u(Doub_I, Doub_I, Doub_I) const;
  double Jump_betau_xi(Doub_I, Doub_I, Doub_I, Doub_I, Doub_I, Doub_I) const;
  double Jump_u_eta(Doub_I, Doub_I, Doub_I, Doub_I, Doub_I, Doub_I) const;
  double Jump_u_tau(Doub_I, Doub_I, Doub_I, Doub_I, Doub_I, Doub_I) const;
  double Jump_betau_eta(Doub_I, Doub_I, Doub_I, Doub_I, Doub_I, Doub_I) const;
  double Jump_betau_tau(Doub_I, Doub_I, Doub_I, Doub_I, Doub_I, Doub_I) const;
};
