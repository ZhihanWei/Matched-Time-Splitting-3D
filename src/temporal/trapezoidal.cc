#include <iomanip>
#include <iostream>

#include "trapezoidal.h"
/******************************************************************************************
 Constructor

 INPUT
 mesh     : mesh object
 inter    : object of all intersections
 beta     : object of variable coefficients
 time     : vector of 3 double values represent beginning time, finishing time
 and time step
 accuracy : accuracy of scheme
 *******************************************************************************************/
TS::TS(Intersections &inter, Mesh &mesh, Beta &beta, VecDoub_I time) {
  dt = time[2];

  nx = mesh.nx;
  ny = mesh.ny;
  nz = mesh.nz;

  dx = mesh.dx;
  dy = mesh.dy;
  dz = mesh.dz;

  xi = &mesh.xi[0];
  yi = &mesh.yi[0];
  zi = &mesh.zi[0];
  indicator = &mesh.mesh_value[0];
}

/*********************************************************************************
 Trapezoidal Splitting solver at each time step

 INPUT
 eq      : equation object at current time step
 eq_half : equation object at middle time step
 eq_one  : euqation object at next time step
 inter   : object of all intersections
 beta    : object of variable coefficient

 OUTPUT
 uh : three-dimensional solution at current time step to next time step
 ********************************************************************************/
void TS::Solve_2nd(Equation &eq, Equation &eq_half, Equation &eq_one,
                   Intersections &inter, CubicDoub &uh, Beta &beta) {
  VecDoub ax, bx, cx, rx, utx;
  VecDoub ay, by, cy, ry, uty;
  VecDoub az, bz, cz, rz, utz;
  CubicDoub v1, v2, v3, v4, v5, v6, v7, v8, uh1;
  CubicDoub src;

  // Initialize size for solution
  v1.resize(nx);
  v2.resize(nx);
  v3.resize(nx);
  v4.resize(nx);
  v5.resize(nx);
  v6.resize(nx);
  v7.resize(nx);
  v8.resize(nx);
  src.resize(nx);
  for (int ix = 0; ix < nx; ix++) {
    v1[ix].resize(ny);
    v2[ix].resize(ny);
    v3[ix].resize(ny);
    v4[ix].resize(ny);
    v5[ix].resize(ny);
    v6[ix].resize(ny);
    v7[ix].resize(ny);
    v8[ix].resize(ny);
    src[ix].resize(ny);
    for (int iy = 0; iy < ny; iy++) {
      v1[ix][iy].resize(nz);
      v2[ix][iy].resize(nz);
      v3[ix][iy].resize(nz);
      v4[ix][iy].resize(nz);
      v5[ix][iy].resize(nz);
      v6[ix][iy].resize(nz);
      v7[ix][iy].resize(nz);
      v8[ix][iy].resize(nz);
      src[ix][iy].resize(nz);
    }
  }

  //-------------- First half of Trapezoidal Splitting method [n] -> [n+1/2]
  //--------------------
  inter.Refresh_Jump(eq_one, uh, beta);

  /*******************************************************************************************************
   First Step:
   v1^{n+1/2} = ( 1 + dt/2 * beta * D_xx ) v1^{n}   with v1^{n} = u^{n}
   ******************************************************************************************************/
  // Set up RHS
  D_xx_r(eq_one, inter, uh, v1, beta);

  // Set up boundary conditions for v1
  Set_bc(eq_half, v1);

  /*******************************************************************************************************
   Second Step:
   v2^{n+1/2} = ( 1 + dt/2 * beta * D_yy ) v2^{n}   with v2^{n} = v1^{n+1/2}
   ******************************************************************************************************/
  // Set up RHS
  D_yy_r(eq_one, inter, v1, v2, beta);

  // Set up boundary conditions for v2
  Set_bc(eq_half, v2);

  /*******************************************************************************************************
   Thrid Step:
   v3^{n+1/2} = ( 1 + dt/2 * beta * D_zz ) v3^{n}   with v3^{n} = v2^{n+1/2}
   ******************************************************************************************************/
  // Set up RHS
  D_zz_r(eq_one, inter, v2, v3, beta);

  // Set up boundary conditions for v3
  Set_bc(eq_half, v3);

  /*******************************************************************************************************
   Fourth Step:
   v4^{n+1/2} = v4^{n} + dt/2 * f(t^{n}, v4^{n})  with v4^{n} = v3^{n+1/2}
   ******************************************************************************************************/
  // Initialize source term at time step n
  Src_2nd(eq, inter, src);

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        v4[ix][iy][iz] = v3[ix][iy][iz] + 0.5 * dt * src[ix][iy][iz];
      }
    }
  }

  //------------ Second half of Trapezoidal Splitting method [n+1/2] -> [n+1]
  //---------------------
  // inter.Refresh_Jump(eq_one,v4,beta);

  /*******************************************************************************************************
   Fifth Step:
   v5^{n+1} = v5^{n+1/2} + dt/2 * f(t^{n+1}, v5^{n+1})  with v5^{n} = v4^{n+1/2}
   ******************************************************************************************************/
  Src_2nd(eq_one, inter, src);

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        v5[ix][iy][iz] = v4[ix][iy][iz] + 0.5 * dt * src[ix][iy][iz];
      }
    }
  }

  /*******************************************************************************************************
   Sixth Step:
   (1 - dt/2 * beat * D_zz) v6^{n+1} = v6^{n+1/2}  with v6^{n+1/2} = v5^{n+1}
   ******************************************************************************************************/
  // Set up RHS
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        v6[ix][iy][iz] = v5[ix][iy][iz];
      }
    }
  }

  // Set up boundary conditions for v6
  Set_bc(eq_one, v6);

  // Set up LHS
  az.resize(nz);
  bz.resize(nz);
  cz.resize(nz);
  rz.resize(nz);
  utz.resize(nz);

  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      D_zz_l(ix, iy, eq_one, inter, v6, beta, az, bz, cz, rz);

      // Thomas Algorithm
      TDMA(az, bz, cz, rz, utz);

      for (int iz = 0; iz < nz; iz++) {
        v6[ix][iy][iz] = utz[iz];
      }
    }
  }

  /*******************************************************************************************************
   Seventh Step:
   (1 - dt/2 * beat * D_yy) v7^{n+1} = v7^{n+1/2}  with v7^{n+1/2} = v6^{n+1}
   ******************************************************************************************************/
  // Set up RHS
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        v7[ix][iy][iz] = v6[ix][iy][iz];
      }
    }
  }

  // Set up boundary conditions for v7
  Set_bc(eq_one, v7);

  // Set up LHS
  ay.resize(ny);
  by.resize(ny);
  cy.resize(ny);
  ry.resize(ny);
  uty.resize(ny);

  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      D_yy_l(ix, iz, eq_one, inter, v7, beta, ay, by, cy, ry);

      // Thomas Algorithm
      TDMA(ay, by, cy, ry, uty);

      for (int iy = 0; iy < ny; iy++) {
        v7[ix][iy][iz] = uty[iy];
      }
    }
  }

  /*******************************************************************************************************
   Eighth Step:
   (1 - dt/2 * beat * D_xx) v8^{n+1} = v8^{n+1/2}  with v8^{n+1/2} = v7^{n+1}
   ******************************************************************************************************/
  // Set up RHS
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        v8[ix][iy][iz] = v7[ix][iy][iz];
      }
    }
  }

  // Set up boundary conditions for v8
  Set_bc(eq_one, v8);

  // Set up LHS
  ax.resize(nx);
  bx.resize(nx);
  cx.resize(nx);
  rx.resize(nx);
  utx.resize(nx);

  for (int iy = 1; iy < ny - 1; iy++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      D_xx_l(iy, iz, eq_one, inter, v8, beta, ax, bx, cx, rx);

      // Thomas Algorithm
      TDMA(ax, bx, cx, rx, utx);

      for (int ix = 0; ix < nx; ix++) {
        v8[ix][iy][iz] = utx[ix];
      }
    }
  }

  /*******************************************************************************************************
   Update UH to next time step,  [t^{n},t^{n+1}]
   ******************************************************************************************************/
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        uh[ix][iy][iz] = v8[ix][iy][iz];
      }
    }
  }
}

/*********************************************************************************
 Set up boundary condition

 INPUT
 eq    : equation object at current time step

 OUTPUT
 uc    : cubic matrix
 *********************************************************************************/
void TS::Set_bc(Equation &eq, CubicDoub &uc) {
  // Set boundary conditions for uc
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      uc[ix][iy][0] = eq.Outer_u(xi[ix], yi[iy], zi[0]);
      uc[ix][iy][nz - 1] = eq.Outer_u(xi[ix], yi[iy], zi[nz - 1]);
    }
  }
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      uc[ix][0][iz] = eq.Outer_u(xi[ix], yi[0], zi[iz]);
      uc[ix][ny - 1][iz] = eq.Outer_u(xi[ix], yi[ny - 1], zi[iz]);
    }
  }
  for (int iy = 0; iy < ny; iy++) {
    for (int iz = 0; iz < nz; iz++) {
      uc[0][iy][iz] = eq.Outer_u(xi[0], yi[iy], zi[iz]);
      uc[nx - 1][iy][iz] = eq.Outer_u(xi[nx - 1], yi[iy], zi[iz]);
    }
  }
}

/*********************************************************************************
 1 + 1/2*Dt*beta*D_{xx} operator for right hand side vector

 INPUT
 iy    : coordinate index on y-direction
 iz    : coordinate index on z-direction
 eq    : equation object at current time step
 inter : object of all intersections
 beta  : object of variable coefficient
 uc1   : right hand side solution at current time step

 OUTPUT
 uc2   : right hand side solution at current time step
 ********************************************************************************/
void TS::D_xx_r(Equation &eq, Intersections &inter, CubicDoub &uc1,
                CubicDoub &uc2, Beta &beta) {
  int ix;
  size_t ip;
  double coef, sum;
  double jump_u, jump_ul, jump_ur;
  double jump_betaux, jump_betauxl, jump_betauxr;
  VecDoub vdex;
  Intersection_Data data, datal, datar;

  coef = 0.5; // coefficient of Dt*beta*D_{xx}

  vdex.resize(3);

  // Step I: Initialize without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          sum += uc1[ix + i][iy][iz] * vdex[i + 1];
        }

        uc2[ix][iy][iz] = coef * dt * sum;
      }
    }
  }

  // Step II: Apply MIB to operators
  for (int iy = 1; iy < ny - 1; iy++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      ip = 0;
      while (ip < inter.ifpx[iy][iz].size()) {
        // An irregular interface point
        if (inter.ifpx[iy][iz][ip].ID > 0) {
          data = inter.ifpx[iy][iz][ip];

          jump_u = eq.Jump_u(data.coord.x_value, data.coord.y_value,
                             data.coord.z_value);
          if (JP == 'r') {
            jump_betaux = eq.Jump_betau_x(
                data.coord.x_value, data.coord.y_value, data.coord.z_value);
          } else {
            jump_betaux = inter.ifpx[iy][iz][ip].jump.u_dir;
          }

          // Approximate IX
          ix = data.left_loc;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdex[i + 1] * uc1[ix + i][iy][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdex[2] * uc1[data.left_loc - 1 + i][iy][iz] *
                   data.wei.weir[0][i];
          }
          sum += vdex[2] * (jump_u * data.wei.weir[0][4] +
                            jump_betaux * data.wei.weir[0][5]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IX+1
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdex[i + 1] * uc1[ix + i][iy][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdex[0] * uc1[data.left_loc - 1 + i][iy][iz] *
                   data.wei.weil[0][i];
          }
          sum += vdex[0] * (jump_u * data.wei.weil[0][4] +
                            jump_betaux * data.wei.weil[0][5]);
          uc2[ix][iy][iz] = coef * dt * sum;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpx[iy][iz][ip];
          datar = inter.ifpx[iy][iz][ip + 1];

          jump_ul = eq.Jump_u(datal.coord.x_value, datal.coord.y_value,
                              datal.coord.z_value);
          jump_ur = eq.Jump_u(datar.coord.x_value, datar.coord.y_value,
                              datar.coord.z_value);
          if (JP == 'r') {
            jump_betauxl = eq.Jump_betau_x(
                datal.coord.x_value, datal.coord.y_value, datal.coord.z_value);
            jump_betauxr = eq.Jump_betau_x(
                datar.coord.x_value, datar.coord.y_value, datar.coord.z_value);
          } else {
            jump_betauxl = inter.ifpx[iy][iz][ip].jump.u_dir;
            jump_betauxr = inter.ifpx[iy][iz][ip + 1].jump.u_dir;
          }

          // Approximate IX
          ix = datal.left_loc;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdex[i + 1] * uc1[ix + i][iy][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdex[2] * uc1[datal.left_loc - 1 + i][iy][iz] *
                   datal.wei.weir[0][i];
          }
          sum += vdex[2] * (jump_ul * datal.wei.weir[0][5] +
                            jump_betauxl * datal.wei.weir[0][6] +
                            jump_ur * datal.wei.weir[0][7] +
                            jump_betauxr * datal.wei.weir[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IX+1, the inside corner point
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          sum = 0;
          // IX-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdex[0] * uc1[datal.left_loc - 1 + i][iy][iz] *
                   datal.wei.weil[0][i];
          }
          sum += vdex[0] * (jump_ul * datal.wei.weil[0][5] +
                            jump_betauxl * datal.wei.weil[0][6] +
                            jump_ur * datal.wei.weil[0][7] +
                            jump_betauxr * datal.wei.weil[0][8]);
          // IX,  no FP
          sum += vdex[1] * uc1[ix][iy][iz];
          // IX+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdex[2] * uc1[datar.left_loc - 2 + i][iy][iz] *
                   datar.wei.weir[0][i];
          }
          sum += vdex[2] * (jump_ul * datar.wei.weir[0][5] +
                            jump_betauxl * datar.wei.weir[0][6] +
                            jump_ur * datar.wei.weir[0][7] +
                            jump_betauxr * datar.wei.weir[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IX+2, right FP
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdex[i + 1] * uc1[ix + i][iy][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdex[0] * uc1[datar.left_loc - 2 + i][iy][iz] *
                   datar.wei.weil[0][i];
          }
          sum += vdex[0] * (jump_ul * datar.wei.weil[0][5] +
                            jump_betauxl * datar.wei.weil[0][6] +
                            jump_ur * datar.wei.weil[0][7] +
                            jump_betauxr * datar.wei.weil[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          ip += 2;
        }
      }
    }
  }

  // Step III: add U_{i,j,k}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uc2[ix][iy][iz] += uc1[ix][iy][iz];
      }
    }
  }
}

/*********************************************************************************
 1 + 1/2*Dt*beta*D_{yy} operator for right hand side vector

 INPUT
 iy    : coordinate index on y-direction
 iz    : coordinate index on z-direction
 eq    : equation object at current time step
 inter : object of all intersections
 beta  : object of variable coefficient
 uc1   : right hand side solution at current time step

 OUTPUT
 uc2   : right hand side solution at current time step
 ********************************************************************************/
void TS::D_yy_r(Equation &eq, Intersections &inter, CubicDoub &uc1,
                CubicDoub &uc2, Beta &beta) {
  int iy;
  size_t ip;
  double coef, sum;
  double jump_u, jump_ul, jump_ur;
  double jump_betauy, jump_betauyl, jump_betauyr;
  VecDoub vdey;
  Intersection_Data data, datal, datar;

  coef = 0.5; // coefficient of Dt*beta*D_{yy}

  vdey.resize(3);

  // Step I: Initialize without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          sum += uc1[ix][iy + i][iz] * vdey[i + 1];
        }

        uc2[ix][iy][iz] = coef * dt * sum;
      }
    }
  }

  // Step II: Apply MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      ip = 0;
      while (ip < inter.ifpy[ix][iz].size()) {
        // An irregular interface point
        if (inter.ifpy[ix][iz][ip].ID > 0) {
          data = inter.ifpy[ix][iz][ip];

          jump_u = eq.Jump_u(data.coord.x_value, data.coord.y_value,
                             data.coord.z_value);
          if (JP == 'r') {
            jump_betauy = eq.Jump_betau_y(
                data.coord.x_value, data.coord.y_value, data.coord.z_value);
          } else {
            jump_betauy = inter.ifpy[ix][iz][ip].jump.u_dir;
          }

          // Approximate IY
          iy = data.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdey[i + 1] * uc1[ix][iy + i][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdey[2] * uc1[ix][data.left_loc - 1 + i][iz] *
                   data.wei.weir[0][i];
          }
          sum += vdey[2] * (jump_u * data.wei.weir[0][4] +
                            jump_betauy * data.wei.weir[0][5]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IY+1
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdey[i + 1] * uc1[ix][iy + i][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdey[0] * uc1[ix][data.left_loc - 1 + i][iz] *
                   data.wei.weil[0][i];
          }
          sum += vdey[0] * (jump_u * data.wei.weil[0][4] +
                            jump_betauy * data.wei.weil[0][5]);
          uc2[ix][iy][iz] = coef * dt * sum;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpy[ix][iz][ip];
          datar = inter.ifpy[ix][iz][ip + 1];

          jump_ul = eq.Jump_u(datal.coord.x_value, datal.coord.y_value,
                              datal.coord.z_value);
          jump_ur = eq.Jump_u(datar.coord.x_value, datar.coord.y_value,
                              datar.coord.z_value);
          if (JP == 'r') {
            jump_betauyl = eq.Jump_betau_y(
                datal.coord.x_value, datal.coord.y_value, datal.coord.z_value);
            jump_betauyr = eq.Jump_betau_y(
                datar.coord.x_value, datar.coord.y_value, datar.coord.z_value);
          } else {
            jump_betauyl = inter.ifpy[ix][iz][ip].jump.u_dir;
            jump_betauyr = inter.ifpy[ix][iz][ip + 1].jump.u_dir;
          }

          // Approximate IY
          iy = datal.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdey[i + 1] * uc1[ix][iy + i][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdey[2] * uc1[ix][datal.left_loc - 1 + i][iz] *
                   datal.wei.weir[0][i];
          }
          sum += vdey[2] * (jump_ul * datal.wei.weir[0][5] +
                            jump_betauyl * datal.wei.weir[0][6] +
                            jump_ur * datal.wei.weir[0][7] +
                            jump_betauyr * datal.wei.weir[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IY+1, the inside corner point
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          sum = 0;
          // IY-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdey[0] * uc1[ix][datal.left_loc - 1 + i][iz] *
                   datal.wei.weil[0][i];
          }
          sum += vdey[0] * (jump_ul * datal.wei.weil[0][5] +
                            jump_betauyl * datal.wei.weil[0][6] +
                            jump_ur * datal.wei.weil[0][7] +
                            jump_betauyr * datal.wei.weil[0][8]);
          // IY,  no FP
          sum += vdey[1] * uc1[ix][iy][iz];
          // IY+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdey[2] * uc1[ix][datar.left_loc - 2 + i][iz] *
                   datar.wei.weir[0][i];
          }
          sum += vdey[2] * (jump_ul * datar.wei.weir[0][5] +
                            jump_betauyl * datar.wei.weir[0][6] +
                            jump_ur * datar.wei.weir[0][7] +
                            jump_betauyr * datar.wei.weir[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IY+2, right FP
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdey[i + 1] * uc1[ix][iy + i][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdey[0] * uc1[ix][datar.left_loc - 2 + i][iz] *
                   datar.wei.weil[0][i];
          }
          sum += vdey[0] * (jump_ul * datar.wei.weil[0][5] +
                            jump_betauyl * datar.wei.weil[0][6] +
                            jump_ur * datar.wei.weil[0][7] +
                            jump_betauyr * datar.wei.weil[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          ip += 2;
        }
      }
    }
  }

  // Step III: add U_{i,j,k}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uc2[ix][iy][iz] += uc1[ix][iy][iz];
      }
    }
  }
}

/*********************************************************************************
 1 + 1/2*Dt*beta*D_{zz} operator for right hand side vector

 INPUT
 iy    : coordinate index on y-direction
 iz    : coordinate index on z-direction
 eq    : equation object at current time step
 inter : object of all intersections
 beta  : object of variable coefficient
 uc1   : right hand side solution at current time step

 OUTPUT
 uc2   : right hand side solution at current time step
 ********************************************************************************/
void TS::D_zz_r(Equation &eq, Intersections &inter, CubicDoub &uc1,
                CubicDoub &uc2, Beta &beta) {
  int iz;
  size_t ip;
  double coef, sum;
  double jump_u, jump_ul, jump_ur;
  double jump_betauz, jump_betauzl, jump_betauzr;
  VecDoub vdez;
  Intersection_Data data, datal, datar;

  coef = 0.5; // coefficient of Dt*beta*D_{zz}

  vdez.resize(3);

  // Step I: Initialize without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          sum += uc1[ix][iy][iz + i] * vdez[i + 1];
        }
        uc2[ix][iy][iz] = coef * dt * sum;
      }
    }
  }

  // Step II-2: Apply MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      ip = 0;
      while (ip < inter.ifpz[ix][iy].size()) {
        // An irregular interface point
        if (inter.ifpz[ix][iy][ip].ID > 0) {
          data = inter.ifpz[ix][iy][ip];

          jump_u = eq.Jump_u(data.coord.x_value, data.coord.y_value,
                             data.coord.z_value);
          if (JP == 'r') {
            jump_betauz = eq.Jump_betau_z(
                data.coord.x_value, data.coord.y_value, data.coord.z_value);
          } else {
            jump_betauz = inter.ifpz[ix][iy][ip].jump.u_dir;
          }

          // Approximate IZ
          iz = data.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdez[i + 1] * uc1[ix][iy][iz + i];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdez[2] * uc1[ix][iy][data.left_loc - 1 + i] *
                   data.wei.weir[0][i];
          }
          sum += vdez[2] * (jump_u * data.wei.weir[0][4] +
                            jump_betauz * data.wei.weir[0][5]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IZ+1
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdez[i + 1] * uc1[ix][iy][iz + i];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdez[0] * uc1[ix][iy][data.left_loc - 1 + i] *
                   data.wei.weil[0][i];
          }
          sum += vdez[0] * (jump_u * data.wei.weil[0][4] +
                            jump_betauz * data.wei.weil[0][5]);
          uc2[ix][iy][iz] = coef * dt * sum;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpz[ix][iy][ip];
          datar = inter.ifpz[ix][iy][ip + 1];

          jump_ul = eq.Jump_u(datal.coord.x_value, datal.coord.y_value,
                              datal.coord.z_value);
          jump_ur = eq.Jump_u(datar.coord.x_value, datar.coord.y_value,
                              datar.coord.z_value);
          if (JP == 'r') {
            jump_betauzl = eq.Jump_betau_z(
                datal.coord.x_value, datal.coord.y_value, datal.coord.z_value);
            jump_betauzr = eq.Jump_betau_z(
                datar.coord.x_value, datar.coord.y_value, datar.coord.z_value);
          } else {
            jump_betauzl = inter.ifpz[ix][iy][ip].jump.u_dir;
            jump_betauzr = inter.ifpz[ix][iy][ip + 1].jump.u_dir;
          }

          // Approximate IZ, left FP
          iz = datal.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdez[i + 1] * uc1[ix][iy][iz + i];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdez[2] * uc1[ix][iy][datal.left_loc - 1 + i] *
                   datal.wei.weir[0][i];
          }
          sum += vdez[2] * (jump_ul * datal.wei.weir[0][5] +
                            jump_betauzl * datal.wei.weir[0][6] +
                            jump_ur * datal.wei.weir[0][7] +
                            jump_betauzr * datal.wei.weir[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IZ+1, the inside corner point
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          sum = 0;
          // IZ-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdez[0] * uc1[ix][iy][datal.left_loc - 1 + i] *
                   datal.wei.weil[0][i];
          }
          sum += vdez[0] * (jump_ul * datal.wei.weil[0][5] +
                            jump_betauzl * datal.wei.weil[0][6] +
                            jump_ur * datal.wei.weil[0][7] +
                            jump_betauzr * datal.wei.weil[0][8]);
          // IZ,  no FP
          sum += vdez[1] * uc1[ix][iy][iz];
          // IZ+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdez[2] * uc1[ix][iy][datar.left_loc - 2 + i] *
                   datar.wei.weir[0][i];
          }
          sum += vdez[2] * (jump_ul * datar.wei.weir[0][5] +
                            jump_betauzl * datar.wei.weir[0][6] +
                            jump_ur * datar.wei.weir[0][7] +
                            jump_betauzr * datar.wei.weir[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          // Approximate IZ+2, right FP
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdez[i + 1] * uc1[ix][iy][iz + i];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdez[0] * uc1[ix][iy][datar.left_loc - 2 + i] *
                   datar.wei.weil[0][i];
          }
          sum += vdez[0] * (jump_ul * datar.wei.weil[0][5] +
                            jump_betauzl * datar.wei.weil[0][6] +
                            jump_ur * datar.wei.weil[0][7] +
                            jump_betauzr * datar.wei.weil[0][8]);
          uc2[ix][iy][iz] = coef * dt * sum;

          ip += 2;
        }
      }
    }
  }

  // Step III: add U_{i,j,k}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uc2[ix][iy][iz] += uc1[ix][iy][iz];
      }
    }
  }
}

/*********************************************************************************
 1 - 1/2*Dt*beta*D_{xx} operator for left hand side matrix

 INPUT
 iy    : coordinate index on y-direction
 iz    : coordinate index on z-direction
 eq    : equation object at current time step
 inter : object of all intersections
 uc    : right hand side solution at current time step
 beta  : object of variable coefficient

 OUTPUT
 a     : first hypotenuse of tridiagonal matrix
 b     : second hypotenuse of tridiagonal matrix
 c     : third hypotenuse of tridiagonal matrix
 r     : right hand side vector
 ********************************************************************************/
void TS::D_xx_l(Int_I iy, Int_I iz, Equation &eq, Intersections &inter,
                CubicDoub &uc, Beta &beta, VecDoub_O &a, VecDoub_O &b,
                VecDoub_O &c, VecDoub_O &r) {
  int ix;
  size_t ip;
  double coef;
  double jump_u, jump_ul, jump_ur;
  double jump_betaux, jump_betauxl, jump_betauxr;
  VecDoub vdex;
  MatrixDoub irr_row, cor_row;
  Intersection_Data data, datal, datar;

  coef = 0.5; // coefficient of Dt*beta*D_{xx}

  vdex.resize(3);
  // set up irr_row and cor_row
  irr_row.resize(2);
  for (int i = 0; i < 2; i++) {
    irr_row[i].resize(5);
  }
  cor_row.resize(3);
  for (int i = 0; i < 3; i++) {
    cor_row[i].resize(6);
  }

  // Step I: LHS without MIB
  a[0] = 0.0;
  b[0] = 1.0;
  c[0] = 0.0;
  for (int ix = 1; ix < nx - 1; ix++) {
    Operator_weights_x(beta, vdex, ix, iy, iz, dx);

    a[ix] = -coef * dt * vdex[0];
    b[ix] = 1 - coef * dt * vdex[1];
    c[ix] = -coef * dt * vdex[2];
  }
  a[nx - 1] = 0.0;
  b[nx - 1] = 1.0;
  c[nx - 1] = 0.0;

  for (int ix = 0; ix < nx; ix++) {
    r[ix] = uc[ix][iy][iz];
  }

  // Step II: LHS with MIB
  ip = 0;
  while (ip < inter.ifpx[iy][iz].size()) {
    // An irregular interface
    if (inter.ifpx[iy][iz][ip].ID > 0) {
      data = inter.ifpx[iy][iz][ip];

      jump_u =
          eq.Jump_u(data.coord.x_value, data.coord.y_value, data.coord.z_value);
      if (JP == 'r') {
        jump_betaux = eq.Jump_betau_x(data.coord.x_value, data.coord.y_value,
                                      data.coord.z_value);
      } else {
        jump_betaux = inter.ifpx[iy][iz][ip].jump.u_dir;
      }

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 5; j++) {
          irr_row[i][j] = 0;
        }
      }

      // Approximate IX
      ix = data.left_loc;
      Operator_weights_x(beta, vdex, ix, iy, iz, dx);
      // use right FP
      irr_row[0][0] = -coef * dt * vdex[0];
      irr_row[0][1] = 1 - coef * dt * vdex[1];
      for (int i = 0; i < 4; i++) {
        irr_row[0][i] += -coef * dt * vdex[2] * data.wei.weir[0][i];
      }
      irr_row[0][4] = uc[ix][iy][iz] +
                      coef * dt * vdex[2] * (jump_u * data.wei.weir[0][4] +
                                             jump_betaux * data.wei.weir[0][5]);

      // Approximate IX+1
      ix += 1;
      Operator_weights_x(beta, vdex, ix, iy, iz, dx);
      // use left FP
      irr_row[1][2] = 1 - coef * dt * vdex[1];
      irr_row[1][3] = -coef * dt * vdex[2];
      for (int i = 0; i < 4; i++) {
        irr_row[1][i] += -coef * dt * vdex[0] * data.wei.weil[0][i];
      }
      irr_row[1][4] = uc[ix][iy][iz] +
                      coef * dt * vdex[0] * (jump_u * data.wei.weil[0][4] +
                                             jump_betaux * data.wei.weil[0][5]);

      Convert2Tri_irr(data.left_loc, irr_row, a, b, c, r);

      ip += 1;
    }

    // Corner interface
    else {
      datal = inter.ifpx[iy][iz][ip];
      datar = inter.ifpx[iy][iz][ip + 1];

      jump_ul = eq.Jump_u(datal.coord.x_value, datal.coord.y_value,
                          datal.coord.z_value);
      jump_ur = eq.Jump_u(datar.coord.x_value, datar.coord.y_value,
                          datar.coord.z_value);
      if (JP == 'r') {
        jump_betauxl = eq.Jump_betau_x(datal.coord.x_value, datal.coord.y_value,
                                       datal.coord.z_value);
        jump_betauxr = eq.Jump_betau_x(datar.coord.x_value, datar.coord.y_value,
                                       datar.coord.z_value);
      } else {
        jump_betauxl = inter.ifpx[iy][iz][ip].jump.u_dir;
        jump_betauxr = inter.ifpx[iy][iz][ip + 1].jump.u_dir;
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
          cor_row[i][j] = 0;
        }
      }

      // Approximate IX
      ix = datal.left_loc;
      Operator_weights_x(beta, vdex, ix, iy, iz, dx);
      // IX+1 cross left interface, use left interface's left FP
      cor_row[0][0] = -coef * dt * vdex[0];
      cor_row[0][1] = 1 - coef * dt * vdex[1];
      for (int i = 0; i < 5; i++) {
        cor_row[0][i] -= coef * dt * vdex[2] * datal.wei.weir[0][i];
      }
      cor_row[0][5] =
          uc[ix][iy][iz] +
          coef * dt * vdex[2] * (jump_ul * datal.wei.weir[0][5] +
                                 jump_betauxl * datal.wei.weir[0][6] +
                                 jump_ur * datal.wei.weir[0][7] +
                                 jump_betauxr * datal.wei.weir[0][8]);

      // Approximate IX+1
      ix += 1;
      Operator_weights_x(beta, vdex, ix, iy, iz, dx);
      cor_row[1][2] = 1 - coef * dt * vdex[1];
      // IX-1 cross left interface, use left interface's left FP
      for (int i = 0; i < 5; i++) {
        cor_row[1][i] -= coef * dt * vdex[0] * datal.wei.weil[0][i];
      }
      cor_row[1][5] =
          uc[ix][iy][iz] +
          coef * dt * vdex[0] * (jump_ul * datal.wei.weil[0][5] +
                                 jump_betauxl * datal.wei.weil[0][6] +
                                 jump_ur * datal.wei.weil[0][7] +
                                 jump_betauxr * datal.wei.weil[0][8]);
      // IX+1 cross right interface, use right interface's right FP
      for (int i = 0; i < 5; i++) {
        cor_row[1][i] -= coef * dt * vdex[2] * datar.wei.weir[0][i];
      }
      cor_row[1][5] +=
          coef * dt * vdex[2] * (jump_ul * datar.wei.weir[0][5] +
                                 jump_betauxl * datar.wei.weir[0][6] +
                                 jump_ur * datar.wei.weir[0][7] +
                                 jump_betauxr * datar.wei.weir[0][8]);

      // Approximate IX+2
      ix += 1;
      Operator_weights_x(beta, vdex, ix, iy, iz, dx);
      // IX-1 cross right interface, use right interface's left FP
      cor_row[2][3] = 1 - coef * dt * vdex[1];
      cor_row[2][4] = -coef * dt * vdex[2];
      for (int i = 0; i < 5; i++) {
        cor_row[2][i] -= coef * dt * vdex[0] * datar.wei.weil[0][i];
      }
      cor_row[2][5] =
          uc[ix][iy][iz] +
          coef * dt * vdex[0] * (jump_ul * datar.wei.weil[0][5] +
                                 jump_betauxl * datar.wei.weil[0][6] +
                                 jump_ur * datar.wei.weil[0][7] +
                                 jump_betauxr * datar.wei.weil[0][8]);

      Convert2Tri_cor(datal.left_loc, cor_row, a, b, c, r);

      ip += 2;
    }
  }
}

/*********************************************************************************
 1 - 1/2*Dt*beta*D_{yy} operator for left hand side matrix

 INPUT
 iy    : coordinate index on y-direction
 iz    : coordinate index on z-direction
 eq    : equation object at current time step
 inter : object of all intersections
 uc    : right han side solution at current time step
 beta  : object of variable coefficient

 OUTPUT
 a     : first hypotenuse of tridiagonal matrix
 b     : second hypotenuse of tridiagonal matrix
 c     : third hypotenuse of tridiagonal matrix
 r     : right hand side vector
 ********************************************************************************/
void TS::D_yy_l(Int_I ix, Int_I iz, Equation &eq, Intersections &inter,
                CubicDoub &uc, Beta &beta, VecDoub_O &a, VecDoub_O &b,
                VecDoub_O &c, VecDoub_O &r) {
  int iy;
  size_t ip;
  double coef;
  double jump_u, jump_ul, jump_ur;
  double jump_betauy, jump_betauyl, jump_betauyr;
  VecDoub vdey;
  MatrixDoub irr_row, cor_row;
  Intersection_Data data, datal, datar;

  coef = 0.5; // coefficient of Dt*beta*D_{yy}

  vdey.resize(3);
  // set up irr_row and cor_row
  irr_row.resize(2);
  for (int i = 0; i < 2; i++) {
    irr_row[i].resize(5);
  }
  cor_row.resize(3);
  for (int i = 0; i < 3; i++) {
    cor_row[i].resize(6);
  }

  // Step I: LHS without MIB
  a[0] = 0.0;
  b[0] = 1.0;
  c[0] = 0.0;
  for (int iy = 1; iy < ny - 1; iy++) {
    Operator_weights_y(beta, vdey, ix, iy, iz, dy);

    a[iy] = -coef * dt * vdey[0];
    b[iy] = 1 - coef * dt * vdey[1];
    c[iy] = -coef * dt * vdey[2];
  }
  a[ny - 1] = 0.0;
  b[ny - 1] = 1.0;
  c[ny - 1] = 0.0;

  for (int iy = 0; iy < ny; iy++) {
    r[iy] = uc[ix][iy][iz];
  }

  // Step II: LHS with MIB
  ip = 0;
  while (ip < inter.ifpy[ix][iz].size()) {
    // An irregular interface
    if (inter.ifpy[ix][iz][ip].ID > 0) {
      data = inter.ifpy[ix][iz][ip];

      jump_u =
          eq.Jump_u(data.coord.x_value, data.coord.y_value, data.coord.z_value);
      if (JP == 'r') {
        jump_betauy = eq.Jump_betau_y(data.coord.x_value, data.coord.y_value,
                                      data.coord.z_value);
      } else {
        jump_betauy = inter.ifpy[ix][iz][ip].jump.u_dir;
      }

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 5; j++) {
          irr_row[i][j] = 0;
        }
      }

      // Approximate IY
      iy = data.left_loc;
      Operator_weights_y(beta, vdey, ix, iy, iz, dy);
      // use right FP
      irr_row[0][0] = -coef * dt * vdey[0];
      irr_row[0][1] = 1 - coef * dt * vdey[1];
      for (int i = 0; i < 4; i++) {
        irr_row[0][i] += -coef * dt * vdey[2] * data.wei.weir[0][i];
      }
      irr_row[0][4] = uc[ix][iy][iz] +
                      coef * dt * vdey[2] * (jump_u * data.wei.weir[0][4] +
                                             jump_betauy * data.wei.weir[0][5]);

      // Approximate IY+1
      iy += 1;
      Operator_weights_y(beta, vdey, ix, iy, iz, dy);
      // use left FP
      irr_row[1][2] = 1 - coef * dt * vdey[1];
      irr_row[1][3] = -coef * dt * vdey[2];
      for (int i = 0; i < 4; i++) {
        irr_row[1][i] += -coef * dt * vdey[0] * data.wei.weil[0][i];
      }
      irr_row[1][4] = uc[ix][iy][iz] +
                      coef * dt * vdey[0] * (jump_u * data.wei.weil[0][4] +
                                             jump_betauy * data.wei.weil[0][5]);

      Convert2Tri_irr(data.left_loc, irr_row, a, b, c, r);

      ip += 1;
    }
    // Corner interface
    else {
      datal = inter.ifpy[ix][iz][ip];
      datar = inter.ifpy[ix][iz][ip + 1];

      jump_ul = eq.Jump_u(datal.coord.x_value, datal.coord.y_value,
                          datal.coord.z_value);
      jump_ur = eq.Jump_u(datar.coord.x_value, datar.coord.y_value,
                          datar.coord.z_value);
      if (JP == 'r') {
        jump_betauyl = eq.Jump_betau_y(datal.coord.x_value, datal.coord.y_value,
                                       datal.coord.z_value);
        jump_betauyr = eq.Jump_betau_y(datar.coord.x_value, datar.coord.y_value,
                                       datar.coord.z_value);
      } else {
        jump_betauyl = inter.ifpy[ix][iz][ip].jump.u_dir;
        jump_betauyr = inter.ifpy[ix][iz][ip + 1].jump.u_dir;
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
          cor_row[i][j] = 0;
        }
      }

      // Approximate IY
      iy = datal.left_loc;
      Operator_weights_y(beta, vdey, ix, iy, iz, dy);
      // IY+1 cross left interface, use left interface's right FP
      cor_row[0][0] = -coef * dt * vdey[0];
      cor_row[0][1] = 1 - coef * dt * vdey[1];
      for (int i = 0; i < 5; i++) {
        cor_row[0][i] -= coef * dt * vdey[2] * datal.wei.weir[0][i];
      }
      cor_row[0][5] =
          uc[ix][iy][iz] +
          coef * dt * vdey[2] * (jump_ul * datal.wei.weir[0][5] +
                                 jump_betauyl * datal.wei.weir[0][6] +
                                 jump_ur * datal.wei.weir[0][7] +
                                 jump_betauyr * datal.wei.weir[0][8]);

      // Approximate IY+1
      iy += 1;
      Operator_weights_y(beta, vdey, ix, iy, iz, dy);
      cor_row[1][2] = 1 - coef * dt * vdey[1];
      // IY-1 cross left interface, use left interface's left FP
      for (int i = 0; i < 5; i++) {
        cor_row[1][i] -= coef * dt * vdey[0] * datal.wei.weil[0][i];
      }
      cor_row[1][5] =
          uc[ix][iy][iz] +
          coef * dt * vdey[0] * (jump_ul * datal.wei.weil[0][5] +
                                 jump_betauyl * datal.wei.weil[0][6] +
                                 jump_ur * datal.wei.weil[0][7] +
                                 jump_betauyr * datal.wei.weil[0][8]);
      // IY+1 cross right interface, use right interface's right FP
      for (int i = 0; i < 5; i++) {
        cor_row[1][i] -= coef * dt * vdey[2] * datar.wei.weir[0][i];
      }
      cor_row[1][5] +=
          coef * dt * vdey[2] * (jump_ul * datar.wei.weir[0][5] +
                                 jump_betauyl * datar.wei.weir[0][6] +
                                 jump_ur * datar.wei.weir[0][7] +
                                 jump_betauyr * datar.wei.weir[0][8]);

      // Approximate IY+2
      iy += 1;
      Operator_weights_y(beta, vdey, ix, iy, iz, dy);
      // IY-1 cross right interface, use right interface's left FP
      cor_row[2][3] = 1 - coef * dt * vdey[1];
      cor_row[2][4] = -coef * dt * vdey[2];
      for (int i = 0; i < 5; i++) {
        cor_row[2][i] -= coef * dt * vdey[0] * datar.wei.weil[0][i];
      }
      cor_row[2][5] =
          uc[ix][iy][iz] +
          coef * dt * vdey[0] * (jump_ul * datar.wei.weil[0][5] +
                                 jump_betauyl * datar.wei.weil[0][6] +
                                 jump_ur * datar.wei.weil[0][7] +
                                 jump_betauyr * datar.wei.weil[0][8]);

      Convert2Tri_cor(datal.left_loc, cor_row, a, b, c, r);

      ip += 2;
    }
  }
}

/*********************************************************************************
 1 - 1/2*Dt*beta*D_{zz} operator for left hand side matrix

 INPUT
 iy    : coordinate index on y-direction
 iz    : coordinate index on z-direction
 eq    : equation object at current time step
 inter : object of all intersections
 uc    : right han side solution at current time step
 beta  : object of variable coefficient

 OUTPUT
 a     : first hypotenuse of tridiagonal matrix
 b     : second hypotenuse of tridiagonal matrix
 c     : third hypotenuse of tridiagonal matrix
 r     : right hand side vector
 ********************************************************************************/
void TS::D_zz_l(Int_I ix, Int_I iy, Equation &eq, Intersections &inter,
                CubicDoub &uc, Beta &beta, VecDoub_O &a, VecDoub_O &b,
                VecDoub_O &c, VecDoub_O &r) {
  int iz;
  size_t ip;
  double coef;
  double jump_u, jump_ul, jump_ur;
  double jump_betauz, jump_betauzl, jump_betauzr;
  VecDoub vdez;
  MatrixDoub irr_row, cor_row;
  Intersection_Data data, datal, datar;

  coef = 0.5; // coefficient of Dt*beta*D_{zz}

  vdez.resize(3);
  // set up irr_row and cor_row
  irr_row.resize(2);
  for (int i = 0; i < 2; i++) {
    irr_row[i].resize(5);
  }
  cor_row.resize(3);
  for (int i = 0; i < 3; i++) {
    cor_row[i].resize(6);
  }

  // Step I: LHS without MIB
  a[0] = 0.0;
  b[0] = 1.0;
  c[0] = 0.0;
  for (int iz = 1; iz < nz - 1; iz++) {
    Operator_weights_z(beta, vdez, ix, iy, iz, dz);

    a[iz] = -coef * dt * vdez[0];
    b[iz] = 1 - coef * dt * vdez[1];
    c[iz] = -coef * dt * vdez[2];
  }
  a[nz - 1] = 0.0;
  b[nz - 1] = 1.0;
  c[nz - 1] = 0.0;

  for (int iz = 0; iz < nz; iz++) {
    r[iz] = uc[ix][iy][iz];
  }

  // Step II: LHS with MIB
  ip = 0;
  while (ip < inter.ifpz[ix][iy].size()) {
    // An irregular interface
    if (inter.ifpz[ix][iy][ip].ID > 0) {
      data = inter.ifpz[ix][iy][ip];

      jump_u =
          eq.Jump_u(data.coord.x_value, data.coord.y_value, data.coord.z_value);
      if (JP == 'r') {
        jump_betauz = eq.Jump_betau_z(data.coord.x_value, data.coord.y_value,
                                      data.coord.z_value);
      } else {
        jump_betauz = inter.ifpz[ix][iy][ip].jump.u_dir;
      }

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 5; j++) {
          irr_row[i][j] = 0;
        }
      }

      // Approximate IZ
      iz = data.left_loc;
      Operator_weights_z(beta, vdez, ix, iy, iz, dz);
      // use right FP
      irr_row[0][0] = -coef * dt * vdez[0];
      irr_row[0][1] = 1 - coef * dt * vdez[1];
      for (int i = 0; i < 4; i++) {
        irr_row[0][i] += -coef * dt * vdez[2] * data.wei.weir[0][i];
      }
      irr_row[0][4] = uc[ix][iy][iz] +
                      coef * dt * vdez[2] * (jump_u * data.wei.weir[0][4] +
                                             jump_betauz * data.wei.weir[0][5]);

      // Approximate IZ+1
      iz += 1;
      Operator_weights_z(beta, vdez, ix, iy, iz, dz);
      // use left FP
      irr_row[1][2] = 1 - coef * dt * vdez[1];
      irr_row[1][3] = -coef * dt * vdez[2];
      for (int i = 0; i < 4; i++) {
        irr_row[1][i] += -coef * dt * vdez[0] * data.wei.weil[0][i];
      }
      irr_row[1][4] = uc[ix][iy][iz] +
                      coef * dt * vdez[0] * (jump_u * data.wei.weil[0][4] +
                                             jump_betauz * data.wei.weil[0][5]);

      Convert2Tri_irr(data.left_loc, irr_row, a, b, c, r);

      ip += 1;
    }
    // Corner interface
    else {
      datal = inter.ifpz[ix][iy][ip];
      datar = inter.ifpz[ix][iy][ip + 1];

      jump_ul = eq.Jump_u(datal.coord.x_value, datal.coord.y_value,
                          datal.coord.z_value);
      jump_ur = eq.Jump_u(datar.coord.x_value, datar.coord.y_value,
                          datar.coord.z_value);
      if (JP == 'r') {
        jump_betauzl = eq.Jump_betau_z(datal.coord.x_value, datal.coord.y_value,
                                       datal.coord.z_value);
        jump_betauzr = eq.Jump_betau_z(datar.coord.x_value, datar.coord.y_value,
                                       datar.coord.z_value);
      } else {
        jump_betauzl = inter.ifpz[ix][iy][ip].jump.u_dir;
        jump_betauzr = inter.ifpz[ix][iy][ip + 1].jump.u_dir;
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
          cor_row[i][j] = 0;
        }
      }

      // Approximate IZ
      iz = datal.left_loc;
      Operator_weights_z(beta, vdez, ix, iy, iz, dz);
      // IZ+1 cross left interface, use left interface's left FP
      cor_row[0][0] = -coef * dt * vdez[0];
      cor_row[0][1] = 1 - coef * dt * vdez[1];
      for (int i = 0; i < 5; i++) {
        cor_row[0][i] -= coef * dt * vdez[2] * datal.wei.weir[0][i];
      }
      cor_row[0][5] =
          uc[ix][iy][iz] +
          coef * dt * vdez[2] * (jump_ul * datal.wei.weir[0][5] +
                                 jump_betauzl * datal.wei.weir[0][6] +
                                 jump_ur * datal.wei.weir[0][7] +
                                 jump_betauzr * datal.wei.weir[0][8]);
      ;

      // Approximate IZ+1
      iz += 1;
      Operator_weights_z(beta, vdez, ix, iy, iz, dz);
      cor_row[1][2] = 1 - coef * dt * vdez[1];
      // IZ-1 cross left interface, use left interface's left FP
      for (int i = 0; i < 5; i++) {
        cor_row[1][i] -= coef * dt * vdez[0] * datal.wei.weil[0][i];
      }
      cor_row[1][5] =
          uc[ix][iy][iz] +
          coef * dt * vdez[0] * (jump_ul * datal.wei.weil[0][5] +
                                 jump_betauzl * datal.wei.weil[0][6] +
                                 jump_ur * datal.wei.weil[0][7] +
                                 jump_betauzr * datal.wei.weil[0][8]);
      // IZ+1 cross right interface, use right interface's right FP
      for (int i = 0; i < 5; i++) {
        cor_row[1][i] -= coef * dt * vdez[2] * datar.wei.weir[0][i];
      }
      cor_row[1][5] +=
          coef * dt * vdez[2] * (jump_ul * datar.wei.weir[0][5] +
                                 jump_betauzl * datar.wei.weir[0][6] +
                                 jump_ur * datar.wei.weir[0][7] +
                                 jump_betauzr * datar.wei.weir[0][8]);

      // Approximate IZ+2
      iz += 1;
      Operator_weights_z(beta, vdez, ix, iy, iz, dz);
      // IZ-1 cross right interface, use right interface's left FP
      cor_row[2][3] = 1 - coef * dt * vdez[1];
      cor_row[2][4] = -coef * dt * vdez[2];
      for (int i = 0; i < 5; i++) {
        cor_row[2][i] -= coef * dt * vdez[0] * datar.wei.weil[0][i];
      }
      cor_row[2][5] =
          uc[ix][iy][iz] +
          coef * dt * vdez[0] * (jump_ul * datar.wei.weil[0][5] +
                                 jump_betauzl * datar.wei.weil[0][6] +
                                 jump_ur * datar.wei.weil[0][7] +
                                 jump_betauzr * datar.wei.weil[0][8]);

      Convert2Tri_cor(datal.left_loc, cor_row, a, b, c, r);

      ip += 2;
    }
  }
}

/*********************************************************************************
 Source term initialization for Douglas ADI solver at each time step

 INPUT
 eq    : euqation object at next time step
 inter : object of all intersections

 OUTPUT
 Update Source terms in Douglas ADI at each time step
 ********************************************************************************/
void TS::Src_2nd(Equation &eq, Intersections &inter, CubicDoub &src) {
  int indx;

  // Initialize source term for all grid nodes
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        indx = To1d(ix, iy, iz);
        // Outside
        if (indicator[indx] == 1) {
          src[ix][iy][iz] = eq.Outer_f(xi[ix], yi[iy], zi[iz]);
        }
        // Inside
        else if (indicator[indx] == -1) {
          src[ix][iy][iz] = eq.Inner_f(xi[ix], yi[iy], zi[iz]);
        } else {
          cout << "Check mesh indx!" << endl;
          exit(0);
        }
      }
    }
  }
}

/**********************************************************************************************************
 Calculate weights with variable diffusion coefficients for 2nd order,
 beta*Delta_{xx}

 INPUT
 beta  : variable coefficients
 ix    : grid No. on x-axis
 iy    : grid No. on y-axis
 iz    : grid No. on z-axis
 dx    : grid mesh size

 OUTPUT
 vdex : weights for 2nd order 3-point stencil
 **********************************************************************************************************/
void TS::Operator_weights_x(Beta &beta, VecDoub_O &vdex, Int_I ix, Int_I iy,
                            Int_I iz, Doub_I dx) {
  int indx;
  double beta1, beta2;

  indx = To1d(ix, iy, iz);

  // Discrete outside grid node
  if (indicator[indx] > 0) {
    beta1 = beta.Outside(xi[ix] - dx / 2, yi[iy], zi[iz]);
    beta2 = beta.Outside(xi[ix] + dx / 2, yi[iy], zi[iz]);
  }
  // Discrete inside grid node
  else {
    beta1 = beta.Inside(xi[ix] - dx / 2, yi[iy], zi[iz]);
    beta2 = beta.Inside(xi[ix] + dx / 2, yi[iy], zi[iz]);
  }

  vdex[0] = beta1 / (dx * dx);
  vdex[1] = -(beta1 + beta2) / (dx * dx);
  vdex[2] = beta2 / (dx * dx);
}

/**********************************************************************************************************
 Calculate weights with variable diffusion coefficients for 2nd order,
 beta*Delta_{yy}

 INPUT
 beta  : variable coefficients
 ix    : grid No. on x-axis
 iy    : grid No. on y-axis
 iz    : grid No. on z-axis
 dy    : grid mesh size

 OUTPUT
 vdey : weights for 2nd order 3-point stencil
 **********************************************************************************************************/
void TS::Operator_weights_y(Beta &beta, VecDoub_O &vdey, Int_I ix, Int_I iy,
                            Int_I iz, Doub_I dy) {
  int indx;
  double beta1, beta2;

  indx = To1d(ix, iy, iz);

  // Discrete outside grid node
  if (indicator[indx] > 0) {
    beta1 = beta.Outside(xi[ix], yi[iy] - dy / 2, zi[iz]);
    beta2 = beta.Outside(xi[ix], yi[iy] + dy / 2, zi[iz]);
  }
  // Discrete inside grid node
  else {
    beta1 = beta.Inside(xi[ix], yi[iy] - dy / 2, zi[iz]);
    beta2 = beta.Inside(xi[ix], yi[iy] + dy / 2, zi[iz]);
  }

  vdey[0] = beta1 / (dy * dy);
  vdey[1] = -(beta1 + beta2) / (dy * dy);
  vdey[2] = beta2 / (dy * dy);
}

/**********************************************************************************************************
 Calculate weights with variable diffusion coefficients for 2nd order,
 beta*Delta_{zz}

 INPUT
 beta  : variable coefficients
 ix    : grid No. on x-axis
 iy    : grid No. on y-axis
 iz    : grid No. on z-axis
 dz    : grid mesh size

 OUTPUT
 vdez : weights for 2nd order 3-point stencil
 **********************************************************************************************************/
void TS::Operator_weights_z(Beta &beta, VecDoub_O &vdez, Int_I ix, Int_I iy,
                            Int_I iz, Doub_I dz) {
  int indx;
  double beta1, beta2;

  indx = To1d(ix, iy, iz);

  // Discrete outside grid node
  if (indicator[indx] > 0) {
    beta1 = beta.Outside(xi[ix], yi[iy], zi[iz] - dz / 2);
    beta2 = beta.Outside(xi[ix], yi[iy], zi[iz] + dz / 2);
  }
  // Discrete inside grid node
  else {
    beta1 = beta.Inside(xi[ix], yi[iy], zi[iz] - dz / 2);
    beta2 = beta.Inside(xi[ix], yi[iy], zi[iz] + dz / 2);
  }

  vdez[0] = beta1 / (dz * dz);
  vdez[1] = -(beta1 + beta2) / (dz * dz);
  vdez[2] = beta2 / (dz * dz);
}

/**********************************************************************************************************
 Convert matrix formed by irregular interface to tridiagonal matrix

 INPUT
 ix  : coordinate location of starting row
 row : matrix formed by corner interface

 OUTPUT
 a   : first hypotenuse of tridiagonal matrix
 b   : second hypotenuse of tridiagonal matrix
 c   : third hypotenuse of tridiagonal matrix
 r   : right hand side vector
 **********************************************************************************************************/
void TS::Convert2Tri_irr(Int_I ix, MatrixDoub_I &row, VecDoub_O &a,
                         VecDoub_O &b, VecDoub_O &c, VecDoub_O &r) {
  double rate;

  // Eliminate row[1][0]
  rate = row[1][0] / row[0][0];
  for (int i = 0; i < 5; i++) {
    row[1][i] -= rate * row[0][i];
  }

  // Eliminate row[0][3]
  rate = row[0][3] / row[1][3];
  for (int i = 0; i < 5; i++) {
    row[0][i] -= rate * row[1][i];
  }

  // Give back as tridiagonal matrix
  for (int i = 0; i < 2; i++) {
    a[ix + i] = row[i][i];
    b[ix + i] = row[i][i + 1];
    c[ix + i] = row[i][i + 2];
    r[ix + i] = row[i][4];
  }
}

/**********************************************************************************************************
 Convert matrix formed by corner interface to tridiagonal matrix

 INPUT
 ix  : coordinate location of starting row
 row : matrix formed by corner interface
 a   : first hypotenuse of tridiagonal matrix
 b   : second hypotenuse of tridiagonal matrix
 c   : third hypotenuse of tridiagonal matrix
 r   : right hand side vector
 **********************************************************************************************************/
void TS::Convert2Tri_cor(Int_I ix, MatrixDoub_I &row, VecDoub_O &a,
                         VecDoub_O &b, VecDoub_O &c, VecDoub_O &r) {
  double rate[3];

  // Eliminate row[1][0] and row[2][0]
  rate[1] = row[1][0] / row[0][0];
  rate[2] = row[2][0] / row[0][0];
  for (int i = 0; i < 6; i++) {
    row[1][i] -= rate[1] * row[0][i];
    row[2][i] -= rate[2] * row[0][i];
  }
  // Eliminate row[0][4] and row[1][4]
  rate[0] = row[0][4] / row[2][4];
  rate[1] = row[1][4] / row[2][4];
  for (int i = 0; i < 6; i++) {
    row[0][i] -= rate[0] * row[2][i];
    row[1][i] -= rate[1] * row[2][i];
  }

  // Eliminate row[2][1]
  rate[2] = row[2][1] / row[1][1];
  for (int i = 0; i < 6; i++) {
    row[2][i] -= rate[2] * row[1][i];
  }

  // Eliminate row[0][3]
  rate[0] = row[0][3] / row[1][3];
  for (int i = 0; i < 6; i++) {
    row[0][i] -= rate[0] * row[1][i];
  }

  // Give back as tridiagonal matrix
  for (int i = 0; i < 3; i++) {
    a[ix + i] = row[i][i];
    b[ix + i] = row[i][i + 1];
    c[ix + i] = row[i][i + 2];
    r[ix + i] = row[i][5];
  }
}

/**********************************************************************
 Tridiagonal matrix fast solver

 INPUT
 a : first hypotenuse of tridiagonal matrix, a[0] = 0
 b : second hypotenuse of tridiagonal matrix
 c : third hypotenuse of tridiagonal matrix, b[n-1] = 0
 r : right hand side vector

 OUTPUT
 u : solution of tridiagonal matrix system, Ax = b
 *********************************************************************/
void TS::TDMA(VecDoub_I &a, VecDoub_I &b, VecDoub_I &c, VecDoub_I &r,
              VecDoub_O &u) {
  int n = (int)a.size();
  double bet;
  VecDoub gam(n); // One vector of workspace, gam, is needed.

  if (b[0] == 0.0) {
    throw("Error 1 in tridag");
  }
  // If this happens, then you should rewrite your equations as a set of order
  // N-1, with u1 trivially eliminated.
  u[0] = r[0] / (bet = b[0]);
  // Decomposition and forward substitution.
  for (int j = 1; j < n; j++) {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0) {
      throw("Error 2 in tridag"); // Algorithm fails
    }
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }
  for (int j = (n - 2); j >= 0; j--)
    u[j] -= gam[j + 1] * u[j + 1]; // Backsubstitution.
}

/******************************************************************************
 Calculate weights of Lagrange polynomial interpolation

 INPUT
 z: location where approximatetions are to be accurated
 x: x(n) grid points relative locations
 n: dimension of vector x(n), numbers of grid points
 m: highest derivative for which weghts are sought

 OUTPUT
 c: c(0:n-1,0:m) weights at grid locations x(n) for derivatives of order 0:m
 ******************************************************************************/
void TS::Weights(Doub_I z, VecDoub_I &x, Int_I n, Int_I m, MatrixDoub_O &c) {
  double c1, c2, c3, c4, c5;
  int mn;

  c1 = 1;
  c4 = x[0] - z;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m + 1; j++) {
      c[i][j] = 0;
    }
  }

  c[0][0] = 1;

  for (int i = 1; i < n; i++) {
    mn = min(i, m);
    c2 = 1;
    c5 = c4;
    c4 = x[i] - z;
    for (int j = 0; j < i; j++) {
      c3 = x[i] - x[j];
      c2 = c2 * c3;
      if (j == i - 1) {
        for (int k = mn; k > 0; k--) {
          c[i][k] = c1 * (k * c[i - 1][k - 1] - c5 * c[i - 1][k]) / c2;
        }
        c[i][0] = -c1 * c5 * c[i - 1][0] / c2;
      }
      for (int k = mn; k > 0; k--) {
        c[j][k] = (c4 * c[j][k] - k * c[j][k - 1]) / c3;
      }
      c[j][0] = c4 * c[j][0] / c3;
    }
    c1 = c2;
  }
}

/**********************************************************
 Convert 3D coefficients to 1D coefficients

 INPUT
 ix : index of grid node on x-axis
 iy : index of grid node on y-axis
 iz : index of grid node on z-axis

 OUTPUT
 1D coefficients converted from 3D
 **********************************************************/
int TS::To1d(Int_I ix, Int_I iy, Int_I iz) {
  int temp;

  temp = (iz * nx * ny) + (iy * nx) + ix;

  return temp;
}

/************************************************************************
 Initialization of analytical solution

 INPUT
 eq : equation object
 uh : three-dimensional analytical solution
 ************************************************************************/
void TS::Initialization(Equation &eq, CubicDoub &uh) {
  int indx;

  uh.resize(nx);
  for (int ix = 0; ix < nx; ix++) {
    uh[ix].resize(ny);
    for (int iy = 0; iy < ny; iy++) {
      uh[ix][iy].resize(nz);
    }
  }

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        indx = To1d(ix, iy, iz);

        // Outside
        if (indicator[indx] == 1) {
          uh[ix][iy][iz] = eq.Outer_u(xi[ix], yi[iy], zi[iz]);
        }
        // Inside
        else if (indicator[indx] == -1) {
          uh[ix][iy][iz] = eq.Inner_u(xi[ix], yi[iy], zi[iz]);
        } else {
          cout << "Check mesh indx!" << endl;
          exit(0);
        }
      }
    }
  }
}

/************************************************************************
 Error of numerical solution, show errors in L^{2} and L^{infinity} norm

 INPUT
 eq : equation object
 uh : three-dimensional analytical solution
 ************************************************************************/
void TS::Error(Equation &eq, CubicDoub &uh, ofstream &out_file) {
  int indx;
  double l2, lmax, temp;

  int xmax, ymax, zmax;

  lmax = 0.0;
  l2 = 0.0;

  xmax = 0;
  ymax = 0;
  zmax = 0;

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        indx = To1d(ix, iy, iz);

        // Outside
        if (indicator[indx] == 1) {
          temp = abs(uh[ix][iy][iz] - eq.Outer_u(xi[ix], yi[iy], zi[iz]));
        }
        // Inside
        else if (indicator[indx] == -1) {
          temp = abs(uh[ix][iy][iz] - eq.Inner_u(xi[ix], yi[iy], zi[iz]));
        } else {
          cout << "Check mesh indx!" << endl;
          exit(0);
        }

        if (temp > lmax) {
          lmax = temp;
          xmax = ix;
          ymax = iy;
          zmax = iz;
        }
        l2 += temp * temp;
      }
    }
  }
  l2 = sqrt(l2 / (1.0 * nx * ny * nz));

  out_file << setprecision(3) << scientific;
  out_file << "[Xmax,Ymax,Zmax] = [" << xmax << "," << ymax << "," << zmax
           << "]" << endl;
  out_file << " Lmax = " << lmax << ";  L2 = " << l2 << endl;
  out_file << fixed << endl;
}
