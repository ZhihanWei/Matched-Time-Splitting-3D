#include "Douglas_ADI.h"
#include <iomanip>
#include <iostream>

/******************************************************************************************
                                 Constructor

 INPUT
 mesh     : mesh object
 inter    : object of all intersections
 in_beta  : vector of 2 double values represent beta^{-} and beta^{+}
 time     : vector of 3 double values represent beginning time, finishing time and time step
 accuracy : accuracy of scheme
 *******************************************************************************************/
Douglas_ADI::Douglas_ADI(Intersections& inter, Mesh& mesh, Beta& beta,
                         VecDoub_I time, Int_I accuracy) {
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

  // Initialize size for solution
  uhs.resize(nx);
  uhs2.resize(nx);
  uhss.resize(nx);
  uh1.resize(nx);
  src.resize(nx);
    
  for (int ix = 0; ix < nx; ix++) {
    uhs[ix].resize(ny);
    uhs2[ix].resize(ny);
    uhss[ix].resize(ny);
    uh1[ix].resize(ny);
    src[ix].resize(ny);
    for (int iy = 0; iy < ny; iy++) {
      uhs[ix][iy].resize(nz);
      uhs2[ix][iy].resize(nz);
      uhss[ix][iy].resize(nz);
      uh1[ix][iy].resize(nz);
      src[ix][iy].resize(nz);
    }
  }
}

/*********************************************************************************
 Douglas ADI solver at each time step

 INPUT
 eq    : euqation object at next time step
 inter : object of all intersections

 OUTPUT
 uh : three-dimensional solution at current time step to next time step
 ********************************************************************************/
void Douglas_ADI::Solve_2nd(Equation& eq, Intersections& inter, CubicDoub& uh,
                            Beta& beta) {
  int ip;
  int ix, iy, iz;
  double sum;
  VecDoub vdex, vdey, vdez;
  VecDoub a, b, c, r, ut;
  MatrixDoub irr_row, cor_row;
  Intersection_Data data, datal, datar;

  inter.Refresh_Jump(eq, uh, beta);
  Src_2nd(eq, inter, beta, uh);

  // Initialize Central finite difference weights for 3-points stencil
  vdex.resize(3);
  vdey.resize(3);
  vdez.resize(3);

  /*******************************************************************************************************
   First Step:
   (1 - Dt*beta*D_{xx})U^{*} = Dt*beta*D_{yy}*U^{N} + Dt*beta*D_{zz}*U^{N} + U^{N} + Dt*SRC^{N+1}
   ******************************************************************************************************/
  // Generate RHS

  // Step I: Dt*beta*D_{yy}*U^{N}
  // Step I-1: Dt*beta*D_{yy}*U^{N} without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          sum += uh[ix][iy + i][iz] * vdey[i + 1];
        }

        uhs[ix][iy][iz] = dt * sum;
      }
    }
  }

  // Step I-2: Apply MIB to the first part of RHS
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      ip = 0;
      while (ip < inter.ifpy[ix][iz].size()) {
        // An irregular interface point
        if (inter.ifpy[ix][iz][ip].ID > 0) {
          data = inter.ifpy[ix][iz][ip];

          // Approximate IY
          iy = data.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdey[2] * uh[ix][data.left_loc - 1 + i][iz] *
                   data.wei.weir[0][i];
          }
          uhs[ix][iy][iz] = sum * dt;

          // Approximate IY+1
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdey[0] * uh[ix][data.left_loc - 1 + i][iz] *
                   data.wei.weil[0][i];
          }
          uhs[ix][iy][iz] = sum * dt;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpy[ix][iz][ip];
          datar = inter.ifpy[ix][iz][ip + 1];

          // Approximate IY
          iy = datal.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdey[2] * uh[ix][datal.left_loc - 1 + i][iz] *
                   datal.wei.weir[0][i];
          }
          uhs[ix][iy][iz] = sum * dt;

          // Approximate IY+1, the inside corner point
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          sum = 0;
          // IY-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdey[0] * uh[ix][datal.left_loc - 1 + i][iz] *
                   datal.wei.weil[0][i];
          }
          // IY,  no FP
          sum += vdey[1] * uh[ix][iy][iz];
          // IY+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdey[2] * uh[ix][datar.left_loc - 2 + i][iz] *
                   datar.wei.weir[0][i];
          }
          uhs[ix][iy][iz] = sum * dt;

          // Approximate IY+2, right FP
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdey[0] * uh[ix][datar.left_loc - 2 + i][iz] *
                   datar.wei.weil[0][i];
          }
          uhs[ix][iy][iz] = sum * dt;

          ip += 2;
        }
      }
    }
  }

  // Step II: Dt*beta*D_{zz}*U^{N}
  // Step II-1: Dt*beta*D_{zz}*U^{N} without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          sum += uh[ix][iy][iz + i] * vdez[i + 1];
        }
        uhs2[ix][iy][iz] = dt * sum;
      }
    }
  }

  // Step II-2: Apply MIB to the second part of RHS
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      ip = 0;
      while (ip < inter.ifpz[ix][iy].size()) {
        // An irregular interface point
        if (inter.ifpz[ix][iy][ip].ID > 0) {
          data = inter.ifpz[ix][iy][ip];

          // Approximate IZ
          iz = data.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdez[2] * uh[ix][iy][data.left_loc - 1 + i] *
                   data.wei.weir[0][i];
          }
          uhs2[ix][iy][iz] = sum * dt;

          // Approximate IZ+1
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdez[0] * uh[ix][iy][data.left_loc - 1 + i] *
                   data.wei.weil[0][i];
          }
          uhs2[ix][iy][iz] = sum * dt;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpz[ix][iy][ip];
          datar = inter.ifpz[ix][iy][ip + 1];

          // Approximate IZ, left FP
          iz = datal.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdez[2] * uh[ix][iy][datal.left_loc - 1 + i] *
                   datal.wei.weir[0][i];
          }
          uhs2[ix][iy][iz] = sum * dt;

          // Approximate IZ+1, the inside corner point
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          sum = 0;
          // IZ-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdez[0] * uh[ix][iy][datal.left_loc - 1 + i] *
                   datal.wei.weil[0][i];
          }
          // IZ,  no FP
          sum += vdez[1] * uh[ix][iy][iz];
          // IZ+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdez[2] * uh[ix][iy][datar.left_loc - 2 + i] *
                   datar.wei.weir[0][i];
          }
          uhs2[ix][iy][iz] = sum * dt;

          // Approximate IZ+2, right FP
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdez[0] * uh[ix][iy][datar.left_loc - 2 + i] *
                   datar.wei.weil[0][i];
          }
          uhs2[ix][iy][iz] = sum * dt;

          ip += 2;
        }
      }
    }
  }

  // Add Step I and Step II
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uhs[ix][iy][iz] += uhs2[ix][iy][iz];
      }
    }
  }

  // Step III: U^{N}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uhs[ix][iy][iz] += uh[ix][iy][iz];
      }
    }
  }

  // Step IV: Dt*SRC^{N+1}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uhs[ix][iy][iz] += dt * src[ix][iy][iz];
      }
    }
  }

  // Step V: set boundary conditions for uhs
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      uhs[ix][iy][0] = eq.Outer_u(xi[ix], yi[iy], zi[0]);
      uhs[ix][iy][nz - 1] = eq.Outer_u(xi[ix], yi[iy], zi[nz - 1]);
    }
  }
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      uhs[ix][0][iz] = eq.Outer_u(xi[ix], yi[0], zi[iz]);
      uhs[ix][ny - 1][iz] = eq.Outer_u(xi[ix], yi[ny - 1], zi[iz]);
    }
  }
  for (int iy = 0; iy < ny; iy++) {
    for (int iz = 0; iz < nz; iz++) {
      uhs[0][iy][iz] = eq.Outer_u(xi[0], yi[iy], zi[iz]);
      uhs[nx - 1][iy][iz] = eq.Outer_u(xi[nx - 1], yi[iy], zi[iz]);
    }
  }

  // Generate LHS
  a.resize(nx);
  b.resize(nx);
  c.resize(nx);
  r.resize(nx);
  ut.resize(nx);
  // set up irr_row and cor_row
  irr_row.resize(2);
  for (int i = 0; i < 2; i++) {
    irr_row[i].resize(5);
  }
  cor_row.resize(3);
  for (int i = 0; i < 3; i++) {
    cor_row[i].resize(6);
  }

  for (int iy = 1; iy < ny - 1; iy++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      // Step I: LHS without MIB
      a[0] = 0.0;
      b[0] = 1.0;
      c[0] = 0.0;
      for (int ix = 1; ix < nx - 1; ix++) {
        Operator_weights_x(beta, vdex, ix, iy, iz, dx);

        a[ix] = -dt * vdex[0];
        b[ix] = 1 - dt * vdex[1];
        c[ix] = -dt * vdex[2];
      }
      a[nx - 1] = 0.0;
      b[nx - 1] = 1.0;
      c[nx - 1] = 0.0;

      for (int ix = 0; ix < nx; ix++) {
        r[ix] = uhs[ix][iy][iz];
      }

      // Step II: LHS with MIB
      ip = 0;
      while (ip < inter.ifpx[iy][iz].size()) {
        // An irregular interface
        if (inter.ifpx[iy][iz][ip].ID > 0) {
          data = inter.ifpx[iy][iz][ip];

          for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 5; j++) {
              irr_row[i][j] = 0;
            }
          }

          // Approximate IX
          ix = data.left_loc;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // use right FP
          irr_row[0][0] = -dt * vdex[0];
          irr_row[0][1] = 1 - dt * vdex[1];
          for (int i = 0; i < 4; i++) {
            irr_row[0][i] += -dt * vdex[2] * data.wei.weir[0][i];
          }
          irr_row[0][4] = uhs[ix][iy][iz];

          // Approximate IX+1
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // use left FP
          irr_row[1][2] = 1 - dt * vdex[1];
          irr_row[1][3] = -dt * vdex[2];
          for (int i = 0; i < 4; i++) {
            irr_row[1][i] += -dt * vdex[0] * data.wei.weil[0][i];
          }
          irr_row[1][4] = uhs[ix][iy][iz];

          Convert2Tri_irr(data.left_loc, irr_row, a, b, c, r);

          ip += 1;
        }
        // Corner interface
        else {
          datal = inter.ifpx[iy][iz][ip];
          datar = inter.ifpx[iy][iz][ip + 1];

          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 6; j++) {
              cor_row[i][j] = 0;
            }
          }

          // Approximate IX
          ix = datal.left_loc;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX+1 cross left interface, use left interface's left FP
          cor_row[0][0] = -dt * vdex[0];
          cor_row[0][1] = 1 - dt * vdex[1];
          for (int i = 0; i < 5; i++) {
            cor_row[0][i] -= dt * vdex[2] * datal.wei.weir[0][i];
          }
          cor_row[0][5] = uhs[ix][iy][iz];

          // Approximate IX+1
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          cor_row[1][2] = 1 - dt * vdex[1];
          // IX-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            cor_row[1][i] -= dt * vdex[0] * datal.wei.weil[0][i];
          }
          // IX+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            cor_row[1][i] -= dt * vdex[2] * datar.wei.weir[0][i];
          }
          cor_row[1][5] = uhs[ix][iy][iz];

          // Approximate IX+2
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX-1 cross right interface, use right interface's left FP
          cor_row[2][3] = 1 - dt * vdex[1];
          cor_row[2][4] = -dt * vdex[2];
          for (int i = 0; i < 5; i++) {
            cor_row[2][i] -= dt * vdex[0] * datar.wei.weil[0][i];
          }
          cor_row[2][5] = uhs[ix][iy][iz];

          Convert2Tri_cor(datal.left_loc, cor_row, a, b, c, r);

          ip += 2;
        }
      }

      // Thomas Algorithm
      TDMA(a, b, c, r, ut);

      for (int ix = 0; ix < nx; ix++) {
        uhs[ix][iy][iz] = ut[ix];
      }

      /*
      ofstream out_file;
      string out_file_name;

      out_file_name = "result/<first_step>.txt";
      out_file.open(out_file_name, ios::out | ios::app);
      if(out_file.good())
      {
          for(int ix = 0; ix < nx; ix++)
          {
              out_file << ut[ix] << endl;

          }
      }
      out_file.close();
       */
    }
  }

  /*******************************************************************************************************
   Second Step:
   (1 - Dt*beta*D_{yy})U^{**} = U^{*} - Dt*beta*D_{yy}*U^{N}
   ******************************************************************************************************/
  // Generate RHS

  // Step I: -Dt*beta*D_{yy}*U^{N}
  // Step I-1: -Dt*beta*D_{yy}*U^{N} without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);

          sum += uh[ix][iy + i][iz] * vdey[i + 1];
        }
        uhss[ix][iy][iz] = -dt * sum;
      }
    }
  }

  // Step I-2: Apply MIB to the first part of RHS
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      ip = 0;
      while (ip < inter.ifpy[ix][iz].size()) {
        // An irregular interface point
        if (inter.ifpy[ix][iz][ip].ID > 0) {
          data = inter.ifpy[ix][iz][ip];

          // Approximate IY
          iy = data.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdey[2] * uh[ix][data.left_loc - 1 + i][iz] *
                   data.wei.weir[0][i];
          }
          uhss[ix][iy][iz] = -dt * sum;

          // Approximate IY+1
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdey[0] * uh[ix][data.left_loc - 1 + i][iz] *
                   data.wei.weil[0][i];
          }
          uhss[ix][iy][iz] = -dt * sum;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpy[ix][iz][ip];
          datar = inter.ifpy[ix][iz][ip + 1];

          // Approximate IY, left FP
          iy = datal.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdey[2] * uh[ix][datal.left_loc - 1 + i][iz] *
                   datal.wei.weir[0][i];
          }
          uhss[ix][iy][iz] = -dt * sum;

          // Approximate IY+1, the inside corner point
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          sum = 0;
          // IY-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdey[0] * uh[ix][datal.left_loc - 1 + i][iz] *
                   datal.wei.weil[0][i];
          }
          // IY,  no FP
          sum += vdey[1] * uh[ix][iy][iz];
          // IY+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdey[2] * uh[ix][datar.left_loc - 2 + i][iz] *
                   datar.wei.weir[0][i];
          }
          uhss[ix][iy][iz] = -dt * sum;

          // Approximate IY+2, right FP
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdey[i + 1] * uh[ix][iy + i][iz];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdey[0] * uh[ix][datar.left_loc - 2 + i][iz] *
                   datar.wei.weil[0][i];
          }
          uhss[ix][iy][iz] = -dt * sum;

          ip += 2;
        }
      }
    }
  }

  // Step II: U^{*}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uhss[ix][iy][iz] += uhs[ix][iy][iz];
      }
    }
  }

  // Step V: set boundary conditions for uhss
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      uhss[ix][iy][0] = eq.Outer_u(xi[ix], yi[iy], zi[0]);
      uhss[ix][iy][nz - 1] = eq.Outer_u(xi[ix], yi[iy], zi[nz - 1]);
    }
  }
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      uhss[ix][0][iz] = eq.Outer_u(xi[ix], yi[0], zi[iz]);
      uhss[ix][ny - 1][iz] = eq.Outer_u(xi[ix], yi[ny - 1], zi[iz]);
    }
  }
  for (int iy = 0; iy < ny; iy++) {
    for (int iz = 0; iz < nz; iz++) {
      uhss[0][iy][iz] = eq.Outer_u(xi[0], yi[iy], zi[iz]);
      uhss[nx - 1][iy][iz] = eq.Outer_u(xi[nx - 1], yi[iy], zi[iz]);
    }
  }

  // Generate LHS
  a.resize(ny);
  b.resize(ny);
  c.resize(ny);
  r.resize(ny);
  ut.resize(ny);

  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      // Step I: LHS without MIB
      a[0] = 0.0;
      b[0] = 1.0;
      c[0] = 0.0;
      for (int iy = 1; iy < ny - 1; iy++) {
        Operator_weights_y(beta, vdey, ix, iy, iz, dy);

        a[iy] = -dt * vdey[0];
        b[iy] = 1 - dt * vdey[1];
        c[iy] = -dt * vdey[2];
      }
      a[ny - 1] = 0.0;
      b[ny - 1] = 1.0;
      c[ny - 1] = 0.0;

      for (int iy = 0; iy < ny; iy++) {
        r[iy] = uhss[ix][iy][iz];
      }

      // Step II: LHS with MIB
      ip = 0;
      while (ip < inter.ifpy[ix][iz].size()) {
        // An irregular interface
        if (inter.ifpy[ix][iz][ip].ID > 0) {
          data = inter.ifpy[ix][iz][ip];

          for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 5; j++) {
              irr_row[i][j] = 0;
            }
          }

          // Approximate IY
          iy = data.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // use right FP
          irr_row[0][0] = -dt * vdey[0];
          irr_row[0][1] = 1 - dt * vdey[1];
          for (int i = 0; i < 4; i++) {
            irr_row[0][i] += -dt * vdey[2] * data.wei.weir[0][i];
          }
          irr_row[0][4] = uhss[ix][iy][iz];

          // Approximate IY+1
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // use left FP
          irr_row[1][2] = 1 - dt * vdey[1];
          irr_row[1][3] = -dt * vdey[2];
          for (int i = 0; i < 4; i++) {
            irr_row[1][i] += -dt * vdey[0] * data.wei.weil[0][i];
          }
          irr_row[1][4] = uhss[ix][iy][iz];

          Convert2Tri_irr(data.left_loc, irr_row, a, b, c, r);

          ip += 1;
        }
        // Corner interface
        else {
          datal = inter.ifpy[ix][iz][ip];
          datar = inter.ifpy[ix][iz][ip + 1];

          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 6; j++) {
              cor_row[i][j] = 0;
            }
          }

          // Approximate IY
          iy = datal.left_loc;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY+1 cross left interface, use left interface's right FP
          cor_row[0][0] = -dt * vdey[0];
          cor_row[0][1] = 1 - dt * vdey[1];
          for (int i = 0; i < 5; i++) {
            cor_row[0][i] -= dt * vdey[2] * datal.wei.weir[0][i];
          }
          cor_row[0][5] = uhss[ix][iy][iz];

          // Approximate IY+1
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          cor_row[1][2] = 1 - dt * vdey[1];
          // IY-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            cor_row[1][i] -= dt * vdey[0] * datal.wei.weil[0][i];
          }
          // IY+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            cor_row[1][i] -= dt * vdey[2] * datar.wei.weir[0][i];
          }
          cor_row[1][5] = uhss[ix][iy][iz];

          // Approximate IY+2
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross right interface, use right interface's left FP
          cor_row[2][3] = 1 - dt * vdey[1];
          cor_row[2][4] = -dt * vdey[2];
          for (int i = 0; i < 5; i++) {
            cor_row[2][i] -= dt * vdey[0] * datar.wei.weil[0][i];
          }
          cor_row[2][5] = uhss[ix][iy][iz];

          Convert2Tri_cor(datal.left_loc, cor_row, a, b, c, r);

          ip += 2;
        }
      }

      // Thomas Algorithm
      TDMA(a, b, c, r, ut);

      for (int iy = 0; iy < ny; iy++) {
        uhss[ix][iy][iz] = ut[iy];
      }

      /*
      ofstream out_file;
      string out_file_name;

      out_file_name = "result/<first_step>.txt";
      out_file.open(out_file_name, ios::out | ios::app);
      if(out_file.good())
      {
          for(int iy = 0; iy < ny; iy++)
          {
              out_file << ut[iy] << endl;

          }
      }
      out_file.close();
       */
    }
  }

  /*******************************************************************************************************
   Third Step:
   (1 - Dt*beta*D_{zz})U^{N+1} = U^{**} - Dt*beta*D_{zz}*U^{N}
   ******************************************************************************************************/
  // Generate RHS

  // Step II: Dt*beta*D_{zz}*U^{N}
  // Step II-1: Dt*beta*D_{zz}*U^{N} without MIB
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        sum = 0;
        for (int i = -1; i < 2; i++) {
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);

          sum += uh[ix][iy][iz + i] * vdez[i + 1];
        }
        uh1[ix][iy][iz] = -dt * sum;
      }
    }
  }

  // Step II-2: Apply MIB to the second part of RHS
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      ip = 0;
      while (ip < inter.ifpz[ix][iy].size()) {
        // An irregular interface point
        if (inter.ifpz[ix][iy][ip].ID > 0) {
          data = inter.ifpz[ix][iy][ip];

          // Approximate IZ
          iz = data.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross interface, use right FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdez[2] * uh[ix][iy][data.left_loc - 1 + i] *
                   data.wei.weir[0][i];
          }
          uh1[ix][iy][iz] = -dt * sum;

          // Approximate IZ+1
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross interface, use left FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 4; i++) {
            sum += vdez[0] * uh[ix][iy][data.left_loc - 1 + i] *
                   data.wei.weil[0][i];
          }
          uh1[ix][iy][iz] = -dt * sum;

          ip += 1;
        }
        // A corner interface point
        else {
          datal = inter.ifpz[ix][iy][ip];
          datar = inter.ifpz[ix][iy][ip + 1];

          // Approximate IZ, left FP
          iz = datal.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross left interface, use middle FP
          sum = 0;
          for (int i = -1; i < 1; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdez[2] * uh[ix][iy][datal.left_loc - 1 + i] *
                   datal.wei.weir[0][i];
          }
          uh1[ix][iy][iz] = -dt * sum;

          // Approximate IZ+1, the inside corner point
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          sum = 0;
          // IZ-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            sum += vdez[0] * uh[ix][iy][datal.left_loc - 1 + i] *
                   datal.wei.weil[0][i];
          }
          // IZ,  no FP
          sum += vdez[1] * uh[ix][iy][iz];
          // IZ+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            sum += vdez[2] * uh[ix][iy][datar.left_loc - 2 + i] *
                   datar.wei.weir[0][i];
          }
          uh1[ix][iy][iz] = -dt * sum;

          // Approximate IZ+2, right FP
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross right interface, use middle FP
          sum = 0;
          for (int i = 0; i < 2; i++) {
            sum += vdez[i + 1] * uh[ix][iy][iz + i];
          }
          for (int i = 0; i < 5; i++) {
            sum += vdez[0] * uh[ix][iy][datar.left_loc - 2 + i] *
                   datar.wei.weil[0][i];
          }
          uh1[ix][iy][iz] = -dt * sum;

          ip += 2;
        }
      }
    }
  }

  // Step III: U^{**}
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      for (int iz = 1; iz < nz - 1; iz++) {
        uh1[ix][iy][iz] += uhss[ix][iy][iz];
      }
    }
  }

  // Step V: set boundary conditions for uh1
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      uh1[ix][iy][0] = eq.Outer_u(xi[ix], yi[iy], zi[0]);
      uh1[ix][iy][nz - 1] = eq.Outer_u(xi[ix], yi[iy], zi[nz - 1]);
    }
  }
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      uh1[ix][0][iz] = eq.Outer_u(xi[ix], yi[0], zi[iz]);
      uh1[ix][ny - 1][iz] = eq.Outer_u(xi[ix], yi[ny - 1], zi[iz]);
    }
  }
  for (int iy = 0; iy < ny; iy++) {
    for (int iz = 0; iz < nz; iz++) {
      uh1[0][iy][iz] = eq.Outer_u(xi[0], yi[iy], zi[iz]);
      uh1[nx - 1][iy][iz] = eq.Outer_u(xi[nx - 1], yi[iy], zi[iz]);
    }
  }

  // Generate LHS
  a.resize(nz);
  b.resize(nz);
  c.resize(nz);
  r.resize(nz);
  ut.resize(nz);

  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iy = 1; iy < ny - 1; iy++) {
      // Step I: LHS without MIB
      a[0] = 0.0;
      b[0] = 1.0;
      c[0] = 0.0;
      for (int iz = 1; iz < nz - 1; iz++) {
        Operator_weights_z(beta, vdez, ix, iy, iz, dz);

        a[iz] = -dt * vdez[0];
        b[iz] = 1 - dt * vdez[1];
        c[iz] = -dt * vdez[2];
      }
      a[nz - 1] = 0.0;
      b[nz - 1] = 1.0;
      c[nz - 1] = 0.0;

      for (int iz = 0; iz < nz; iz++) {
        r[iz] = uh1[ix][iy][iz];
      }

      // Step II: LHS with MIB
      ip = 0;
      while (ip < inter.ifpz[ix][iy].size()) {
        // An irregular interface
        if (inter.ifpz[ix][iy][ip].ID > 0) {
          data = inter.ifpz[ix][iy][ip];

          for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 5; j++) {
              irr_row[i][j] = 0;
            }
          }

          // Approximate IZ
          iz = data.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // use right FP
          irr_row[0][0] = -dt * vdez[0];
          irr_row[0][1] = 1 - dt * vdez[1];
          for (int i = 0; i < 4; i++) {
            irr_row[0][i] += -dt * vdez[2] * data.wei.weir[0][i];
          }
          irr_row[0][4] = uh1[ix][iy][iz];

          // Approximate IZ+1
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // use left FP
          irr_row[1][2] = 1 - dt * vdez[1];
          irr_row[1][3] = -dt * vdez[2];
          for (int i = 0; i < 4; i++) {
            irr_row[1][i] += -dt * vdez[0] * data.wei.weil[0][i];
          }
          irr_row[1][4] = uh1[ix][iy][iz];

          Convert2Tri_irr(data.left_loc, irr_row, a, b, c, r);

          ip += 1;
        }
        // Corner interface
        else {
          datal = inter.ifpz[ix][iy][ip];
          datar = inter.ifpz[ix][iy][ip + 1];

          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 6; j++) {
              cor_row[i][j] = 0;
            }
          }

          // Approximate IZ
          iz = datal.left_loc;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ+1 cross left interface, use left interface's left FP
          cor_row[0][0] = -dt * vdez[0];
          cor_row[0][1] = 1 - dt * vdez[1];
          for (int i = 0; i < 5; i++) {
            cor_row[0][i] -= dt * vdez[2] * datal.wei.weir[0][i];
          }
          cor_row[0][5] = uh1[ix][iy][iz];

          // Approximate IZ+1
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          cor_row[1][2] = 1 - dt * vdez[1];
          // IZ-1 cross left interface, use left interface's left FP
          for (int i = 0; i < 5; i++) {
            cor_row[1][i] -= dt * vdez[0] * datal.wei.weil[0][i];
          }
          // IZ+1 cross right interface, use right interface's right FP
          for (int i = 0; i < 5; i++) {
            cor_row[1][i] -= dt * vdez[2] * datar.wei.weir[0][i];
          }
          cor_row[1][5] = uh1[ix][iy][iz];

          // Approximate IZ+2
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross right interface, use right interface's left FP
          cor_row[2][3] = 1 - dt * vdez[1];
          cor_row[2][4] = -dt * vdez[2];
          for (int i = 0; i < 5; i++) {
            cor_row[2][i] -= dt * vdez[0] * datar.wei.weil[0][i];
          }
          cor_row[2][5] = uh1[ix][iy][iz];

          Convert2Tri_cor(datal.left_loc, cor_row, a, b, c, r);

          ip += 2;
        }
      }

      // Thomas Algorithm
      TDMA(a, b, c, r, ut);

      for (int iz = 0; iz < nz; iz++) {
        uh1[ix][iy][iz] = ut[iz];
      }

      /*
      ofstream out_file;
      string out_file_name;

      out_file_name = "result/<first_step>.txt";
      out_file.open(out_file_name, ios::out | ios::app);
      if(out_file.good())
      {
          for(int iz = 0; iz < nz; iz++)
          {
              out_file << ut[iz] << endl;

          }
      }
      out_file.close();
       */
    }
  }

  // Update UH to next time step
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        uh[ix][iy][iz] = uh1[ix][iy][iz];
      }
    }
  }
}

/*********************************************************************************
 Source term initialization for Douglas ADI solver at each time step

 INPUT
 eq    : euqation object at next time step
 inter : object of all intersections
 uh    : three-dimensional solution at current time step to next time step

 OUTPUT
 Update Source terms in Douglas ADI at each time step
 ********************************************************************************/
void Douglas_ADI::Src_2nd(Equation& eq, Intersections& inter, Beta& beta,
                          CubicDoub_I& uh) {
  int indx, ip;
  int ix, iy, iz;
  VecDoub vdex, vdey, vdez;
  Intersection_Data data, datal, datar;

  // Initialize Central finite difference weights for 3-points stencil
  vdex.resize(3);
  vdey.resize(3);
  vdez.resize(3);

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

  // Step I: Apply MIB to source term in X-direction
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
          src[ix][iy][iz] += vdex[2] * (jump_u * data.wei.weir[0][4] +
                                        jump_betaux * data.wei.weir[0][5]);

          // Approximate IX+1
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX-1 cross interface, use left FP
          src[ix][iy][iz] += vdex[0] * (jump_u * data.wei.weil[0][4] +
                                        jump_betaux * data.wei.weil[0][5]);

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
          src[ix][iy][iz] += vdex[2] * (jump_ul * datal.wei.weir[0][5] +
                                        jump_betauxl * datal.wei.weir[0][6] +
                                        jump_ur * datal.wei.weir[0][7] +
                                        jump_betauxr * datal.wei.weir[0][8]);

          // Approximate IX+1, the inside corner point
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX-1 cross left interface, use left interface's left FP
          src[ix][iy][iz] += vdex[0] * (jump_ul * datal.wei.weil[0][5] +
                                        jump_betauxl * datal.wei.weil[0][6] +
                                        jump_ur * datal.wei.weil[0][7] +
                                        jump_betauxr * datal.wei.weil[0][8]);
          // IX+1 cross right interface, use right interface's right FP
          src[ix][iy][iz] += vdex[2] * (jump_ul * datar.wei.weir[0][5] +
                                        jump_betauxl * datar.wei.weir[0][6] +
                                        jump_ur * datar.wei.weir[0][7] +
                                        jump_betauxr * datar.wei.weir[0][8]);

          // Approximate IX+2, right FP
          ix += 1;
          Operator_weights_x(beta, vdex, ix, iy, iz, dx);
          // IX-1 cross right interface, use middle FP
          src[ix][iy][iz] += vdex[0] * (jump_ul * datar.wei.weil[0][5] +
                                        jump_betauxl * datar.wei.weil[0][6] +
                                        jump_ur * datar.wei.weil[0][7] +
                                        jump_betauxr * datar.wei.weil[0][8]);

          ip += 2;
        }
      }
    }
  }

  // Step II: Apply MIB to source term in Y-direction
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
          src[ix][iy][iz] += vdey[2] * (jump_u * data.wei.weir[0][4] +
                                        jump_betauy * data.wei.weir[0][5]);

          // Approximate IY+1
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross interface, use left FP
          src[ix][iy][iz] += vdey[0] * (jump_u * data.wei.weil[0][4] +
                                        jump_betauy * data.wei.weil[0][5]);

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
          src[ix][iy][iz] += vdey[2] * (jump_ul * datal.wei.weir[0][5] +
                                        jump_betauyl * datal.wei.weir[0][6] +
                                        jump_ur * datal.wei.weir[0][7] +
                                        jump_betauyr * datal.wei.weir[0][8]);

          // Approximate IY+1, the inside corner point
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross left interface, use left interface's left FP
          src[ix][iy][iz] += vdey[0] * (jump_ul * datal.wei.weil[0][5] +
                                        jump_betauyl * datal.wei.weil[0][6] +
                                        jump_ur * datal.wei.weil[0][7] +
                                        jump_betauyr * datal.wei.weil[0][8]);
          // IY+1 cross right interface, use right interface's right FP
          src[ix][iy][iz] += vdey[2] * (jump_ul * datar.wei.weir[0][5] +
                                        jump_betauyl * datar.wei.weir[0][6] +
                                        jump_ur * datar.wei.weir[0][7] +
                                        jump_betauyr * datar.wei.weir[0][8]);

          // Approximate IY+2, right FP
          iy += 1;
          Operator_weights_y(beta, vdey, ix, iy, iz, dy);
          // IY-1 cross right interface, use middle FP
          src[ix][iy][iz] += vdey[0] * (jump_ul * datar.wei.weil[0][5] +
                                        jump_betauyl * datar.wei.weil[0][6] +
                                        jump_ur * datar.wei.weil[0][7] +
                                        jump_betauyr * datar.wei.weil[0][8]);

          ip += 2;
        }
      }
    }
  }

  // Step III: Apply MIB to source term in Z-direction
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
          src[ix][iy][iz] += vdez[2] * (jump_u * data.wei.weir[0][4] +
                                        jump_betauz * data.wei.weir[0][5]);

          // Approximate IZ+1
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross interface, use left FP
          src[ix][iy][iz] += vdez[0] * (jump_u * data.wei.weil[0][4] +
                                        jump_betauz * data.wei.weil[0][5]);

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
          src[ix][iy][iz] += vdez[2] * (jump_ul * datal.wei.weir[0][5] +
                                        jump_betauzl * datal.wei.weir[0][6] +
                                        jump_ur * datal.wei.weir[0][7] +
                                        jump_betauzr * datal.wei.weir[0][8]);

          // Approximate IZ+1, the inside corner point
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross left interface, use left interface's left FP
          src[ix][iy][iz] += vdez[0] * (jump_ul * datal.wei.weil[0][5] +
                                        jump_betauzl * datal.wei.weil[0][6] +
                                        jump_ur * datal.wei.weil[0][7] +
                                        jump_betauzr * datal.wei.weil[0][8]);
          // IZ+1 cross right interface, use right interface's right FP
          src[ix][iy][iz] += vdez[2] * (jump_ul * datar.wei.weir[0][5] +
                                        jump_betauzl * datar.wei.weir[0][6] +
                                        jump_ur * datar.wei.weir[0][7] +
                                        jump_betauzr * datar.wei.weir[0][8]);

          // Approximate IZ+2, right FP
          iz += 1;
          Operator_weights_z(beta, vdez, ix, iy, iz, dz);
          // IZ-1 cross right interface, use middle FP
          src[ix][iy][iz] += vdez[0] * (jump_ul * datar.wei.weil[0][5] +
                                        jump_betauzl * datar.wei.weil[0][6] +
                                        jump_ur * datar.wei.weil[0][7] +
                                        jump_betauzr * datar.wei.weil[0][8]);

          ip += 2;
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
void Douglas_ADI::Operator_weights_x(Beta& beta, VecDoub_O& vdex, Int_I ix,
                                     Int_I iy, Int_I iz, Doub_I dx) {
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
void Douglas_ADI::Operator_weights_y(Beta& beta, VecDoub_O& vdey, Int_I ix,
                                     Int_I iy, Int_I iz, Doub_I dy) {
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
void Douglas_ADI::Operator_weights_z(Beta& beta, VecDoub_O& vdez, Int_I ix,
                                     Int_I iy, Int_I iz, Doub_I dz) {
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
void Douglas_ADI::Convert2Tri_irr(Int_I ix, MatrixDoub_I& row, VecDoub_O& a,
                                  VecDoub_O& b, VecDoub_O& c, VecDoub_O& r) {
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

/************************************************************************************************
 Fast algorithm to convert rectrangular matrix to banded matrix

 INPUT
 lhs : left hand side of equation.  Matrix size: n*2n
 rhs : right hand side of equation. Vector size: n
 ************************************************************************************************/
void Douglas_ADI::Diagonal(MatrixDoub_I& lhs, VecDoub_I& rhs) {
  int length, width;
  double ratio;

  width = (int)lhs.size();
  length = (int)lhs[0].size();

  if (width != (int)rhs.size()) {
    cout << "Left hand side and right hand side should be same size!" << endl;
    exit(0);
  }

  if ((width % 2 != 0) || ((width * 2) != length)) {
    cout << "Bad size for matrix in Diagonal Subroutine!" << endl;
    exit(0);
  }

  // Diagonal matrix
  // make lower triangular zero
  for (int k = 0; k < width - 1; k++) {
    for (int i = k + 1; i < width; i++) {
      ratio = -lhs[i][k] / lhs[k][k];
      for (int j = 0; j < length; j++) {
        lhs[i][j] += lhs[k][j] * ratio;
      }
      rhs[i] += rhs[k] * ratio;
    }
  }
  // make upper triangular zero
  for (int k = 0; k < width - 1; k++) {
    for (int i = width - 2 - k; i >= 0; i--) {
      ratio = -lhs[i][length - 1 - k] / lhs[width - 1 - k][length - 1 - k];
      for (int j = 0; j < length; j++) {
        lhs[i][j] += lhs[width - 1 - k][j] * ratio;
      }
      rhs[i] += rhs[width - 1 - k] * ratio;
    }
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
void Douglas_ADI::Convert2Tri_cor(Int_I ix, MatrixDoub_I& row, VecDoub_O& a,
                                  VecDoub_O& b, VecDoub_O& c, VecDoub_O& r) {
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

/************************************************************************
 Initialization of analytical solution

 INPUT
 eq : equation object
 uh : three-dimensional analytical solution
 ************************************************************************/
void Douglas_ADI::Initialization(Equation& eq, CubicDoub& uh) {
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
void Douglas_ADI::Error(Equation& eq, CubicDoub& uh, ofstream& out_file) {
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
void Douglas_ADI::TDMA(VecDoub_I& a, VecDoub_I& b, VecDoub_I& c, VecDoub_I& r,
                       VecDoub_O& u) {
  int n = (int)a.size();
  double bet;
  VecDoub gam(n);  // One vector of workspace, gam, is needed.

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
      throw("Error 2 in tridag");  // Algorithm fails
    }
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }
  for (int j = (n - 2); j >= 0; j--)
    u[j] -= gam[j + 1] * u[j + 1];  // Backsubstitution.
}

/*************************************************************
 Pentadiagonal matrix fast solver

 INPUT
 a : first hypotenuse of pentadiagonal matrix, a[0],a[1] = 0
 b : second hypotenuse of pentadiagonal matrix, b[0] = 0
 c : third hypotenuse of pentadiagonal matrix
 d : fourth hypotenuse of pentadiagonal matrix, d[n-1] = 0
 e : fifth hypotenuse of pentadiagonal matrix, d[n-2],d[n-1] = 0
 r : right hand side vector

 OUTPUT
 u : solution of pentadiagonal matrix system, Ax = b
 *************************************************************/
void Douglas_ADI::PDMA(VecDoub_I& a, VecDoub_I& b, VecDoub_I& c, VecDoub_I& d,
                       VecDoub_I& e, VecDoub_I& r, VecDoub_O& u) {
  int n = (int)a.size();
  VecDoub p, q;
  double bet, den;

  p.resize(n);
  q.resize(n);

  if (c[0] == 0) {
    cout << "Eliminate u2 trivially" << endl;
    exit(1);
  }

  bet = 1.0 / c[0];
  p[0] = -d[0] * bet;
  q[0] = -e[0] * bet;
  u[0] = r[0] * bet;

  bet = c[1] + b[1] * p[0];
  if (bet == 0) {
    cout << "Singular 1 in PDMA" << endl;
    exit(1);
  }
  bet = -1.0 / bet;
  p[1] = (d[1] + b[1] * q[0]) * bet;
  q[1] = e[1] * bet;
  u[1] = (b[1] * u[0] - r[1]) * bet;

  for (int i = 2; i < n; i++) {
    bet = b[i] + a[i] * p[i - 2];
    den = c[i] + a[i] * q[i - 2] + bet * p[i - 1];
    if (den == 0) {
      cout << "Singular 2 in PDMA" << endl;
      exit(1);
    }
    den = -1.0 / den;
    p[i] = (d[i] + bet * q[i - 1]) * den;
    q[i] = e[i] * den;
    u[i] = (a[i] * u[i - 2] + bet * u[i - 1] - r[i]) * den;
  }

  u[n - 2] = u[n - 2] + p[n - 2] * u[n - 1];
  for (int i = n - 3; i >= 0; i--) {
    u[i] = u[i] + p[i] * u[i + 1] + q[i] * u[i + 2];
  }
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
void Douglas_ADI::Weights(Doub_I z, VecDoub_I& x, Int_I n, Int_I m,
                          MatrixDoub_O& c) {
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
int Douglas_ADI::To1d(Int_I ix, Int_I iy, Int_I iz) {
  int temp;

  temp = (iz * nx * ny) + (iy * nx) + ix;

  return temp;
}