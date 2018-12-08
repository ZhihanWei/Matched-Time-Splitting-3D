# MIB + Fast Time Splitting Algorithm

## Project Name

Algorithm based on finite difference for parabolic interface problem in 3D

## Project Description

A parabolic PDE solver for interface problem in 3D space based on finite difference method. The solver handles interface jump conditions by involving fictitious points. A new developed numerical method - Matched Interface and Boundary method has been applied to smooth the solutions along the interfaces. Taking the advantage of the tridiagonal matrix structure and TDMA algorithm, the solver is able to achieve unconditional stable for various complicated interfaces and jump conditions. Meanwhile, the high dimensional problem are deducted to several one dimensional problems and achieves a flop count of O(N). This solver contains two spacial methods: MIB-L1 and MIB-L2 and four temporal methods: Douglas-ADI, LOD-IE, LOD-CN, Trapezoidal Splitting. 

## Build
```
mkdir build && cd build && cmake ..
```
Default build in release mode, to build with debug, add `-DCMAKE_BUILD_TYPE=debug`, binary `Matched-TS` is in `build/bin`

## Reminds

• Compiled with C++14 support. 

## Classes

• Diffusion coefficients : "Beta"
 
• Analytical solution : "Eq"

• Implicit function of surface in a Cartesian coordinate : "Surface_Cartesian"

• Mesh constraction : "Mesh"

• MIB : "Intersection"

• Temporal methods : "ADI", "LOD", "TS"

• Decomposition : "LU"












