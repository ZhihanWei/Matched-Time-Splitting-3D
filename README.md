# README

## Project Name

Algorithm based on finite difference for parabolic interface problem in 3D

## Project Description

A parabolic PDE solver for interface problem in 3D space based on finite difference method. The solver handles interface jump conditions by involving fictitious points. A new developed numerical method - Matched Interface and Boundary method has been applied to smooth the solutions along the interfaces. Taking the advantage of the tridiagonal matrix structure and TDMA algorithm, the solver is able to achieve unconditional stable for various complicated interfaces and jump conditions. Meanwhile, the high dimensional problem are deducted to several one dimensional problems and achieves a flop count of O(N). This solver contains two spacial methods: MIB-L1 and MIB-L2 and four temporal methods: Douglas-ADI, LOD-IE, LOD-CN, Trapezoidal Splitting. 

## Implementation

Step 1: In terminal(Mac)/Command Line(Windows), go to the folder where the makefile is located, and type "make", it will compile and link to an executable file.

Step 2: Once you see ">>> compiled on (hostname of your PC) with  <<<", this mean the executable file is generated successfully.

Step 3: Type "./1D_MIB_TD" to run the executable file, you'll see results on the screen.

Step 4: If you want to compile with a set of new parameters, type "make clean" to remove the previous object files and executable file. Then go to Step 1.

## Reminds

• The solver has to be compiled by the compiler(like GNU ur LLVM) with C++11 support.

• Created a folder called "result" before running which will contain all the results in txt files.  

## Classes

• Variable coefficients : "Beta"
  1. If extra Beta class is added, include header in main.cpp
  2. Add it to line 164 - 179 in main.cpp
 
• Analytical solution : "Eq"
	1. If extra Beta class is added, include header in main.cpp
	2. Add it to line 183 - 216 in main.cpp

• Implicit function of surface in a Cartesian coordinate : "Surface", "Surface_Cartesian"
	1. If extra Beta class is added, include header in main.cpp
	2. Add it to all ADI, LOD, TS starting and solver file

• Mesh constraction : "Mesh"

• Fictitious points and jump approximation : "Intersections"

• Temporal methods : "ADI", "LOD", "TS"

• LU decomposition : "LU"

## User-specified data

### Data.txt

• xl,xr,yl,yr,zl,zr: leftmost and rightmost domains for x-, y-, z- directions, refer to line 64-79 in "main.cpp" for details for different surface.

• nx,ny,nz: numbers of nodes in x-, y-, z- directions. 

• t_start: strating time, usually set to be 0.

• t_finish: ending time.

• t_step: time step.

• surface: letter for surface refer to line 64-79 in "main.cpp".

• method: A: Douglas-ADI; L: LOD(LOD-IE/LOD-CN); T: Trapezoidal Splitting.

• equation: number for analytical solutions, refer to "Eq_<number>.cpp"s.

• beta: number for diffusion coefficient, refer to "Beta_<number>.cpp"s.

• mib_method: spacial method, 1: MIB-L1; 2: MIB-L2.

• accuracy: equal to 2 always

### Constant.h

• NPRINT: numbers of print out results for one calculation. 

• REG: jump approximation region. i: only inside; others: both inside and outside. 

• JP: jump in temporal method. r: real jump: a: approximated jump;  

### main.cpp

line of 70 - 90  : parameters to determine shape of surfaces 

line of 99 - 102 : switch from running 1 input file to several input files. 















