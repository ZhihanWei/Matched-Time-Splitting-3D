# ----- makefile to compile the program heat2dll ----- #

#++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#               compilation environment                #
#   HOST: Mac,Windows,dmc; OMPILER: g++,clang++        #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#COMPILER  = clang++
COMPILER	 = g++
CFLAGS     = -c -Wall -std=c++11 -O3
VPATH      = ./src
SOURCES    = main.cpp Data.cpp Mesh.cpp Intersections.cpp Douglas_ADI.cpp Surface_Cartesian.cpp \
						 Surface_Cube.cpp Surface_Cylinder.cpp Surface_Ellipsoid.cpp Surface_Pile.cpp Surface_Torus.cpp \
						 Surface_Cone.cpp Surface_Dupin_Cyclide.cpp Surface_Molecular.cpp Surface_Heart.cpp Surface_Tanglecube.cpp \
						 Equation.cpp Eq_0.cpp Eq_1.cpp Eq_2.cpp Eq_3.cpp Eq_4.cpp Eq_5.cpp Eq_6.cpp Eq_7.cpp LU.cpp
OBJECTS    = $(SOURCES:.cpp=.o)
EXECUTABLE = 3D_MADI

#+++++++++++++++++ compilation options ++++++++++++++++#
$(EXECUTABLE): $(OBJECTS)
	$(COMPILER) $(OBJECTS) -o $(EXECUTABLE)
	@echo "     "
	@echo ">>> compiled on `hostname -s` with $(COMPILER) <<<"

%.o: %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

run:
	./$(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)