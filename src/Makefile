DEPENDENCIES =  mainLoop.o plasmaSlabField.o disp_rel_fluid.o
LOOPPPATH = ../Loopp/loopp.a
HEADERSPATH =../Loopp/headers
OUTPUTPATH = ../run

libs = -lm
CXXFLAGS =  -O3 -std=c++14 #-Wall
CXX = g++
MPIXX = mpic++

all: plasmaSlabSim

plasmaSlabSim : $(DEPENDENCIES)
	$(CXX) -I$(HEADERSPATH) -o $@ $^ $(LOOPPPATH)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -I$(HEADERSPATH) -c $< 

clean:
	$(RM) *.o
	$(RM) plasmaSlabSim
