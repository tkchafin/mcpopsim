#Makefile to build mcpopsim
# -- macros
CC= mpicc
CXX=mpicxx
LDFLAGS=-g
CXXFLAGS= -O3 -std=c++11 -m64
OBJECTS= mcpopsim.o deme.o allele_dat.o event.o locus.o
#LDLIBS=

# -- targets
#mcpopsim: mcpopsim.cpp mcpopsim.h $(OBJECTS)
#	$(CXX) $(CXXFLAGS) -c mcpopsim.cpp

all: mcpopsim allele_dat.o deme.o event.o locus.o mcpopsim.o

mcpopsim: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o mcpopsim $(OBJECTS) $(LDLIBS)

mcpopsim.o: mcpopsim.cpp
	$(CXX) $(CXXFLAGS) -c mcpopsim.cpp

allele_dat.o: allele_dat.cpp allele_dat.h locus.h
	$(CXX) $(CXXFLAGS) -c allele_dat.cpp

deme.o: deme.cpp deme.h allele_dat.h locus.h
	$(CXX) $(CXXFLAGS) -c deme.cpp

event.o: event.cpp event.h
	$(CXX) $(CXXFLAGS) -c event.cpp

locus.o: locus.cpp locus.h
	$(CXX) $(CXXFLAGS) -c locus.cpp

# -- remove stuff
clean:
	 rm -f mcpopsim $(OBJECTS)
