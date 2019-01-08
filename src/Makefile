#Makefile to build mcpopsim
# -- macros
CC= mpicc
CXXFLAGS= -O3 -std=c++11 -m64
OBJECTS= deme.o allele_dat.o event.o locus.o

# -- targets
all: mcpopsim
	mcpopsim: $(OBJECTS)
	$(CC) -o mcpopsim $(OBJECTS)

deme.o: deme.c
	$(CC) $(CXXFLAGS) -c deme.c

events.o: events.c
	$(CC) $(CXXFLAGS) -c events.c

allele_dat.o: allele_dat.c
	$(CC) $(CXXFLAGS) -c allele_dat.c

locus.o: locus.c
	$(CC) $(CXXFLAGS) -c locus.c
# -- remove stuff
clean:
	 rm -f mcpopsim $(OBJECTS)
