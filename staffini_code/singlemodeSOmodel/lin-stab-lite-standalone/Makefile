HOSTNAME := $(shell hostname)

a2ps=a2pscolor


include Makefile.$(HOSTNAME)


#DEBUG=-O3 -fopenmp
DEBUG=-fbounds-check

WARN=-Wall 

LDFLAGS= $(DEBUG) $(NAGDIR)  $(LAPACKDIR) $(ARCHLDFLAGS)
FFLAGS= $(NAGINC)  $(DEBUG) $(WARN) $(ARCHFFLAGS)
LOADLIBES=$(NAGLIB)

FC=gfortran
LD=gfortran

PROGRAMS=find_transition

all: $(PROGRAMS)
	@echo "DONE"

TEMPLATE=
TEMPLATE: $(TEMPLATE)
	$(LD) $(LDFLAGS) $(TEMPLATE) -o $@ $(LOADLIBES) 

find_transition= global.o overlaps.o linear_stability.o
find_transition: $(find_transition)
	$(LD) $(LDFLAGS) $(find_transition) -o $@ $(LOADLIBES) 

print: 
	$(a2ps) $(patsubst %.o, %.f90, $(find_transition))


%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

%.mod: %.o

clean:
	rm -f *.o *.mod $(PROGRAMS)

realclean: clean
	rm -f *~
