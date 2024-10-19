DEBUG=-O3

WARN=-Wall -Wno-unused -Wno-tabs

HOSTNAME := $(shell hostname)

include Makefile.$(HOSTNAME)

FFLAGS= $(NAGINC)  $(DEBUG) $(WARN) 
LOADLIBES=$(NAGLIB)

PROGRAMS=main 

all: $(PROGRAMS)
	@echo "DONE"

main=global.o eom_so.o evolution.o sweeps.o overlaps.o linearstability.o main.o
main: $(main)
	$(LD) $(LDFLAGS) $(main) -o $@ $(LOADLIBES) 	

TEMPLATE=
TEMPLATE: $(TEMPLATE)
	$(LD) $(LDFLAGS) $(TEMPLATE) -o $@ $(LOADLIBES) 	

%.o: %.f95
	$(FC) $(FFLAGS) -c $*.f95

clean:
	rm -f *.o *.mod $(PROGRAMS)

realclean: clean
	rm -f *~
