# This is a make file for my EVB_AddOn code
# Start of makefile
#=====================================================================

SHELL=/bin/sh

.SUFFIXES:
.SUFFIXES: .f90 .o

BINROOT=.
EX=PSO.Z
EXE=$(BINROOT)/$(EX)

TYPE=master

FC=undefined
LD=undefined

# Define object files
#=====================================================================

OBJ_ALL = \
	Soln.o \
	\
	Startup.o \
	functions.o \
	\
	velocity.o \
	\
	PSOTest.o
	
# Examine targets manually
#=====================================================================
all:
	@echo
	@echo
	@echo "--------------------------"
	@echo
	@echo

# Clean up the source directory
#=====================================================================

clean:
	rm -f $(OBJ_ALL) *.mod

#====================== Generic f95 compilers - DEBUG ================
jake:
	$(MAKE) LD="gfortran -o" \
	LDFLAGS="-O3" \
	FC="gfortran -c" \
	FCFLAGS="-O3" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

# Declare rules
#=====================================================================

# Master
#=====================================================================

master: message $(OBJ_ALL)
	$(LD) $(EXE) $(LDFLAGS) $(OBJ_ALL)

# Message
message:
	@echo "PSO compiler"
	@echo
	@echo


.f90.o:
	$(FC) $(FCFLAGS) $*.f90
