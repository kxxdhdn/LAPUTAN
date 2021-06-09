#*****************************************************************************
#*
#*                        Makefile for the programs
#*
#*****************************************************************************

# Compiler
FC = h5fc
FFLAGS = -Wall -O3 -g # Standard (mono-processor)
# FFLAGS = -Wall -O3 -g -fcheck=all # debug
#FFLAGS = -Wall -O3 -ftree-parallelize-loops=4 # Mac (parallelization)
#FFLAGS = -warn all -O3 -heap-arrays -traceback # irfucoast intel compiler
#FFLAGS = -warn all -check all -g -debug -heap-arrays -traceback # debug
LDFLAGS = $(FFLAGS)

# Modules
MODS = utilities.f90 \
       constants.f90 \
       arrays.f90 \
       inout.f90 \
       linear_system.f90 \
       interpolation.f90 \
       integration.f90 \
       adaptative_grid.f90 \
       instrument_filters.f90

# Name of the library and its paths
LIB = swing
LIBAR = lib$(LIB).a

# Banner
define BANNER

*******************************************************************************
                      MAKING SwING LIBRARY snippet
*******************************************************************************

endef
$(info $(BANNER))

# Programs
SRCS = synthetic_photometry.f90
OBJS = $(MODS:.f90=.o)
EXEC = $(SRCS:.f90=)

# Generate all executables
all : $(EXEC) lib

# Make the library archive
lib : $(LIBAR)
$(LIBAR) : $(OBJS)
	@echo "\nCREATING THE STATIC LIBRARY...\n"
	ar cr $(LIBAR) $^

# Linking
% : $(OBJS) %.o
	$(FC) $(LDFLAGS) -o $@ $^

# Individual file compilation
%.o : ../*/%.f90
	$(FC) $(FFLAGS) -c $?
.PRECIOUS: %.o # keep intermediate object files

# House cleaning
clean :
	@echo "CLEANING THE LIBRARY...\n"
	@rm -f *.o
mrproper : clean
	@echo "MR. PROPERING THE LIBRARY...\n"
	@rm -f *.mod $(EXEC) $(LIBAR)
