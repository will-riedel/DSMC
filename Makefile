######################################################################
####### Novice users should modify entries in this block only ########
######################################################################

# FORTRAN_SRC1 = Main.F90 Readin.F90 Matrix_Assem+Solver.F90 Equation_Library.F90 Equation_Library_Explicit.F90 Visualization.F90
FORTRAN_SRC1 = Main.f90 Initialize.f90 Save_Data.f90 Run_Collisions.f90 Initialize_Source.f90 Reflection_Handling.f90 Update_Cell_Index.f90 Read_Input_Data.f90 Cylindrical_Functions.f90 Randn.f90

EXC1 = Main

######################################################################
#### More advanced users may want to modify entries in this block ####
######################################################################

# FC = mpif90
FC = gfortran
# FC = ifort

# FCFLAGS = -I$(HPC_PETSC_DIR)/include

# for gfortran: -Ofast optimizes for speed, -Og optimizes for debugging
# FCFLAGS = 
# gfortran
# FCFLAGS = -O3 -fbacktrace -Wall
FCFLAGS = -O3
# ifort
# FCFLAGS = -fast -assume buffered_io -xCORE-AVX2
# FCFLAGS = -fast -assume buffered_io
# FCFLAGS = -fast

# LDFLAGS = -shared-intel
LDFLAGS = 

PETSC_LIBS = -L$(HPC_PETSC_LIB) -lpetsc

# LIBS = $(PETSC_LIBS)
LIBS = 

F77 = ifort

######################################################################
################## Experts only beyond this point!! ##################
######################################################################

OBJF1 = $(FORTRAN_SRC1:.f=.o)

OBJ1 = $(OBJF1)

##BINDIR = /home/acm/dinesh/fin

######################################################################
############################## Targets ###############################
######################################################################

$(EXC1): $(OBJ1)
	$(FC) $(LDFLAGS) $(FCFLAGS) -o $(EXC1) $(OBJ1) $(LIBS)

depend: $(FORTRAN_SRC1)
	$(BINDIR)/fmakedepend $(FORTRAN_SRC1)

.f.o:  
	$(FC) $(FCFLAGS) $(LIBS) $< -o $(<:.f=.o)

.c.o:
	$(CC) $(CFLAGS) $<

clean: 
	rm *.o

distclean:
	rm -f main *.o

######################################################################
## Automatic dependencies from make depend follow ##
######################################################################

## DO NOT ADD, MODIFY, OR DELETE THIS LINE -- Make depends on it ##

