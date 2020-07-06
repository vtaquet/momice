MODUL=modules

#-- Building Rules
.SUFFIXES:
.SUFFIXES: .f90 .o
.SUFFIXES: .f .o
.SUFFIXES: .c .o
.f.o :
	$(FC) $(FLAG90) $(INC90) -c $*.f
.f90.o :
	$(FC) $(FLAG90) $(INC90) -c $*.f90

 SRC = variables.f90 \
	intlib.f90 \
	allocation.f90 \
	input.f90 \
	opkda1.f \
	opkda2.f \
	opkdmain.f \
	physprop.f90 \
	chemprop.f90 \
	odes.f90 \
	numprop.f90 \
	output.f90 \
	main.f90

 OBJ = variables.o \
	intlib.o \
	allocation.o \
	input.o \
	opkda1.o \
	opkda2.o \
	opkdmain.o \
	physprop.o \
	chemprop.o \
	odes.o \
	numprop.o \
	output.o \
	main.o

 OBJ2 = variables.o \
	allocation.o \
	input.o \
	physprop.o \
	chemprop.o \
	numprop.o \
	odes.o \
	output.o \
	main.o

main:
	make `hostname -s`

ifort:
	make momice \
	"OUT=momice.ifort" \
	"FC=ifort" \
	"FLAG90=-debug -O3" \
	"INC90  = -I$(MODUL)" \
	"LIB90= "

gfort:
	make momice \
	"OUT = momice.gfort" \
	"FC=gfortran" \
	"FLAG90 = -fbounds-check -J$(MODUL) -O3" \
	"LIB90 = "
	#"LIB90  = -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial"  -debug -fstack-arrays 


f2py:
	make momice \
	"OUT = momice.f2py" \
	"FC=f2py" \
	"FLAG90 = " \
	"LIB90 = "
	#"LIB90  = -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial"  -debug -fstack-arrays 

momice: $(OBJ)
	$(FC) $(FLAG90) $(OBJ) $(INC90) $(LIB90) -o $(OUT) #-pg

cleanall:
	\rm -f *.o $(MODUL)/*.mod

clean: 
	\rm -f $(OBJ2) $(MODUL)/$(OBJ2).mod
