FC = gfortran
objects = help.o test.o full.o lowrank.o matprod.o makeGFull.o solveIntFull.o plotField.o hmatrices.o wrongArg.o
modules = helpmod.mod testmod.mod fullmod.mod lowrankmod.mod matprodmod.mod makeGFullmod.mod solveIntFullmod.mod plotFieldmod.mod wrongArgmod.mod

# Eventuele compilatievlaggen per compiler
FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
# FFLAGS_gfortran = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -g -debug full -implicitnone -check -warn -free -Tf
FFLAGS_nagfor   = -g -C=all -gline -u -info -colour -kind=byte

# Selecteer de juiste vlaggen voor de huidige compiler
FFLAGS=$(FFLAGS_$(FC))

hmatrices: $(objects)
	$(FC) -o hmatrices $(options) $(objects)

hmatrices.o: hmatrices.f90 $(modules)
	$(FC) -c $(FFLAGS) hmatrices.f90
help.o: help.f90
	$(FC) -c $(FFLAGS) help.f90
helpmod.mod: help.f90 help.o
	$(FC) -c $(FFLAGS) help.f90
test.o: test.f90
	$(FC) -c $(FFLAGS) test.f90
testmod.mod: test.f90 test.o
	$(FC) -c $(FFLAGS) test.f90
full.o: full.f90
	$(FC) -c $(FFLAGS) full.f90
fullmod.mod: full.f90 full.o
	$(FC) -c $(FFLAGS) full.f90
lowrank.o: lowrank.f90
	$(FC) -c $(FFLAGS) lowrank.f90
lowrankmod.mod: lowrank.f90 lowrank.o
	$(FC) -c $(FFLAGS) lowrank.f90
matprod.o: matprod.f90
	$(FC) -c $(FFLAGS) matprod.f90
matprodmod.mod: matprod.f90 matprod.o
	$(FC) -c $(FFLAGS) matprod.f90
makeGFull.o: makeGFull.f90
	$(FC) -c $(FFLAGS) makeGFull.f90
makeGFullmod.mod: makeGFull.f90 makeGFull.o
	$(FC) -c $(FFLAGS) makeGFull.f90
solveIntFull.o: solveIntFull.f90
	$(FC) -c $(FFLAGS) solveIntFull.f90
solveIntFullmod.mod: solveIntFull.f90 solveIntFull.o
	$(FC) -c $(FFLAGS) solveIntFull.f90
plotField.o: plotField.f90
	$(FC) -c $(FFLAGS) plotField.f90
plotFieldmod.mod: plotField.f90 plotField.o
	$(FC) -c $(FFLAGS) plotField.f90
wrongArg.o: wrongArg.f90
	$(FC) -c $(FFLAGS) wrongArg.f90
wrongArgmod.mod: wrongArg.f90 wrongArg.o
	$(FC) -c $(FFLAGS) wrongArg.f90

clean:
	rm $(objects) $(modules)
