FC = ifort

src = ./src

objects = help.o test.o full.o lowrank.o matprod.o makeGFull.o makeGHmat.o solveIntFull.o plotField.o vecProdHmat.o hmatrices.o matrixConverter.o
modules = helpmod.mod testmod.mod fullmod.mod lowrankmod.mod matprodmod.mod makegfullmod.mod makeghmatmod.mod solveintfullmod.mod plotfieldmod.mod vecprodhmatmod.mod matrixconverter.mod

# Eventuele compilatievlaggen per compiler
FFLAGS_g95      = -O3
FFLAGS_gfortran = -O3 #-g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -O3 #-g -debug full -implicitnone -check -warn -free -Tf
FFLAGS_nagfor   = -O3 #-g -C=all -gline -u -info -colour -kind=byte

options = -lblas -llapack 

# Selecteer de juiste vlaggen voor de huidige compiler
FFLAGS=$(FFLAGS_$(FC))

hmatrices: $(objects)
	$(FC) -o hmatrices $(objects) $(options)
#	make clean

hmatrices.o: $(src)/hmatrices.f90 $(modules)
	$(FC) -c $(FFLAGS) $(src)/hmatrices.f90

help.o: $(src)/help.f90
	$(FC) -c $(FFLAGS) $(src)/help.f90
helpmod.mod: $(src)/help.f90 help.o
	$(FC) -c $(FFLAGS) $(src)/help.f90

test.o: $(src)/test.f90 matrixconverter.mod fullmod.mod solveintfullmod.mod makegfullmod.mod lowrankmod.mod matprodmod.mod plotfieldmod.mod vecprodhmatmod.mod makeghmatmod.mod
	$(FC) -c $(FFLAGS) $(src)/test.f90
testmod.mod: $(src)/test.f90 test.o matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/test.f90

full.o: $(src)/full.f90 matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/full.f90
fullmod.mod: $(src)/full.f90 full.o
	$(FC) -c $(FFLAGS) $(src)/full.f90

lowrank.o: $(src)/lowrank.f90 matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/lowrank.f90
lowrankmod.mod: $(src)/lowrank.f90 lowrank.o
	$(FC) -c $(FFLAGS) $(src)/lowrank.f90

matprod.o: $(src)/matprod.f90 matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/matprod.f90
matprodmod.mod: $(src)/matprod.f90 matprod.o matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/matprod.f90

makeGFull.o: $(src)/makeGFull.f90 matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/makeGFull.f90
makegfullmod.mod: $(src)/makeGFull.f90 makeGFull.o
	$(FC) -c $(FFLAGS) $(src)/makeGFull.f90

makeGHmat.o: $(src)/makeGHmat.f90 matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/makeGHmat.f90
makeghmatmod.mod: $(src)/makeGHmat.f90 makeGHmat.o
	$(FC) -c $(FFLAGS) $(src)/makeGHmat.f90

solveIntFull.o: $(src)/solveIntFull.f90 matrixconverter.mod makegfullmod.mod
	$(FC) -c $(FFLAGS) $(src)/solveIntFull.f90
solveintfullmod.mod: $(src)/solveIntFull.f90 solveIntFull.o
	$(FC) -c $(FFLAGS) $(src)/solveIntFull.f90

plotField.o: $(src)/plotField.f90 matrixconverter.mod
	$(FC) -c $(FFLAGS) $(src)/plotField.f90
plotfieldmod.mod: $(src)/plotField.f90 plotField.o
	$(FC) -c $(FFLAGS) $(src)/plotField.f90

matrixConverter.o: $(src)/matrixConverter.f90
	$(FC) -c $(FFLAGS) $(src)/matrixConverter.f90
matrixconverter.mod: $(src)/matrixConverter.f90 matrixConverter.o
	$(FC) -c $(FFLAGS) $(src)/matrixConverter.f90

vecProdHmat.o: $(src)/vecProdHmat.f90 matrixconverter.mod makeghmatmod.mod
	$(FC) -c $(FFLAGS) $(src)/vecProdHmat.f90
vecprodhmatmod.mod: $(src)/vecProdHmat.f90 vecProdHmat.o
	$(FC) -c $(FFLAGS) $(src)/vecProdHmat.f90

# clean:
# 	rm $(objects) $(modules)
