FC=mpif90

OPT+=-DPID

OPTIONS=$(OPT)

#FFTWLIB=/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/
#FFTWINC=/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/

FFTWINC=/usr/local/include/
FFTWLIB=/usr/local/lib/

XFLAG=-O3 -cpp -mcmodel=medium -fcoarray=single -fopenmp -lfftw3f -I$(FFTWINC) -L$(FFTWLIB) $(OPTIONS)

#XFLAG=-O3 -fpp -shared-intel -openmp -mcmodel=medium -traceback -xHost -coarray=distributed
OFLAG=$(XFLAG) -c


MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o)))

all: cafcube.x
	@echo "done"
cafcube.x: $(OBJFILE)
	@echo $(value OBJFILE)
	$(FC) $(XFLAG) $(OBJFILE) -o $@  -lfftw3f -lm -ldl

## For addtional subroutines, add its dependency on variables.o, and dependency of cafcube.o on it.

variables.o: parameters.o
initialize.o: variables.o cubefft.o penfft.o redshifts.txt
finalize.o: variables.o cubefft.o penfft.o
cafcube.o: variables.o cubefft.o penfft.o
cubefft.o: variables.o
penfft.o: variables.o
kernel_f.o: variables.o
kernel_c.o: variables.o
timestep.o: variables.o
particle_initialization.o: variables.o
buffer_density.o: variables.o
buffer_x.o: variables.o
buffer_v.o: variables.o
pm.o: variables.o cubefft.o penfft.o
update_particle.o: variables.o
checkpoint.o: variables.o
projection.o: variables.o
cafcube.o: parameters.o variables.o initialize.o cubefft.o penfft.o kernel_f.o kernel_c.o timestep.o particle_initialization.o buffer_density.o buffer_x.o buffer_v.o pm.o update_particle.o checkpoint.o projection.o finalize.o

%.o: %.f90
	@echo "Files to be linked:"
	$(FC) $(OFLAG) $< -o $@

clean:
	rm -f *.mod *.o *.out *.x *~
