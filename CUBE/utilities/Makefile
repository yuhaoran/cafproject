OPTIONS=-DFFTFINE
#OPTIONS+=-Wall
#OPTIONS+=-DPID

all: ic.x dsp.x cicpower.x #convert.x
	@echo "done"
ic.x: parameters.o pencil_fft.o initial_conditions.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
dsp.x: parameters.o pencil_fft.o powerspectrum.o displacement.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
cicpower.x: parameters.o pencil_fft.o powerspectrum.o cicpower.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
convert.x: parameters.o convert_zip_to_xv.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@


parameters.o: ../main/parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
pencil_fft.o: ../main/pencil_fft.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
powerspectrum.o: powerspectrum.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
initial_conditions.o: initial_conditions.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
displacement.o: displacement.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
cicpower.o: cicpower.f90 parameters.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
convert_zip_to_xv.o: convert_zip_to_xv.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
*.o: Makefile

clean:
	rm -f *.mod *.o *.out *.x *~
