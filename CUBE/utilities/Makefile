OPTIONS=-DFFTFINE
#OPTIONS=-DFFTHALF
#OPTIONS+=-Wall
#OPTIONS+=-DPID

all: ic.x dsp.x cicpower.x qspace.x corr.x potential.x Jmatrix.x detJ.x #convert.x
	@echo "done"
ic.x: parameters.o pencil_fft.o initial_conditions.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
dsp.x: parameters.o pencil_fft.o powerspectrum.o displacement.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
cicpower.x: parameters.o pencil_fft.o powerspectrum.o cicpower.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
corr.x: parameters.o pencil_fft.o powerspectrum.o correlation_qspace.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
potential.x: parameters.o pencil_fft.o potential_qspace.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
Jmatrix.x: parameters.o pencil_fft.o powerspectrum.o Jmatrix_qspace.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
convert.x: parameters.o convert_zip_to_xv.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
qspace.x: parameters.o xva_qspace.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
detJ.x: parameters.o detJ.o
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
xva_qspace.o: xva_qspace.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
correlation_qspace.o: correlation_qspace.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
potential_qspace.o: potential_qspace.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
Jmatrix_qspace.o: Jmatrix_qspace.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
cicpower.o: cicpower.f90 parameters.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
convert_zip_to_xv.o: convert_zip_to_xv.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
detJ.o: detJ.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
*.o: Makefile

clean:
	rm -f *.mod *.o *.out *.x *~
