FC := pgf90
FCFLAGS = -fast -O3 
FLFLAGS = -lm -L/panfs/pfs.local/software/7/install/fftw3/3.3.8/pgi/19/lib -lfftw3

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

PROGRAM = sfg_fluc_serial
OBJECTS = data_mods.f90 spectra_fluc.f90 read_input.f90 read_field.f90 read_eners.f90 tcf.f90 spec_calc.f90 hist_calc.f90

$(PROGRAM): $(OBJECTS)
	$(FC) $(FCFLAGS) $(FLFLAGS) $(OBJECTS) -o  $(PROGRAM) 


