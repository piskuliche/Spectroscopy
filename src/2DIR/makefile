FC := gfortran
FCFLAGS = -O3 -fopenmp -g -fbounds-check
FLFLAGS = -lm -I/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/gcc-10.2.0/hdf5-1.10.6-ofyjlzowxfz3you35yhc2qd4wv5wwgj4/include -L/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/gcc-10.2.0/hdf5-1.10.6-ofyjlzowxfz3you35yhc2qd4wv5wwgj4/lib -lhdf5_fortran -L/cm/shared/apps/spack/gpu/opt/spack/linux-centos8-skylake_avx512/gcc-8.3.1/fftw-3.3.8-pef6n347azkie3trtz4mliy33kltj7my/lib/ -lfftw3

%.o: %.f90
        $(FC) $(FCFLAGS) -c $<

PROGRAM = 2dir_fluc_serial
OBJECTS = data_mods.f90 2dir_fluc.f90 read_input.f90 read_field.f90 read_eners.f90 tcf.f90 spec_calc.f90

$(PROGRAM): $(OBJECTS)
        $(FC) $(FCFLAGS) $(FLFLAGS) $(OBJECTS) -o  $(PROGRAM)