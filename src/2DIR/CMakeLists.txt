add_executable(calc_2dir_spectrum 2dir_calc.f90 read_field.f90 spec_calc.f90 tcf.f90)

target_link_libraries(calc_2dir_spectrum
    fftw3
    shared_module
    ${HDF5_LIBRARIES}
)