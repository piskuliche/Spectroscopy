add_executable(calc_ir_spectrum ir_calc.f90 read_field.f90 read_input.f90 spec_calc.f90 spectral_histogram.f90 tcf.f90)

target_link_libraries(calc_ir_spectrum
    shared_module
    fftw3
    ${HDF5_LIBRARIES}
)