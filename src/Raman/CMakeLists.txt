add_executable(calc_raman_spectrum raman_calc.f90 spec_calc.f90 tcf.f90)

target_link_libraries(calc_raman_spectrum
    shared_module
    fftw3
    ${HDF5_LIBRARIES}
)