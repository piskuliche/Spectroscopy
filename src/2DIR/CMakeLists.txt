add_executable(calc_2dir_spectrum 2ir_fluc.f90 read_input.f90 read_field.f90 spec_calc.f90 read_eners.f90 tcf.f90)

target_link_libraries(calc_2dir_spectrum
    fftw3
    ${HDF5_LIBRARIES}
)