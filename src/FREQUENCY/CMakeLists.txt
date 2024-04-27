add_executable(calc_freq freq_calc.f90)

target_link_libraries(calc_freq
    shared_module
    ${HDF5_LIBRARIES}
)