add_executable(calc_ffcf ffcf_calc.f90 ffcf.f90 read_field.f90)

target_link_libraries(calc_ffcf
    shared_module
    ${HDF5_LIBRARIES}
)