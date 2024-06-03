add_library(shared_module "")

target_sources(shared_module 
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/data_mods.f90
        ${CMAKE_CURRENT_LIST_DIR}/input_module.f90
        ${CMAKE_CURRENT_LIST_DIR}/command_line_args.f90
        ${CMAKE_CURRENT_LIST_DIR}/spectral_density.f90
)