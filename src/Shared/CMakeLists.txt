add_library(shared_module "")

target_sources(shared_module 
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/data_mods.f90
        ${CMAKE_CURRENT_LIST_DIR}/mapping_module.f90
)