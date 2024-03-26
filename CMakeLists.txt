cmake_minimum_required(VERSION 3.13)

project(electric-field LANGUAGES Fortran)
find_package(HDF5 REQUIRED COMPONENTS Fortran)

message(STATUS ${HDF5_Fortran_INCLUDE_DIRS})

set(CMAKE_Fortran_FLAGS "-O3")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include(GNUInstallDirs)
include_directories(${HDF5_INCLUDE_DIR})

message(STATUS ${CMAKE_BINARY_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# defines targets and sources
add_subdirectory(src)