### CMakeLists.txt
#
# Main build configuration file for poissbox.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

## Basic CMake configuration
cmake_minimum_required(VERSION 3.10)
project(poissbox
  LANGUAGES Fortran)
set(version 0.0)
set(CMAKE_Fortran_PREPROCESS ON)

set(CMAKE_Fortran_FLAGS "-std=f2018 -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -fcheck=all -fimplicit-none -Wall -Wpedantic -Werror -ffpe-trap=invalid,zero,overflow")

include(GNUInstallDirs)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

## Find packages
find_package(PkgConfig REQUIRED)

# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

# PETSc - recommended to
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)

## Build poissbox
add_library(poissbox
  src/constants.f90
  src/coefficients.f90
  src/tridsol.f90
  src/compact_schemes.f90
  src/poissbox.f90)

target_link_libraries(poissbox MPI::MPI_Fortran)
target_link_libraries(poissbox PkgConfig::PETSC)

add_executable(poissbox_demo src/example.f90)
target_link_libraries(poissbox_demo poissbox)
add_executable(compare_fd_compact src/compare_fd_compact.f90)
target_link_libraries(compare_fd_compact poissbox)

## Testing
include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
