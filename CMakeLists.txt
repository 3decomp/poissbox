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

include(GNUInstallDirs)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

## Find packages
find_package(PkgConfig REQUIRED)

# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

# PETSc - recommended to
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)

add_executable(poissbox src/poissbox.f90)
target_link_libraries(poissbox MPI::MPI_Fortran)
target_link_libraries(poissbox PkgConfig::PETSC)

## Testing
include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
