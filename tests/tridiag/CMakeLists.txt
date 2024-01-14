### tests/tridiag/CMakeLists.txt
#
# Builds the test cases for tri-diagonal solvers.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_library(test_tdma_utils OBJECT test_tdma_utils.f90)

add_executable(test_tdma_sweeps test_tdma_sweeps.f90)
target_link_libraries(test_tdma_sweeps poissbox)
target_link_libraries(test_tdma_sweeps test_tdma_utils)
add_test(NAME test_tdma_sweeps COMMAND test_tdma_sweeps)

add_executable(test_tdma test_tdma.f90)
target_link_libraries(test_tdma poissbox)
target_link_libraries(test_tdma test_tdma_utils)
add_test(NAME test_tdma COMMAND test_tdma)
