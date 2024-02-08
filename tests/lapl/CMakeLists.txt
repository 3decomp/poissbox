### tests/lapl/CMakeLists.txt
#
# Builds the test cases for Laplacian evaluation.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_executable(test_lapl_ho test_lapl.f90)
target_link_libraries(test_lapl_ho poissbox)
add_test(NAME test_lapl_ho COMMAND test_lapl_ho)
