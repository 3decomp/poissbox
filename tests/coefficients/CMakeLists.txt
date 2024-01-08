### tests/coefficients/CMakeLists.txt
#
# Builds the test cases for coefficients.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_executable(test_d2dx2 test_d2dx2.f90)
target_link_libraries(test_d2dx2 poissbox)
add_test(NAME test_d2dx2 COMMAND test_d2dx2)

add_executable(test_star test_star.f90)
target_link_libraries(test_star poissbox)
add_test(NAME test_star COMMAND test_star)

add_executable(test_lapl test_lapl.f90)
target_link_libraries(test_lapl poissbox)
add_test(NAME test_lapl COMMAND test_lapl)
