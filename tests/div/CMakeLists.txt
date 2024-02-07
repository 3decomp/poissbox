### tests/div/CMakeLists.txt
#
# Builds the test cases for divergence evaluation.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_executable(test_div_1d test_div_1d.f90)
target_link_libraries(test_div_1d poissbox)
add_test(NAME test_div_1d COMMAND test_div_1d)

add_executable(test_div_3d test_div_3d.f90)
target_link_libraries(test_div_3d poissbox)
add_test(NAME test_div_3d COMMAND test_div_3d)
