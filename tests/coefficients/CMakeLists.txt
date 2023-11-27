### tests/coefficients/CMakeLists.txt
#
# Builds the test cases for coefficients.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_executable(test_d2dx2 test_d2dx2.f90)
add_test(NAME test_d2dx2 COMMAND test_d2dx2)
