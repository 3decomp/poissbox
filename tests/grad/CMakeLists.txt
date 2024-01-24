### tests/grad/CMakeLists.txt
#
# Builds the test cases for gradient evaluation.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_executable(test_grad test_grad.f90)
target_link_libraries(test_grad poissbox)
add_test(NAME test_grad COMMAND test_grad)
