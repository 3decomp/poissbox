### tests/grad/CMakeLists.txt
#
# Builds the test cases for gradient evaluation.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

add_executable(test_grad_1d test_grad_1d.f90)
target_link_libraries(test_grad_1d poissbox)
add_test(NAME test_grad_1d COMMAND test_grad_1d)

add_executable(test_grad_3d test_grad_3d.f90)
target_link_libraries(test_grad_3d poissbox)
add_test(NAME test_grad_3d COMMAND test_grad_3d)
