### tests/CMakeLists.txt
#
# poissbox test driver.
#
## License
#
# SPDX-License-Identifier: BSD-3-Clause

# Coefficient tests
add_subdirectory(coefficients)

# Tridiagonal solver tests
add_subdirectory(tridiag)
