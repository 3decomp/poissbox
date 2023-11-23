# poissbox

Testing solving the Poisson equation using iterative solvers

## Building poissbox

`poissbox` is written in modern Fortran and depends on `MPI` and `PETSc`.
A `CMake` build system is used, this should pick up your `MPI` compilers automatically, however to
find `PETSc` you should set `PKG_CONFIG_PATH` to point to the `PETSc` package configuration
directory, *e.g.*
```
export PKG_CONFIG_PATH=${PETSC_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH}
```
The `poissbox` build can then be configured and run by executing
```
cmake -B build .
cmake --build build/
```
which should build the executable `build/bin/poissbox`.

The executable should run, printing how many ranks it is running on, followed by a "Hello world!"
message from each rank, and the distribution of DoF between ranks (currently hardcoded as 64x64x64)
*e.g.*
```
$ mpirun -np 3 ./build/bin/poissbox
 Running poissbox on            3  ranks
 Hello from            0
 Hello from            1
 Hello from            2
 Rank            0  has        90112  of       262144  expected:       262144
 Rank            1  has        86016  of       262144  expected:       262144
 Rank            2  has        86016  of       262144  expected:       262144
```
this confirms that no errors occur when initialising the `PETSc` data structures.
