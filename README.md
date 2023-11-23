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
message from each rank, *e.g.*
```
$ mpirun -np 2 ./build/bin/poissbox
 Running poissbox on            2  ranks
 Hello from            0
 Hello from            1
```
this confirms that no errors occur when initialising the `PETSc` data structures.
