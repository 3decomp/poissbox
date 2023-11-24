!!! src/poissbox.f90
!
! The poissbox main program.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program poissbox

  use mpi
  use petsc
  
  implicit none

  integer :: nproc ! Number of ranks
  integer :: irank ! My rank
  integer :: ierr  ! PETSc error code

  !! Grid dimensions - hardcoded for simplicity
  integer, parameter :: nx = 64
  integer, parameter :: ny = 64
  integer, parameter :: nz = 64
  integer, dimension(3), parameter :: n = [nx, ny, nz]

  !! DMDA grid
  type(tDM) :: da

  !! Linear system
  type(tMat) :: M
  type(tVec) :: x
  type(tVec) :: b
  
  !! Initialise MPI & PETSc
  call MPI_Init(ierr) ! Could rely on PetscInitialize
  call PetscInitialize(ierr)

  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)

  if (irank == 0) then
     print *, "Running poissbox on ", nproc, " ranks"
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, "Hello from ", irank

  call initialise_grid(n, da)
  call check_grid(n, da)

  call initialise_linear_system(da, M, x, b)
  call check_linear_system(n, M, x, b)
  
  !! Finalise MPI & PETSc
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)

contains

  subroutine initialise_grid(nglobal, da)

    integer, dimension(3), intent(in) :: nglobal
    type(tDM), intent(out) :: da
    
    !! Create a DMDA grid
    call DMDACreate3d(PETSC_COMM_WORLD, &
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
                      DMDA_STENCIL_STAR, &
                      nglobal(1), nglobal(2), nglobal(3), &
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &
                      1, &
                      1, &
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                      da, &
                      ierr)
    call DMSetFromOptions(da, ierr)
    call DMSetUp(da, ierr)

  end subroutine initialise_grid

  subroutine check_grid(nglobal, da)
    !! Basic test: check that the DMDA grid is the correct size, summing the number of DoF on each
    !! rank should give the global DoF count.

    integer, dimension(3), intent(in) :: nglobal
    type(tDM), intent(in) :: da

    integer :: istart, jstart, kstart
    integer :: ni, nj, nk

    integer :: ierr

    integer :: ndof_local
    integer :: ndof_global
    integer :: ndof_expect
    
    call DMDAGetCorners(da, istart, jstart, kstart, ni, nj, nk, ierr)

    ndof_local = ni * nj * nk
    call MPI_Allreduce(ndof_local, ndof_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ndof_expect = product(nglobal)
    print *, "(DMDA): Rank ", irank, " has ", ndof_local, " of ", ndof_global, " expected: ", ndof_expect
    
  end subroutine check_grid

  subroutine initialise_linear_system(da, M, x, b)
    !! Given a grid object, build the linear system matrix, solution and RHS vectors.

    type(tDM), intent(in) :: da
    type(tMat), intent(out) :: M
    type(tVec), intent(out) :: x
    type(tVec), intent(out) :: b

    integer :: ierr

    call DMCreateMatrix(da, M, ierr)
    call MatSetFromOptions(M, ierr)
    call MatSetUp(M, ierr)
    
    call DMCreateGlobalVector(da, x, ierr)
    call VecSetFromOptions(x, ierr)
    call VecSetUp(x, ierr)

    call DMCreateGlobalVector(da, b, ierr)
    call VecSetFromOptions(b, ierr)
    call VecSetUp(b, ierr)
    
  end subroutine initialise_linear_system

  subroutine check_linear_system(nglobal, M, x, b)
    !! Check that the sizes of the linear system components are sensible.

    integer, dimension(3), intent(in) :: nglobal
    type(tMat), intent(in) :: M
    type(tVec), intent(in) :: x
    type(tVec), intent(in) :: b

    integer :: myrow   ! First row index on this rank
    integer :: nextrow ! First row index on next rank

    integer :: ndof_local
    integer :: ndof_global
    integer :: ndof_expect
    
    integer :: ierr

    ndof_expect = product(nglobal)

    call MatGetOwnershipRange(M, myrow, nextrow, ierr)
    ndof_local = nextrow - myrow
    call MPI_Allreduce(ndof_local, ndof_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    print *, "(M): Rank ", irank, " has ", ndof_local, " rows of ", ndof_global, " expected: ", ndof_expect

    call VecGetOwnershipRange(x, myrow, nextrow, ierr)
    ndof_local = nextrow - myrow
    call MPI_Allreduce(ndof_local, ndof_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    print *, "(x): Rank ", irank, " has ", ndof_local, " rows of ", ndof_global, " expected: ", ndof_expect

    call VecGetOwnershipRange(b, myrow, nextrow, ierr)
    ndof_local = nextrow - myrow
    call MPI_Allreduce(ndof_local, ndof_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    print *, "(b): Rank ", irank, " has ", ndof_local, " rows of ", ndof_global, " expected: ", ndof_expect

  end subroutine check_linear_system
  
end program poissbox
