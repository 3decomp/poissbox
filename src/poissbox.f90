!!! src/poissbox.f90
!
! The poissbox library.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

module poissbox

  use mpi
  use petsc

  use constants
  
  implicit none

  integer, public :: nproc ! Number of ranks
  integer, public :: irank ! My rank

  !! DMDA grid
  type(tDM), public :: da

  !! Linear system
  type(tMat), public :: M
  type(tVec), public :: x
  type(tVec), public :: b

  private
  public :: initialise_grid
  public :: initialise_linear_system
  
contains
  
  subroutine initialise_grid(nglobal, da)

    integer, dimension(3), intent(in) :: nglobal
    type(tDM), intent(out) :: da

    integer :: ierr
    
    !! Create a DMDA grid
    call DMDACreate3d(PETSC_COMM_WORLD, &
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
                      DMDA_STENCIL_BOX, &
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
  
end module poissbox
