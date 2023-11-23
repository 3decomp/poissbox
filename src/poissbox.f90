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

  !! Finalise MPI & PETSc
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)
  
end program poissbox
