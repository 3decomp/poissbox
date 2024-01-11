!!! src/example.f90
!
! The poissbox example program.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program poissbox_example

  use mpi
  use petsc

  use constants
  use coefficients
  use poissbox
  use matfree_types
  
  implicit none

  integer :: ierr  ! PETSc error code

  !! Grid dimensions - hardcoded for simplicity
  integer, parameter :: nx = 64
  integer, parameter :: ny = 64
  integer, parameter :: nz = 64
  integer, dimension(3), parameter :: n = [nx, ny, nz]

  !! Problem dimension - hardcoded for simplicity
  real(pb_dp), parameter :: Lx = 1.0_pb_dp
  real(pb_dp), parameter :: Ly = 1.0_pb_dp
  real(pb_dp), parameter :: Lz = 1.0_pb_dp
  real(pb_dp), parameter :: dx = Lx / nx
  real(pb_dp), parameter :: dy = Ly / ny
  real(pb_dp), parameter :: dz = Lz / nz

  !! Local variables
  real(pb_dp) :: error
  type(tVec) :: x2
  type(mat_ctx) :: ctx
  
  !! Initialise MPI & PETSc
  call MPI_Init(ierr) ! Could rely on PetscInitialize
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)

  if (irank == 0) then
     print *, "Running poissbox on ", nproc, " ranks"
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, "Hello from ", irank

  call initialise_grid(n, da)
  call check_grid(n, da)

  ctx%da = da
  ctx%grid_deltas = [dx, dy, dz]
  if (.true.) then ! Set true to use matrix-free system
     call initialise_linear_system(da, ctx, P, A, x, b)
  else
     call initialise_linear_system(da, ctx, P, P, x, b)
     A = P
  end if
  call check_linear_system(n, P, x, b)
  call assemble_laplacian(da, dx, dy, dz, P)

  ! Set a solution, apply the linear operator and check its pointwise approximation of the Laplacian
  call set_solution(da, x)
  print *, "Calling MatMult"
  call MatMult(A, x, b, ierr)
  call check_lapl(da, x, b)

  call check_matrices(A, P, x)

  ! Solve the linear system and compare against the specified solution
  call solve(P, A, x, b)
  call VecDuplicate(x, x2, ierr)
  print *, "Calling MatMult"
  call MatMult(A, x, x2, ierr)
  call VecAXPY(x2, -1.0d0, b, ierr)
  call VecNorm(x2, NORM_2, error, ierr)
  print *, "Solution residual (L2 norm): ", error
  
  !! Finalise MPI & PETSc
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)

contains

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

  subroutine set_solution(da, x)
    !! Sets an initial value in the solution vector

    type(tDM), intent(in) :: da
    type(tVec), intent(inout) :: x

    integer :: istart, jstart, kstart
    integer :: ni, nj, nk

    integer :: i, j, k

    integer :: ierr

    real(pb_dp), dimension(:, :, :), pointer :: xdof
    real(pb_dp) :: x_p ! Value at DoF P
    real(pb_dp) :: xsum, xsum_v ! Vector sum, used to check the vector contents

    xsum = 0.0d0

    call DMDAGetCorners(da, istart, jstart, kstart, ni, nj, nk, ierr)

    call DMDAVecGetArrayF90(da, x, xdof, ierr)
    
    do k = kstart, nk - 1
       do j = jstart, nj - 1
          do i = istart, ni - 1
             call random_number(x_p) ! Get random number 0 <= x <= 1
             x_p = 2 * (0.5d0 - x_p) ! Rescale to -1 <= x <= 1

             xdof(i, j, k) = x_p
             
             xsum = xsum + x_p ! Compute vector sum
          end do
       end do
    end do

    call DMDAVecRestoreArrayF90(da, x, xdof, ierr)
    call VecAssemblyBegin(x, ierr)
    call VecAssemblyEnd(x, ierr)
    
    ! Check vector contents using vector sum
    call VecSum(x, xsum_v, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, xsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    print *, "Rank ", irank, "Delta of XSUM norms computed directly and from X: ", xsum_v - xsum, xsum_v, xsum

  end subroutine set_solution

  subroutine check_lapl(da, x, b)
    !! Checks the Laplacian calculation computed by matrix-vector product vs the pointwise
    !! computation using the stencil operator (the two should be identical).

    use compute_lapl
    
    type(tDM), intent(in) :: da
    type(tVec), intent(in) :: x ! The (specified) solution
    type(tVec), intent(in) :: b ! The Laplacian approximation, computed as b = Mx

    type(tVec) :: b2 ! Copy of the Laplacian approximation
    type(tVec) :: c  ! The Laplacian approximation, computed pointwise
    real(pb_dp) :: residual ! The residual norm between expected and computed Laplacian fields.

    integer :: ierr

    call VecDuplicate(b, b2, ierr)
    call VecCopy(b, b2, ierr)
    
    call VecDuplicate(b, c, ierr)
    call compute_lapl_pointwise(da, [dx, dy, dz], x, c)
    
    ! Compute the norm of the difference between computed Laplacian, and the result computed
    ! pointwise
    call VecAXPY(b2, -1.0d0, c, ierr) ! Computes b = alpha * c + b, alpha = -1
    call VecNorm(b2, NORM_2, residual, ierr)

    print *, "Rank ", irank, "Delta between b=Mx and pointwise calculation: ", residual

    call VecDestroy(b2, ierr)
    call VecDestroy(c, ierr)
    
  end subroutine check_lapl

  subroutine check_matrices(P, A, x)

    type(tMat), intent(in) :: P ! Preconditioner matrix
    type(tMat), intent(in) :: A ! System matrix
    type(tVec), intent(in) :: x ! Solution vector

    type(tVec) :: bP, bA ! Vectors for computing matvecs

    real(pb_dp) :: delta

    integer :: ierr

    call VecDuplicate(x, bP, ierr)
    call VecDuplicate(x, bA, ierr)

    call MatMult(P, x, bP, ierr)
    call MatMult(A, x, bA, ierr)

    call VecAXPY(bA, -1.0d0, bP, ierr)
    call VecNorm(bA, NORM_2, delta, ierr)

    print *, "Ax - Px = ", delta
    
    call VecDestroy(bP, ierr)
    call VecDestroy(bA, ierr)
    
  end subroutine check_matrices
  
end program poissbox_example
