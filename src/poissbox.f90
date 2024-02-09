!!! src/poissbox.f90
!
! The poissbox library.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

module matfree_types

  use petsc

  use constants

  implicit none

  type :: mat_ctx
     type(tDM) :: da
     real(pb_dp), dimension(3) :: grid_deltas
  end type mat_ctx
  
end module matfree_types

module matfree
  
  implicit none

  private

  public :: MatCreateShell
  public :: MatShellGetContext
  public :: MatShellSetContext
  
  interface MatCreateShell
     subroutine MatCreateShell(comm, nrow_l, ncol_l, nrow_g, ncol_g, ctx, M, ierr)
       use petsc
       use matfree_types, only : mat_ctx
       integer :: comm
       integer :: nrow_l ! Local number of rows
       integer :: ncol_l ! Local number of columns
       integer :: nrow_g ! Global number of rows
       integer :: ncol_g ! Global number of columns
       type(mat_ctx) :: ctx ! The shell matrix context
       type(tMat) :: M   ! The matrix object
       integer :: ierr
     end subroutine MatCreateShell
  end interface MatCreateShell

  interface MatShellSetContext
     subroutine MatShellSetContext(M, ctx, ierr)
       use petsc
       use matfree_types, only : mat_ctx
       type(tMat) :: M      ! The matrix object
       type(mat_ctx) :: ctx ! The shell matrix context
       integer :: ierr
     end subroutine MatShellSetContext
  end interface MatShellSetContext

  interface MatShellGetContext
     subroutine MatShellGetContext(M, ctx, ierr)
       use petsc
       use matfree_types, only : mat_ctx
       type(tMat) :: M               ! The matrix object
       type(mat_ctx), pointer :: ctx ! The shell matrix context
       integer :: ierr
     end subroutine MatShellGetContext
  end interface MatShellGetContext
  
end module matfree

module compute_lapl

  use petsc

  use constants
  
  implicit none

  private
  public :: evaluate_laplacian_pointwise
  public :: compute_lapl_pointwise
  public :: compute_lapl_compact
  
contains
  
  subroutine compute_lapl_pointwise(da, grid_deltas, x, b)
    !! Apply the Laplacian stencil to the given solution pointwise - this should be identical to
    !! using the matrix-vector product.

    type(tDM), intent(in) :: da
    real(pb_dp), dimension(3), intent(in) :: grid_deltas
    type(tVec), intent(in) :: x ! The (specified) solution
    type(tVec), intent(in) :: b ! The Laplacian approximation, computed as b_i = stencil_op(x, i)

    integer :: istart, jstart, kstart
    integer :: ni, nj, nk

    integer :: i, j, k
    
    integer :: ierr

    type(tVec) :: xlocal ! Ghosted x
    real(pb_dp), dimension(:, :, :), pointer :: xdof
    real(pb_dp), dimension(:, :, :), pointer :: bdof

    call DMGetLocalVector(da, xlocal, ierr)
    call DMGlobalToLocal(da, x, INSERT_VALUES, xlocal, ierr)
    
    call DMDAGetCorners(da, istart, jstart, kstart, ni, nj, nk, ierr)

    call DMDAVecGetArrayF90(da, xlocal, xdof, ierr)
    call DMDAVecGetArrayF90(da, b, bdof, ierr)

    do k = kstart, kstart + nk - 1
       do j = jstart, jstart + nj - 1
          do i = istart, istart + ni - 1
             bdof(i, j, k) = evaluate_laplacian_pointwise(xdof(i-1:i+1, j-1:j+1, k-1:k+1), &
                                                          grid_deltas)
          end do
       end do
    end do

    call DMDAVecRestoreArrayF90(da, xlocal, xdof, ierr)
    call DMDAVecRestoreArrayF90(da, b, bdof, ierr)
    
    call DMRestoreLocalVector(da, xlocal, ierr)
    
  end subroutine compute_lapl_pointwise

  subroutine compute_lapl_compact(da, x, dx, b)

    use compact_schemes, only : lapl
    
    type(tDM), intent(in) :: da
    type(tVec), intent(in) :: x     ! The (current) solution
    real(pb_dp), dimension(3) :: dx ! The grid spacing
    type(tVec), intent(in) :: b     ! The computed Laplacian

    integer :: istart, jstart, kstart
    integer :: ni, nj, nk
    
    integer :: ierr

    type(tVec) :: xlocal ! Ghosted x
    real(pb_dp), dimension(:, :, :), pointer :: xdof
    real(pb_dp), dimension(:, :, :), pointer :: bdof

    call DMGetLocalVector(da, xlocal, ierr)
    call DMGlobalToLocal(da, x, INSERT_VALUES, xlocal, ierr)
    
    call DMDAGetCorners(da, istart, jstart, kstart, ni, nj, nk, ierr)

    call DMDAVecGetArrayF90(da, xlocal, xdof, ierr)
    call DMDAVecGetArrayF90(da, b, bdof, ierr)
    
    call lapl(xdof(istart:istart + (ni - 1),jstart:jstart + (nj - 1),kstart:kstart + (nk - 1)), &
              dx, &
              bdof(istart:istart + (ni - 1),jstart:jstart + (nj - 1),kstart:kstart + (nk - 1)))
    
    call DMDAVecRestoreArrayF90(da, xlocal, xdof, ierr)
    call DMDAVecRestoreArrayF90(da, b, bdof, ierr)
    
    call DMRestoreLocalVector(da, xlocal, ierr)

  end subroutine compute_lapl_compact
  
  real(pb_dp) pure function evaluate_laplacian_pointwise(f, grid_deltas)
    !! Applies the Laplacian stencil at a point.
    
    use coefficients, only: lapl_star_coeffs
    
    real(pb_dp), dimension(3, 3, 3), intent(in) :: f     ! The field within the stencil region
    real(pb_dp), dimension(3), intent(in) :: grid_deltas ! The X,Y,Z grid spacing

    real(pb_dp) :: dx, dy, dz
    
    real(pb_dp), dimension(3, 3, 3) :: coeffs

    dx = grid_deltas(1)
    dy = grid_deltas(2)
    dz = grid_deltas(3)
    coeffs = lapl_star_coeffs(dx, dy, dz)

    evaluate_laplacian_pointwise = dot_product(reshape(f, [size(f)]), &
                                               reshape(coeffs, [size(coeffs)]))

  end function evaluate_laplacian_pointwise

end module compute_lapl

module poissbox

  use mpi
  use petsc

  use constants
  use compute_lapl
  
  implicit none

  integer, public :: nproc ! Number of ranks
  integer, public :: irank ! My rank

  !! DMDA grid
  type(tDM), public :: da

  !! Linear system
  type(tMat), public :: P ! The preconditioner matrix
  type(tMat), public :: A ! The system matrix
  type(tVec), public :: x
  type(tVec), public :: b

  ! external mfmult
  
  private
  public :: initialise_grid
  public :: initialise_linear_system
  public :: solve
  
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

  subroutine initialise_linear_system(da, ctx, P, A, x, b)
    !! Given a grid object, build the linear system matrix, solution and RHS vectors.

    use matfree_types

    type(tDM), intent(in) :: da
    type(mat_ctx), intent(in) :: ctx
    type(tMat), intent(inout) :: P ! Preconditioner linear system
    type(tMat), intent(inout) :: A ! Main linear system
    type(tVec), intent(out) :: x ! Solution vector
    type(tVec), intent(out) :: b ! RHS vector

    integer :: ierr

    print *, "Initialising linear system"
    
    call DMCreateMatrix(da, P, ierr)
    call MatSetFromOptions(P, ierr)
    call MatSetUp(P, ierr)

    if (A /= P) then
       call initialise_matrix_free(ctx, P, A)
    end if
    
    call DMCreateGlobalVector(da, x, ierr)
    call VecSetFromOptions(x, ierr)
    call VecSetUp(x, ierr)

    call DMCreateGlobalVector(da, b, ierr)
    call VecSetFromOptions(b, ierr)
    call VecSetUp(b, ierr)

    print *, "Done"
    
  end subroutine initialise_linear_system

  subroutine initialise_matrix_free(ctx, P, A)
    !! Create a matrix free object
    
    use matfree_types
    use matfree
    
    type(mat_ctx) :: ctx
    type(tMat), intent(in) :: P
    type(tMat), intent(out) :: A


    integer :: m, n
    
    integer :: ierr

    print *, "- Initialising matrix-free system"
    
    call MatGetLocalSize(P, m, n, ierr)
    
    call MatCreateShell(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE, ctx, A, ierr)
    call MatShellSetContext(A, ctx, ierr) ! Is this necessary?
    call MatShellSetOperation(A, MATOP_MULT, mfmult, ierr)
    
    print *, "- Done"
    
  end subroutine initialise_matrix_free
  
  subroutine solve(P, A, x, b)
    !! Solve the linear system.
    !
    !! XXX: Currently hardcoded to treat periodic/Neumann (singular) problems.
    
    type(tMat), intent(in) :: P    ! The preconditioner matrix
    type(tMat), intent(in) :: A    ! The system matrix
    type(tVec), intent(inout) :: x ! The solution vector
    type(tVec), intent(in) :: b    ! The RHS vector

    type(tKSP) :: ksp
    type(tMatNullSpace) :: nsp

    integer :: ierr

    ! XXX: Periodic/Neumann problems require a null-space to remove the constant solution
    call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_VEC, nsp, ierr)
    call MatSetNullSpace(P, nsp, ierr)
    if (A /= P) then
       print *, "Setting nullspace on A"
       call MatSetNullSpace(A, nsp, ierr)
    end if
    call MatNullSpaceDestroy(nsp, ierr)

    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, A, P, ierr)
    call KSPSetFromOptions(ksp, ierr)
    call KSPSolve(ksp, b, x, ierr)
    
  end subroutine solve

  subroutine mfmult(M, x, f, ierr)
    !! Computes the matrix vector product f = Mx, matrix-free

    use matfree_types
    use matfree

    use compute_lapl
  
    type(tMat) :: M ! The operator
    type(tVec) :: x ! The input vector
    type(tVec) :: f ! The output vector
    integer :: ierr ! Error status (0 indicates success)

    type(mat_ctx), pointer :: ctx
    
    call MatShellGetContext(M, ctx, ierr)
    if (.true.) then
       call compute_lapl_compact(ctx%da, x, ctx%grid_deltas, f)
    else
       call compute_lapl_pointwise(ctx%da, ctx%grid_deltas, x, f)
    end if
    
    ! associate(foo => M)
    ! end associate
    ! call MatMult(P, x, f, ierr)

  end subroutine mfmult
  
end module poissbox
