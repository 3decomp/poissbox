!!! src/coefficients.f90
!
! Coefficients module of Poissbox, provides subroutines for computing coefficients of the Laplacian
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

module coefficients

  use constants
  
  implicit none

  private
  public :: lapl_1d_coeffs
  public :: lapl_star_coeffs
  public :: assemble_laplacian
  
contains

  pure function lapl_1d_coeffs(dx) result(coeffs)

    real(pb_dp), intent(in) :: dx
    real(pb_dp), dimension(3) :: coeffs

    real(pb_dp) :: invdx2
    
    invdx2 = 1.0_pb_dp / dx**2
    
    coeffs(1) = invdx2
    coeffs(2) = -2.0 * invdx2
    coeffs(3) = invdx2
    
  end function lapl_1d_coeffs

  ! Compute the Laplacian coefficients for a 7-point "star", returning a 3x3x3 box of coefficients.
  pure function lapl_star_coeffs(dx, dy, dz) result(coeffs)

    real(pb_dp), intent(in) :: dx, dy, dz
    real(pb_dp), dimension(3, 3, 3) :: coeffs

    coeffs(:, :, :) = 0.0
    coeffs(:, 2, 2) = coeffs(:, 2, 2) + lapl_1d_coeffs(dx)
    coeffs(2, :, 2) = coeffs(2, :, 2) + lapl_1d_coeffs(dy)
    coeffs(2, 2, :) = coeffs(2, 2, :) + lapl_1d_coeffs(dz)
    
  end function lapl_star_coeffs

  subroutine assemble_laplacian(da, dx, dy, dz, M)
    !! Builds the Laplacian operator matrix.
#include "petsc/finclude/petscmat.h"

    use petsc

    type(tDM), intent(in) :: da
    real(pb_dp), intent(in) :: dx, dy, dz
    type(tMat), intent(inout) :: M
    
    real(pb_dp), dimension(3, 3, 3) :: coeffs
    real(pb_dp), dimension(27) :: mat_coeffs
    integer :: coeff_ctr
    MatStencil :: row(4, 1)
    MatStencil :: col(4, 27)
    
    integer :: istart, jstart, kstart
    integer :: ni, nj, nk

    integer :: i, j, k
    integer :: ii, jj, kk
    
    integer :: ierr
    
    call DMDAGetCorners(da, istart, jstart, kstart, ni, nj, nk, ierr)

    do k = kstart, kstart + nk - 1
       do j = jstart, jstart + nj - 1
          do i = istart, istart + ni - 1
             mat_coeffs(:) = 0.0_pb_dp
             
             ! Get coefficients at point
             coeffs = lapl_star_coeffs(dx, dy, dz)

             row(MatStencil_i, 1) = i
             row(MatStencil_j, 1) = j
             row(MatStencil_k, 1) = k

             ! Flatten coefficients and compute stencil indices
             coeff_ctr = 1
             do kk = 1, 3
                do jj = 1, 3
                   do ii = 1, 3
                      mat_coeffs(coeff_ctr) = coeffs(ii, jj, kk)

                      col(MatStencil_i, coeff_ctr) = i + (ii - 2)
                      col(MatStencil_j, coeff_ctr) = j + (jj - 2)
                      col(MatStencil_k, coeff_ctr) = k + (kk - 2)

                      coeff_ctr = coeff_ctr + 1
                   end do
                end do
             end do

             ! Push coefficients to matrix
             call MatSetValuesStencil(M, 1, row, 27, col, mat_coeffs, INSERT_VALUES, ierr)
          end do
       end do
    end do

    call MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr)
    
  end subroutine assemble_laplacian
  
end module coefficients
