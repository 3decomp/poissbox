!!! src/coefficients.f90
!
! Coefficients module of Poissbox, provides subroutines for computing coefficients of the Laplacian
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

module coefficients

  implicit none

  private
  public :: lapl_1d_coeffs
  public :: lapl_star_coeffs
  
contains

  pure function lapl_1d_coeffs(dx) result(coeffs)

    real, intent(in) :: dx
    real, dimension(3) :: coeffs

    real :: invdx2
    
    invdx2 = 1.0 / dx**2
    
    coeffs(1) = invdx2
    coeffs(2) = -2.0 * invdx2
    coeffs(3) = invdx2
    
  end function lapl_1d_coeffs

  ! Compute the Laplacian coefficients for a 7-point "star", returning a 3x3x3 box of coefficients.
  pure function lapl_star_coeffs(dx, dy, dz) result(coeffs)

    real, intent(in) :: dx, dy, dz
    real, dimension(3, 3, 3) :: coeffs

    coeffs(:, :, :) = 0.0
    coeffs(:, 2, 2) = coeffs(:, 2, 2) + lapl_1d_coeffs(dx)
    coeffs(2, :, 2) = coeffs(2, :, 2) + lapl_1d_coeffs(dy)
    coeffs(2, 2, :) = coeffs(2, 2, :) + lapl_1d_coeffs(dz)
    
  end function lapl_star_coeffs
  
end module coefficients
