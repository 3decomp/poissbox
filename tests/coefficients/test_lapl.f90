!!! tests/coefficients/test_star.f90
!
! Test the Laplacian coefficients across the domain by assembling and applying the Laplacian using
! MatVec.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program test_lapl

  implicit none

  real, dimension(:, :, :), allocatable :: f    ! The field
  real :: dx
  real, dimension(:, :, :), allocatable :: lapl ! The computed Laplacian

  stop 1
  call evaluate_laplacian(f, dx, lapl)

contains

  subroutine evaluate_laplacian(f, dx, lapl)

    real, dimension(:, :, :), intent(in) :: f
    real, intent(in) :: dx
    real, dimension(:, :, :), allocatable, intent(out) :: lapl

    associate(foo => dx)
    end associate
    
    ! Initialise the Laplacian result
    allocate(lapl, mold=f)
    lapl(:, :, :) = 0.0

    ! Set Laplacian matrix
    ! Compute Laplacian matvec
    
  end subroutine evaluate_laplacian
  
end program test_lapl
