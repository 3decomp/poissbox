!!! tests/coefficients/test_d2dx2.f90
!
! Test computing the coefficients for the 1-D Laplacian operator d2/dx2.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program test_d2dx2

  implicit none

  real, dimension(3), parameter :: fc = [ 3.14, 3.14, 3.14 ] ! Constant field
  real, dimension(3), parameter :: fg = [ 1.0, 2.0, 3.0 ]    ! Constant gradient field
  real, dimension(3), parameter :: fq = fg**2                ! Quadratic field
  
  ! Laplacian of a constant field is 0
  ! Invariant to scaling
  ! Invariant to shifting

  ! Laplacian of a constant gradient is 0
  ! Invariant to scaling
  ! Invariant to shifting

  ! Laplacian of ax^2 -> 2ax
  ! Proportional to scaling
  ! Invariant to shifting

  ! Laplacian of ax^2 + bx + c -> 2ax
  ! Proportional to scaling of the quadratic term
  ! Invariant to shifting

  stop 1
  
end program test_d2dx2
