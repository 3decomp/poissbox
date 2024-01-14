!!! src/tridsol.f90
!
! Tridiagonal solver module of Poissbox, implements a basic tridiagonal solver and a periodic
! version.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause
module tridsol

  use constants
  
  implicit none

  private
  public :: tdma
  public :: fwd_sweep ! Exported for testing purposes 
  public :: bwd_sweep ! Exported for testing purposes
  
contains
  
  subroutine tdma(a, b, c, d)

    real(pb_dp), dimension(:), allocatable, intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), allocatable, intent(inout) :: b ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(inout) :: d ! RHS/solution

    call fwd_sweep(a, b, c, d)
    call bwd_sweep(b, c, d)
    
  end subroutine tdma

  subroutine fwd_sweep(a, b, c, d)

    real(pb_dp), dimension(:), allocatable, intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), allocatable, intent(inout) :: b ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(inout) :: d ! RHS

    real(pb_dp) :: w
    integer :: i

    integer :: n

    n = size(d)
    
    do i = 2, n
       w = a(i) / b(i - 1)
       b(i) = b(i) - w * c(i - 1)
       d(i) = d(i) - w * d(i - 1)
    end do

  end subroutine fwd_sweep

  subroutine bwd_sweep(b, c, d)

    real(pb_dp), dimension(:), allocatable, intent(in) :: b    ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c    ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(inout) :: d ! Solution

    integer :: i

    integer :: n

    n = size(d)
    
    d(n) = d(n) / b(n)
    do i = n - 1, 1, -1
       d(i) = (d(i) - c(i) * d(i + 1)) / b(i)
    end do
    
  end subroutine bwd_sweep
  
end module tridsol
