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
  public :: tdma, tdma_periodic
  public :: fwd_sweep ! Exported for testing purposes 
  public :: bwd_sweep ! Exported for testing purposes
  
contains
  
  subroutine tdma(a, b, c, d)

    real(pb_dp), dimension(:), intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), intent(inout) :: b ! Super-diagonal
    real(pb_dp), dimension(:), intent(in) :: c ! Diagonal
    real(pb_dp), dimension(:), intent(inout) :: d ! RHS/solution

    call fwd_sweep(a, b, c, d)
    call bwd_sweep(b, c, d)
    
  end subroutine tdma
  
  subroutine tdma_periodic(a, b, c, d)

    real(pb_dp), dimension(:), intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), intent(inout) :: b ! Super-diagonal
    real(pb_dp), dimension(:), intent(in) :: c ! Diagonal
    real(pb_dp), dimension(:), intent(inout) :: d ! RHS/solution

    integer :: n

    ! Auxilliary systems
    real(pb_dp) :: gamma
    real(pb_dp), dimension(:), allocatable :: bmod
    real(pb_dp), dimension(:), allocatable :: u

    n = size(a)
    
    ! Select gamma to increase diagonal dominance
    gamma = -b(1)

    ! Solve auxilliary systems
    bmod = b
    bmod(1) = bmod(1) - gamma
    bmod(n) = bmod(n) - c(n) * a(1) / gamma
    call tdma(a, bmod, c, d)

    bmod = b
    bmod(1) = bmod(1) - gamma
    bmod(n) = bmod(n) - c(n) * a(1) / gamma
    allocate(u(n))
    u(:) = 0.0_pb_dp
    u(1) = gamma
    u(n) = c(n)
    call tdma(a, bmod, c, u)

    ! Assemble solution
    d(:) = d(:) - (u(:) * (d(1) + (a(1) / gamma) * d(n))) &
         / (1.0_pb_dp + (u(1) + (a(1) / gamma) * u(n)))

    deallocate(u)
    
  end subroutine tdma_periodic

  subroutine fwd_sweep(a, b, c, d)

    real(pb_dp), dimension(:), intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), intent(inout) :: b ! Super-diagonal
    real(pb_dp), dimension(:), intent(in) :: c ! Diagonal
    real(pb_dp), dimension(:), intent(inout) :: d ! RHS

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

    real(pb_dp), dimension(:), intent(in) :: b    ! Diagonal
    real(pb_dp), dimension(:), intent(in) :: c    ! Super-diagonal
    real(pb_dp), dimension(:), intent(inout) :: d ! Solution

    integer :: i

    integer :: n

    n = size(d)
    
    d(n) = d(n) / b(n)
    do i = n - 1, 1, -1
       d(i) = (d(i) - c(i) * d(i + 1)) / b(i)
    end do
    
  end subroutine bwd_sweep
  
end module tridsol
