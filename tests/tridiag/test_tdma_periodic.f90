program test_tdma_periodic

  use constants

  use tridsol

  use test_tdma_utils
  
  implicit none

  integer, parameter :: n = 128 ! Problem size
  real(pb_dp), dimension(:), allocatable :: a, b, c, x, d ! Tridiagonal parameters

  logical :: passing

  passing = .true.

  ! Test periodic solver on a periodic system
  call init(n, a, b, c, x, d, opt_periodic=.true.)
  call check_tdma_periodic(a, b, c, x, d)

  ! Test periodic solver on a non-periodic system - this should still work, it just adds cost over
  ! the non-periodic solver.
  call init(n, a, b, c, x, d, opt_periodic=.false.)
  call check_tdma_periodic(a, b, c, x, d)

  call fin()

contains

  subroutine fin()

    deallocate(a, b, c, x, d)

    if (.not. passing) then
       stop 1
    end if
    
  end subroutine fin

  subroutine check_tdma_periodic(a, b, c, x, d)

    real(pb_dp), dimension(:), allocatable, intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: b ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: x ! Solution
    real(pb_dp), dimension(:), allocatable, intent(in) :: d ! RHS

    real(pb_dp), dimension(:), allocatable :: bprime
    real(pb_dp), dimension(:), allocatable :: dprime

    real(pb_dp) :: errrms
    real(pb_dp) :: tol
    
    bprime = b
    dprime = d
    
    call tdma_periodic(a, bprime, c, dprime)

    tol = epsilon(errrms)
    errrms = sum((x - dprime)**2)
    errrms = sqrt(errrms / n)
    if (errrms > (tol * sqrt(sum(x**2) / n))) then
       print *, "Periodic TDMA failed: ", errrms
       passing = .false.
    else
       print *, "Periodic TDMA passed: ", errrms
    end if
    
  end subroutine check_tdma_periodic
  
end program test_tdma_periodic
