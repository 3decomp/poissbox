program test_tdma

  use constants

  use tridsol

  use test_tdma_utils
  
  implicit none

  integer, parameter :: n = 128 ! Problem size
  real(pb_dp), dimension(:), allocatable :: a, b, c, x, d ! Tridiagonal parameters

  logical :: passing

  passing = .true.

  ! Test non-periodic solver on non-periodic system
  call init(n, a, b, c, x, d)
  call check_tdma(a, b, c, x, d, .true.)

  ! Test periodic solver on periodic system (this should fail!)
  call init(n, a, b, c, x, d, opt_periodic=.true.)
  call check_tdma(a, b, c, x, d, .false.)
  
  call fin()

contains

  subroutine fin()

    deallocate(a, b, c, x, d)

    if (.not. passing) then
       stop 1
    end if
    
  end subroutine fin

  subroutine check_tdma(a, b, c, x, d, expect)

    real(pb_dp), dimension(:), allocatable, intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: b ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: x ! Solution
    real(pb_dp), dimension(:), allocatable, intent(in) :: d ! RHS
    logical, intent(in) :: expect ! Expect pass (true) or fail (false)
    
    real(pb_dp), dimension(:), allocatable :: bprime
    real(pb_dp), dimension(:), allocatable :: dprime

    real(pb_dp) :: errrms
    real(pb_dp) :: tol

    logical :: soln_pass
    
    bprime = b
    dprime = d
    
    call tdma(a, bprime, c, dprime)

    tol = epsilon(errrms)
    errrms = sum((x - dprime)**2)
    errrms = sqrt(errrms / n)
    soln_pass = errrms <= (tol * sqrt(sum(x**2) / n))
    if ((expect .and. soln_pass) .or. &
         ((.not. expect) .and. (.not. soln_pass))) then
       print *, "TDMA passed: ", errrms
    else
       print *, "TDMA failed: ", errrms
       passing = .false.
    end if
    
  end subroutine check_tdma
  
end program test_tdma
