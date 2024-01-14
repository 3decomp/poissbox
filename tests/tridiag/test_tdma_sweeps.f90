program test_tdma_sweeps

  use constants

  use tridsol

  use test_tdma_utils
  
  implicit none

  integer, parameter :: n = 128 ! Problem size
  real(pb_dp), dimension(:), allocatable :: a, b, c, x, d ! Tridiagonal parameters

  logical :: passing

  passing = .true.
  
  call init(n, a, b, c, x, d)

  call check_fwd_sweep(a, b, c, x, d)
  call check_bwd_sweep(b, c, x)
  
  call fin()

contains

  subroutine fin()

    deallocate(a, b, c, x, d)

    if (.not. passing) then
       stop 1
    end if
    
  end subroutine fin

  subroutine check_fwd_sweep(a, b, c, x, d)

    real(pb_dp), dimension(:), allocatable, intent(in) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: b ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: x ! Solution
    real(pb_dp), dimension(:), allocatable, intent(in) :: d ! RHS

    real(pb_dp), dimension(:), allocatable :: bprime, dprime

    bprime = b
    dprime = d
    
    ! Forward sweep should convert into upper-triangular form, i.e.
    !   b'(i) + x(i) + c'(i) * x(i + 1) = d'(i)
    ! where primes indicate (potentially) modified values arising from forward sweep.
    call fwd_sweep(a, bprime, c, dprime)

    call check_utri(bprime, c, x, dprime, "fwd")
    
  end subroutine check_fwd_sweep

  subroutine check_bwd_sweep(b, c, x)

    real(pb_dp), dimension(:), allocatable, intent(in) :: b ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: x ! Solution

    real(pb_dp), dimension(:), allocatable :: d, d2 ! Temporary RHS/solution

    integer :: i

    ! The backward sweep solves an upper-triangular system, this test constructs a temporary
    ! upper-triangular system and attempts to solve it.
    allocate(d, mold=x)
    do i = 1, n - 1
       d(i) = b(i) * x(i) + c(i) * x(i + 1)
    end do
    d(n) = b(n) * x(n)
    d2 = d

    call bwd_sweep(b, c, d)

    call check_utri(b, c, d, d2, "bwd")
    
    deallocate(d)
    deallocate(d2)
    
  end subroutine check_bwd_sweep

  subroutine check_utri(b, c, x, d, str)

    real(pb_dp), dimension(:), allocatable, intent(in) :: b ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: c ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(in) :: x ! Solution
    real(pb_dp), dimension(:), allocatable, intent(in) :: d ! RHS
    character(len=*), intent(in) :: str

    real(pb_dp) :: delta
    real(pb_dp) :: r ! Residual
    real(pb_dp) :: t ! Tolerance
    
    integer :: i

    t = epsilon(d) * sqrt(sum(d**2) / n)

    r = 0.0_pb_dp
    do i = 1, n - 1
       delta = d(i) - (b(i) * x(i) + c(i) * x(i + 1))
       r = r + delta**2
    end do
    delta = d(n) - b(n) * x(n)
    r = r + delta**2

    r = sqrt(r / n)
    if (r > t) then
       print *, str//"-sweep not upper-triangular: ", r, " to tolerance ", t
       passing = .false.
    end if

  end subroutine check_utri

end program test_tdma_sweeps
