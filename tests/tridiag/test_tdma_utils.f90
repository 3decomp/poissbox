module test_tdma_utils

  use constants
  
  implicit none

  private
  public :: init
  
contains
  
  subroutine init(n, a, b, c, x, d, opt_periodic)

    integer, intent(in) :: n ! Problem size
    real(pb_dp), dimension(:), allocatable, intent(out) :: a ! Sub-diagonal
    real(pb_dp), dimension(:), allocatable, intent(out) :: b ! Super-diagonal
    real(pb_dp), dimension(:), allocatable, intent(out) :: c ! Diagonal
    real(pb_dp), dimension(:), allocatable, intent(out) :: x ! Solution
    real(pb_dp), dimension(:), allocatable, intent(out) :: d ! RHS
    logical, intent(in), optional :: opt_periodic
    
    integer :: i
    logical :: periodic

    if (present(opt_periodic)) then
       periodic = opt_periodic
    else
       periodic = .false.
    end if
    
    allocate(a(n), b(n), c(n), x(n), d(n))

    ! Randomise system matrix + solution
    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(x)

    if (.not. periodic) then
       a(1) = 0.0_pb_dp
       c(n) = 0.0_pb_dp
    end if
    
    ! Ensure diagonal dominance
    do i = 1, n
       do while(b(i) == 0.0_pb_dp)
          call random_number(b(i))
       end do
       do while (abs(b(i)) < (abs(a(i)) + abs(c(i))))
          b(i) = 10 * b(i)
       end do
    end do

    ! Construct RHS
    d(1) = b(1) * x(1) + c(1) * x(2)
    if (periodic) then
       d(1) = a(1) * x(n) + d(1)
    end if
    do i = 2, n - 1
       d(i) = a(i) * x(i - 1) + b(i) * x(i) + c(i) * x(i + 1)
    end do
    d(n) = a(n) * x(n - 1) + b(n) * x(n)
    if (periodic) then
       d(n) = c(n) * x(1) + d(n)
    end if
    
  end subroutine init
  
end module test_tdma_utils
