module compact_schemes

  use constants
  use tridsol, only: tdma_periodic
  
  implicit none

  private
  public :: grad_1d
  public :: interp_1d
contains

  ! Compute the 1D gradient (staggered) of a field
  subroutine grad_1d(f, dx, df)

    real(pb_dp), dimension(:), intent(in) :: f   ! Field
    real(pb_dp), intent(in) :: dx                ! Grid spacing
    real(pb_dp), dimension(:), intent(out) :: df ! Gradient

    real(pb_dp) :: a, b, alpha ! Scheme parameters

    real(pb_dp), dimension(:), allocatable :: ld, d, ud ! Tridiagonal coefficients

    integer :: n

    n = size(f)
    if (size(df) /= n) then
       print *, "ERROR: periodic gradient is same length as field!"
       stop 7
    end if
    
    !! Assemble system
    allocate(ld(n))
    allocate(d(n))
    allocate(ud(n))
    
    ! Set coefficients
    a = 63.0_pb_dp / 62.0_pb_dp / dx
    b = 17.0_pb_dp / 62.0_pb_dp / (3.0_pb_dp * dx)
    alpha = 9.0_pb_dp / 62.0_pb_dp
    ld(:) = alpha
    d(:) = 1.0_pb_dp
    ud(:) = alpha
    call grad_1d_rhs(a, b, -1, f, df)
    
    !! Solve system
    call tdma_periodic(ld, d, ud, df)
    
    !! Cleanup
    deallocate(ld)
    deallocate(d)
    deallocate(ud)
    
  end subroutine grad_1d

  subroutine interp_1d(f, fi)

    real(pb_dp), dimension(:), intent(in) :: f   ! Field
    real(pb_dp), dimension(:), intent(out) :: fi ! Interpolated field

    real(pb_dp) :: a, b, alpha ! Scheme parameters

    real(pb_dp), dimension(:), allocatable :: ld, d, ud ! Tridiagonal coefficients

    integer :: n

    n = size(f)
    if (size(fi) /= n) then
       print *, "ERROR: periodic gradient is same length as field!"
       stop 7
    end if
    
    !! Assemble system
    allocate(ld(n))
    allocate(d(n))
    allocate(ud(n))
    
    ! Set coefficients
    a = 0.75_pb_dp
    b = 1.0_pb_dp / 20.0_pb_dp
    alpha = 3.0_pb_dp / 10.0_pb_dp
    ld(:) = alpha
    d(:) = 1.0_pb_dp
    ud(:) = alpha
    call grad_1d_rhs(a, b, +1, f, fi)
    
    !! Solve system
    call tdma_periodic(ld, d, ud, fi)
    
    !! Cleanup
    deallocate(ld)
    deallocate(d)
    deallocate(ud)
    
  end subroutine interp_1d
  
  pure subroutine grad_1d_rhs(a, b, opsign, f, rhs)

    real(pb_dp), intent(in) :: a, b ! Scheme parameters
    integer, intent(in) :: opsign   ! Set the sign of the finite difference scheme
                                    ! (-1 difference, +1 interpolation)
    real(pb_dp), dimension(:), intent(in) :: f    ! Field
    real(pb_dp), dimension(:), intent(out) :: rhs ! RHS

    integer :: i
    integer :: n

    n = size(f)

    rhs(1) = a * (f(2) + opsign * f(1)) + b * (f(3) + opsign * f(n))
    do i = 2, n - 2
       rhs(i) = a * (f(i + 1) + opsign * f(i)) + b * (f(i + 2) + opsign * f(i - 1))
    end do
    rhs(n - 1) = a * (f(n) + opsign * f(n - 1)) + b * (f(1) + opsign * f(n - 2))
    rhs(n) = a * (f(1) + opsign * f(n)) + b * (f(2) + opsign * f(n - 1))
    
  end subroutine grad_1d_rhs
  
end module compact_schemes
