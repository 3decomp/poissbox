!!! tests/coefficients/test_compact.f90
!
! Test the coefficients for compact schemes. For a polynomials below the order of the scheme the
! discrete equation should be satisfied exactly.
!
! These tests are for staggered schemes of the form
!
!  alpha f'_{i - 1/2} + f'_{i + 1/2} + alpha f'_{i + 3/2} = a (f_{i + 1} - f_{i - 1}) / 2dx ...
!
! with a similar equation for interpolations.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program test_compact

  use constants
  
  implicit none

  real(pb_dp), parameter :: L = 6.28
  real(pb_dp), parameter :: n = 128
  real(pb_dp), parameter :: dx = L / n

  real(pb_dp), parameter :: a = 3.14
  real(pb_dp), parameter :: b = 0.817
  real(pb_dp), parameter :: c = -7.362
  real(pb_dp), parameter :: d = 8.981
  real(pb_dp), parameter :: e = -10.22
  real(pb_dp), parameter :: g = 0.071
  
  real(pb_dp), dimension(4) :: f0, f1, f2, f3, f4, f5       ! Test functions
  real(pb_dp), dimension(3) :: df0, df1, df2, df3, df4, df5 ! Test derivatives
  real(pb_dp), dimension(3) :: fi0, fi1, fi2, fi3, fi4, fi5 ! Test interpolations

  logical :: passing
  
  call init()
  call check_derivatives()
  call check_interpolation()
  
  if (.not. passing) then
     print *, "FAIL"
     stop 1
  end if
  
contains

  subroutine init()
    !! Initialise the test functions and set their derivatives and interpolants.

    call init_fn(0, a, f0, df0, fi0)
    call init_fn(1, b, f1, df1, fi1)
    call init_fn(2, c, f2, df2, fi2)
    call init_fn(3, d, f3, df3, fi3)
    call init_fn(4, e, f4, df4, fi4)
    call init_fn(5, g, f5, df5, fi5)

    ! Combine polynomial orders for more complex functions
    f1 = f1 + f0
    f2 = f2 + f1
    f3 = f3 + f2
    f4 = f4 + f3
    f5 = f5 + f4

    df1 = df1 + df0
    df2 = df2 + df1
    df3 = df3 + df2
    df4 = df4 + df3
    df5 = df5 + df4

    fi1 = fi1 + fi0
    fi2 = fi2 + fi1
    fi3 = fi3 + fi2
    fi4 = fi4 + fi3
    fi5 = fi5 + fi4

    passing = .true.
    
  end subroutine init

  subroutine init_fn(p, m, f, df, fi)
    !! Initialise a function of a single order, with its associated derivatives and interpolants.
    
    integer, intent(in) :: p ! Polynomial order
    real(pb_dp), intent(in) :: m    ! Scale factor
    real(pb_dp), dimension(:), intent(out) :: f  ! Function
    real(pb_dp), dimension(3), intent(out) :: df ! Derivatives
    real(pb_dp), dimension(3), intent(out) :: fi ! Interpolants

    real(pb_dp) :: x
    integer :: i

    x = 0.0
    do i = 1, size(f)
       f(i) = m * x**p
       x = x + dx
    end do

    ! Evaluate derivative around centre-point
    x = (1.5 * dx)
    df(1) = (p * m) * (x - dx)**(p - 1)
    df(2) = (p * m) * x**(p - 1)
    df(3) = (p * m) * (x + dx)**(p - 1)

    ! Evaluate interpolant around centre-point
    fi(1) = m * (x - dx)**p
    fi(2) = m * x**p
    fi(3) = m * (x + dx)**p
    
  end subroutine init_fn

  subroutine check_derivatives()

    real(pb_dp), dimension(3) :: alpha = [9.0_pb_dp / 62.0_pb_dp, 1.0_pb_dp, 9.0_pb_dp / 62.0_pb_dp]
    real(pb_dp), dimension(4) :: rhs = [-(17.0_pb_dp / 62.0_pb_dp) / (3.0_pb_dp * dx), -(63.0_pb_dp / 62.0_pb_dp) / dx, &
                                         (63.0_pb_dp / 62.0_pb_dp) / dx, (17.0_pb_dp / 62.0_pb_dp) / (3.0_pb_dp * dx)]
    
    call check_scheme(alpha, df0, rhs, f0, "dfdx", "0")
    call check_scheme(alpha, df1, rhs, f1, "dfdx", "1")
    call check_scheme(alpha, df2, rhs, f2, "dfdx", "2")
    call check_scheme(alpha, df3, rhs, f3, "dfdx", "3")
    call check_scheme(alpha, df4, rhs, f4, "dfdx", "4")
    call check_scheme(alpha, df5, rhs, f5, "dfdx", "5")
    
  end subroutine check_derivatives

  subroutine check_interpolation()

    real(pb_dp), dimension(3) :: alpha = [3.0_pb_dp / 10.0_pb_dp, 1.0_pb_dp, 3.0_pb_dp / 10.0_pb_dp]
    real(pb_dp), dimension(4) :: rhs = [(1.0_pb_dp / 20.0_pb_dp), 0.75_pb_dp, &
                                        0.75_pb_dp, (1.0_pb_dp / 20.0_pb_dp)]
    
    call check_scheme(alpha, fi0, rhs, f0, "interp", "0")
    call check_scheme(alpha, fi1, rhs, f1, "interp", "1")
    call check_scheme(alpha, fi2, rhs, f2, "interp", "2")
    call check_scheme(alpha, fi3, rhs, f3, "interp", "3")
    call check_scheme(alpha, fi4, rhs, f4, "interp", "4")
    call check_scheme(alpha, fi5, rhs, f5, "interp", "5")
    
  end subroutine check_interpolation

  subroutine check_scheme(alpha, df, rhs, f, op, pwr)

    real(pb_dp), dimension(3), intent(in) :: alpha
    real(pb_dp), dimension(3), intent(in) :: df
    real(pb_dp), dimension(4), intent(in) :: rhs
    real(pb_dp), dimension(4), intent(in) :: f
    character(len=*), intent(in) :: op
    character(len=*), intent(in) :: pwr

    real(pb_dp), parameter :: tol = 100 * epsilon(alpha)

    real(pb_dp) :: delta

    delta = dot_product(rhs, f) - dot_product(alpha, df)
    if (abs(delta) > tol) then
       passing = .false.
       print *, "Scheme " // op // " failed for p = " // pwr // ", delta = ", delta, "; tol = ", tol
    end if
    
  end subroutine check_scheme
  
end program test_compact
