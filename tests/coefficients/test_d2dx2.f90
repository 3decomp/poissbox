!!! tests/coefficients/test_d2dx2.f90
!
! Test computing the coefficients for the 1-D Laplacian operator d2/dx2.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program test_d2dx2

  use constants
  
  implicit none

  real(pb_dp), parameter :: a = 2.718 ! Parameter for the quadratic formula
  real(pb_dp), parameter :: b = 1.414 ! Parameter for the quadratic formula
  real(pb_dp), parameter :: c = 1.848 ! Parameter for the quadratic formula
  
  real(pb_dp), parameter :: x = 1.618  ! Location for field evaluation
  real(pb_dp), parameter :: dx = 0.155 ! Grid spacing
  
  real(pb_dp), dimension(3), parameter :: fc = [ c, c, c ]                  ! Constant field
  real(pb_dp), dimension(3), parameter :: fg = b * [ x - dx, x, x + dx ]    ! Constant gradient field
  real(pb_dp), dimension(3), parameter :: fq = a * [ x - dx, x, x + dx ]**2 ! Quadratic field

  real(pb_dp), parameter :: shift = 17.29 ! Shift to apply to fields

  logical :: test_pass = .true.
  
  call test_constant_field()
  call test_constant_grad()
  call test_quadratic_field()
  
  ! Laplacian of ax^2 + bx + c -> 2ax
  ! Proportional to scaling of the quadratic term
  ! Invariant to shifting

  if (.not. test_pass) then
     stop 1
  end if
  
contains

  subroutine test_constant_field()

    real(pb_dp):: expected_lapl

    !! Laplacian of a constant field is 0
    expected_lapl = 0.0

    call test_lapl(fc, dx, expected_lapl, "constant")
    call test_scaled_lapl(fc, dx, expected_lapl, "constant")
    call test_shifted_lapl(fc, dx, expected_lapl, "constant")
    call test_spacing_lapl(fc, dx, expected_lapl, "constant")
    
  end subroutine test_constant_field

  subroutine test_constant_grad()

    real(pb_dp):: expected_lapl

    !! Laplacian of a constant gradient field is 0
    expected_lapl = 0.0

    call test_lapl(fg, dx, expected_lapl, "constant gradient")
    call test_scaled_lapl(fg, dx, expected_lapl, "constant gradient")
    call test_shifted_lapl(fg, dx, expected_lapl, "constant gradient")
    call test_spacing_lapl(fg, dx, expected_lapl, "constant gradient")
    
  end subroutine test_constant_grad

  subroutine test_quadratic_field()

    real(pb_dp):: expected_lapl

    !! Laplacian of a quadratic field is 2 * a
    expected_lapl = 2 * a

    call test_lapl(fq, dx, expected_lapl, "quadratic")
    call test_scaled_lapl(fq, dx, expected_lapl, "quadratic")
    call test_shifted_lapl(fq, dx, expected_lapl, "quadratic")
    
  end subroutine test_quadratic_field

  subroutine test_lapl(f, dx, expected_lapl, name)

    real(pb_dp), dimension(3), intent(in) :: f
    real(pb_dp), intent(in) :: dx
    real(pb_dp), intent(in) :: expected_lapl
    character(len=*), intent(in) :: name

    real(pb_dp):: lapl
    
    lapl = evaluate_laplacian(f, dx)
    if (.not. feq(lapl * dx**2, expected_lapl * dx**2)) then
       print *, "Laplacian of " // name // " field: expected ", expected_lapl, " got ", lapl
       test_pass = .false.
    end if
    
  end subroutine test_lapl

  subroutine test_scaled_lapl(f, dx, expected_lapl, name)

    real(pb_dp), dimension(3), intent(in) :: f
    real(pb_dp), intent(in) :: dx
    real(pb_dp), intent(in) :: expected_lapl
    character(len=*), intent(in) :: name

    real(pb_dp):: lapl
    
    lapl = evaluate_laplacian(2 * f, dx)
    if (.not. feq(lapl, 2 * expected_lapl)) then
       print *, "Laplacian of " // name // " field x 2: expected ", 2 * expected_lapl, &
                " got ", lapl
       test_pass = .false.
    end if
    lapl = evaluate_laplacian(f / 2, dx)
    if (.not. feq(lapl * dx**2, expected_lapl * dx**2 / 2)) then
       print *, "Laplacian of " // name // " field / 2: expected ", expected_lapl / 2, &
                " got ", lapl
       test_pass = .false.
    end if
    
  end subroutine test_scaled_lapl

  subroutine test_shifted_lapl(f, dx, expected_lapl, name)

    real(pb_dp), dimension(3), intent(in) :: f
    real(pb_dp), intent(in) :: dx
    real(pb_dp), intent(in) :: expected_lapl
    character(len=*), intent(in) :: name

    real(pb_dp):: lapl

    lapl = evaluate_laplacian(f + shift, dx)
    if (.not. feq(lapl * dx**2, expected_lapl * dx**2)) then
       print *, "Laplacian of " // name // " field shifted up: expected ", expected_lapl, &
                " got ", lapl
       test_pass = .false.
    end if
    lapl = evaluate_laplacian(f - shift, dx)
    if (.not. feq(lapl * dx**2, expected_lapl * dx**2)) then
       print *, "Laplacian of " // name // " field shifted down: expected ", expected_lapl, &
                " got ", lapl
       test_pass = .false.
    end if
    
  end subroutine test_shifted_lapl

  subroutine test_spacing_lapl(f, dx, expected_lapl, name)

    real(pb_dp), dimension(3), intent(in) :: f
    real(pb_dp), intent(in) :: dx
    real(pb_dp), intent(in) :: expected_lapl
    character(len=*), intent(in) :: name

    real(pb_dp):: lapl

    lapl = evaluate_laplacian(f, 2 * dx)
    if (.not. feq(lapl * (2 * dx)**2, expected_lapl * (2 * dx)**2)) then
       print *, "Laplacian of " // name // " field on 2 * dx grid: expected ", expected_lapl, &
                " got ", lapl
       test_pass = .false.
    end if
    lapl = evaluate_laplacian(f, dx / 2)
    if (.not. feq(lapl * (dx / 2)**2, expected_lapl * (dx / 2)**2)) then
       print *, "Laplacian of " // name // " field on dx / 2 grid: expected ", expected_lapl, &
                " got ", lapl
       test_pass = .false.
    end if
    
  end subroutine test_spacing_lapl
  
  real(pb_dp)function evaluate_laplacian(f, dx)

    use coefficients, only : lapl_1d_coeffs

    real(pb_dp), dimension(3), intent(in) :: f
    real(pb_dp), intent(in) :: dx

    real(pb_dp), dimension(3) :: coeffs

    coeffs = lapl_1d_coeffs(dx)

    ! Evaluating by grouping the first and last terms seems to be a bit more accurate than a simple
    ! dot product, i.e. approximate as
    !   (f_i+1 + f_i-1) - 2 * f_i
    ! vs
    !   f_i+1 - 2 * f_i + f_i-1
    evaluate_laplacian = (coeffs(1) * f(1) + coeffs(3) * f(3)) + coeffs(2) * f(2)
    
  end function evaluate_laplacian

  logical function feq(val, ref, opt_tol)

    real(pb_dp), intent(in) :: val ! The value under test
    real(pb_dp), intent(in) :: ref ! The reference value
    real(pb_dp), intent(in), optional :: opt_tol ! The tolerance of equality

    real(pb_dp):: tol
    real(pb_dp):: delta

    if (present(opt_tol)) then
       tol = opt_tol
    else
       tol = 100 * epsilon(ref)
    end if

    delta = abs(val - ref)

    feq = (delta <= tol * abs(ref)) .or. (delta <= tol)

  end function feq

end program test_d2dx2
