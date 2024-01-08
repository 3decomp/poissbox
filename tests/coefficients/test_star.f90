!!! tests/coefficients/test_star.f90
!
! Test the star coefficients calculation.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

program test_star

  use constants
  
  implicit none

  real(pb_dp), parameter :: a = 2.718 ! Parameter for the quadratic formula
  real(pb_dp), parameter :: b = 1.414 ! Parameter for the quadratic formula
  real(pb_dp), parameter :: c = 1.848 ! Parameter for the quadratic formula
  
  real(pb_dp), parameter :: x = 1.618  ! Location for field evaluation
  real(pb_dp), parameter :: dx = 0.155 ! Grid spacing

  real(pb_dp), dimension(3), parameter :: fgx = b * [ x - dx, x, x + dx ]    ! Constant gradient field
  real(pb_dp), dimension(3), parameter :: fqx = a * [ x - dx, x, x + dx ]**2 ! Quadratic field
  
  real(pb_dp), dimension(3, 3, 3) :: fc
  real(pb_dp), dimension(3, 3, 3) :: fg
  real(pb_dp), dimension(3, 3, 3) :: fq

  logical :: test_pass
  
  call init_test()

  call test_constant_field()
  call test_constant_grad()
  call test_quadratic_field()
  
  if (.not. test_pass) then
     stop 1
  end if
  
contains

  subroutine init_test

    integer :: i, j, k

    test_pass = .true.
    
    ! Constant field
    fc(:, :, :) = c

    ! Constant gradient field
    do k = 1, 3
       do j = 1, 3
          fg(:, j, k) = fgx(:)
       end do
    end do
    do k = 1, 3
       do i = 1, 3
          fg(i, :, k) = fg(i, :, k) + fgx(:)
       end do
    end do
    do j = 1, 3
       do i = 1, 3
          fg(i, j, :) = fg(i, j, :) + fgx(:)
       end do
    end do
    
    ! Quadratic field
    do k = 1, 3
       do j = 1, 3
          fq(:, j, k) = fqx(:)
       end do
    end do
    do k = 1, 3
       do i = 1, 3
          fq(i, :, k) = fq(i, :, k) + fqx(:)
       end do
    end do
    do j = 1, 3
       do i = 1, 3
          fq(i, j, :) = fq(i, j, :) + fqx(:)
       end do
    end do

  end subroutine init_test

  subroutine test_constant_field()

    real(pb_dp):: expected_lapl

    expected_lapl = 0.0

    call test_lapl(fc, dx, expected_lapl, "constant")

  end subroutine test_constant_field

  subroutine test_constant_grad()

    real(pb_dp):: expected_lapl

    expected_lapl = 0.0

    call test_lapl(fg, dx, expected_lapl, "constant gradient")

  end subroutine test_constant_grad

  subroutine test_quadratic_field()

    real(pb_dp):: expected_lapl

    expected_lapl = 3 * (2 * a)

    call test_lapl(fq, dx, expected_lapl, "quadratic")

  end subroutine test_quadratic_field
  
  subroutine test_lapl(f, dx, expected_lapl, name)

    real(pb_dp), dimension(3, 3, 3), intent(in) :: f
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
  
  real(pb_dp)function evaluate_laplacian(f, dx)

    use coefficients, only : lapl_star_coeffs

    real(pb_dp), dimension(3, 3, 3), intent(in) :: f
    real(pb_dp), intent(in) :: dx

    real(pb_dp), dimension(3, 3, 3) :: coeffs

    coeffs = lapl_star_coeffs(dx, dx, dx)

    evaluate_laplacian = dot_product(reshape(f, [size(f)]), &
                                     reshape(coeffs, [size(coeffs)]))
    
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
       tol = 100 * (1.1 * epsilon(ref))
    end if

    delta = abs(val - ref)

    print *, tol, tol * abs(ref), delta, abs(ref)
    feq = (delta <= tol * abs(ref)) .or. (delta <= tol)

  end function feq
  
end program test_star

