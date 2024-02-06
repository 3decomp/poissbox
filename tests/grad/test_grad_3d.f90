program test_grad_3d

  use constants
  use compact_schemes
  
  implicit none

  real(pb_dp), parameter :: pi = 4 * atan(1.0_pb_dp)

  integer, parameter :: nx = 64    ! Problem size (nodes)
  integer, parameter :: ny = 64    ! Problem size (nodes)
  integer, parameter :: nz = 64    ! Problem size (nodes)
  real(pb_dp), parameter :: L = pi ! Domain size
  real(pb_dp), parameter :: dx = L / nx ! Grid spacing
  real(pb_dp), parameter :: dy = L / ny ! Grid spacing
  real(pb_dp), parameter :: dz = L / nz ! Grid spacing
  
  real(pb_dp), dimension(:, :, :), allocatable :: f  ! The function
  real(pb_dp), dimension(:, :, :, :), allocatable :: df ! The gradient

  logical :: passing

  passing = .true.
  
  call init()

  call check_constant_field()
  call check_varying_field_x()
  call check_varying_field_y()
  call check_varying_field_z()
  call check_varying_field()
  
  call fin()

  if (.not. passing) then
     print *, "FAIL"
     stop 1
  end if
  
contains

  subroutine init()

    allocate(f(nx, ny, nz))
    allocate(df(nx, ny, nz, 3))
    
  end subroutine init

  subroutine fin()

    if (allocated(f)) then
       deallocate(f)
    end if
    if (allocated(df)) then
       deallocate(df)
    end if
    
  end subroutine fin

  subroutine check_constant_field()

    real(pb_dp) :: rms
    
    f(:, :, :) = 2.8170923_pb_dp ! Arbitrary constant field
    df(:, :, :, :) = 73.29_pb_dp ! Non-zero (i.e. wrong) gradient
    
    call grad(f, [dx, dy, dz], df)

    rms = sqrt(sum(df**2) / nx / ny / nz / 3)
    if (rms > (100 * epsilon(rms))) then
       print *, "FAIL: RMS dfdx = ", rms, "(f = const)"
       passing = .false.
    else
       print *, "PASS: dfdx (f = const)"
    end if

    call interp(f, df(:, :, :, 1))

    rms = sqrt(sum((df(:, :, :, 1) - f)**2) / nx / ny / nz)
    if (rms > (100 * epsilon(rms))) then
       print *, "FAIL: RMS [f] = ", rms, "(f = const)"
       passing = .false.
    else
       print *, "PASS: [f] (f = const)"
    end if
    
  end subroutine check_constant_field

  subroutine check_varying_field()

    real(pb_dp) :: x, y, z

    integer :: i, j, k

    real(pb_dp) :: expect_x, expect_y, expect_z
    real(pb_dp) :: rms, rms_x, rms_y, rms_z

    logical :: test_pass
    
    z = 0.0_pb_dp
    do k = 1, nz
       y = 0.0_pb_dp
       do j = 1, ny
          x = 0.0_pb_dp
          do i = 1, nx
             f(i, j, k) = sin(x) + sin(y) + sin(z)
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do

    call grad(f, [dx, dy, dz], df)

    rms = 0.0_pb_dp
    rms_x = 0.0_pb_dp; rms_y = 0.0_pb_dp; rms_z = 0.0_pb_dp
    z = 0.5_pb_dp * dz
    do k = 1, nz
       y = 0.5_pb_dp * dy
       expect_z = cos(z)
       do j = 1, ny
          x = 0.5_pb_dp * dx
          expect_y = cos(y)
          do i = 1, nx
             expect_x = cos(x)

             rms_x = rms_x + (df(i, j, k, 1) - expect_x)**2
             rms_y = rms_y + (df(i, j, k, 2) - expect_y)**2
             rms_z = rms_z + (df(i, j, k, 3) - expect_z)**2
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    ! print *, rms_x, rms_y, rms_z
    rms = rms_x + rms_y + rms_z
    rms = sqrt(rms / nx / ny / nz / 3)

    test_pass = (rms <= 1.0e-11)
    test_pass = (.not. (rms /= rms)) .and. test_pass
    if (.not. test_pass) then
       print *, "FAIL: RMS dfdx = ", rms, "variable f"
       passing = .false.
    else
       print *, "PASS: RMS dfdx = ", rms, "variable f"
    end if
    passing = .false.
    
    call interp(f, df(:, :, :, 1))

    rms = 0.0_pb_dp
    do k = 1, nz
       y = 0.5_pb_dp * dy
       do j = 1, ny
          x = 0.5_pb_dp * dx
          do i = 1, nx
             x = x + dx

             rms = rms + (df(i, j, k, 1) - (sin(x) + sin(y) + sin(z)))**2
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    rms = sqrt(rms / nx / ny / nz)
    print *, rms

  end subroutine check_varying_field

  subroutine check_varying_field_x()

    real(pb_dp) :: x, y, z

    integer :: i, j, k

    real(pb_dp) :: expect_x, expect_y, expect_z
    real(pb_dp) :: rms, rms_x, rms_y, rms_z

    logical :: test_pass
    
    z = 0.0_pb_dp
    do k = 1, nz
       y = 0.0_pb_dp
       do j = 1, ny
          x = 0.0_pb_dp
          do i = 1, nx
             f(i, j, k) = sin(x)
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do

    call grad(f, [dx, dy, dz], df)

    rms = 0.0_pb_dp
    rms_x = 0.0_pb_dp; rms_y = 0.0_pb_dp; rms_z = 0.0_pb_dp
    z = 0.5_pb_dp * dz
    do k = 1, nz
       y = 0.5_pb_dp * dy
       expect_z = 0.0_pb_dp
       do j = 1, ny
          x = 0.5_pb_dp * dx
          expect_y = 0.0_pb_dp
          do i = 1, nx
             expect_x = cos(x)

             rms_x = rms_x + (df(i, j, k, 1) - expect_x)**2
             rms_y = rms_y + (df(i, j, k, 2) - expect_y)**2
             rms_z = rms_z + (df(i, j, k, 3) - expect_z)**2
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    rms = rms_x + rms_y + rms_z
    rms = sqrt(rms / nx / ny / nz / 3)
    rms_x = sqrt(rms_x / nx)
    rms_y = sqrt(rms_y / ny)
    rms_z = sqrt(rms_z / nz)

    test_pass = (rms <= 1.0e-11)
    test_pass = (.not. (rms /= rms)) .and. test_pass
    if (.not. test_pass) then
       print *, "FAIL: RMS dfdx = ", rms, "f=f(x)"
       print *, rms_x, rms_y, rms_z
       passing = .false.
    else
       print *, "PASS: RMS dfdx = ", rms, "f=f(x)"
    end if
    passing = .false.
    
    call interp(f, df(:, :, :, 1))

    rms = 0.0_pb_dp
    do k = 1, nz
       y = 0.5_pb_dp * dy
       do j = 1, ny
          x = 0.5_pb_dp * dx
          do i = 1, nx
             x = x + dx

             rms = rms + (df(i, j, k, 1) - sin(x))**2
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    rms = sqrt(rms / nx / ny / nz)
    print *, rms
    
  end subroutine check_varying_field_x

  subroutine check_varying_field_y()

    real(pb_dp) :: x, y, z

    integer :: i, j, k

    real(pb_dp) :: expect_x, expect_y, expect_z
    real(pb_dp) :: rms, rms_x, rms_y, rms_z

    logical :: test_pass
    
    z = 0.0_pb_dp
    do k = 1, nz
       y = 0.0_pb_dp
       do j = 1, ny
          x = 0.0_pb_dp
          do i = 1, nx
             f(i, j, k) = sin(y)
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do

    call grad(f, [dx, dy, dz], df)

    rms = 0.0_pb_dp
    rms_x = 0.0_pb_dp; rms_y = 0.0_pb_dp; rms_z = 0.0_pb_dp
    z = 0.5_pb_dp * dz
    do k = 1, nz
       y = 0.5_pb_dp * dy
       expect_z = 0.0_pb_dp
       do j = 1, ny
          x = 0.5_pb_dp * dx
          expect_y = cos(y)
          do i = 1, nx
             expect_x = 0.0_pb_dp

             rms_x = rms_x + (df(i, j, k, 1) - expect_x)**2
             rms_y = rms_y + (df(i, j, k, 2) - expect_y)**2
             rms_z = rms_z + (df(i, j, k, 3) - expect_z)**2
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    ! print *, rms_x, rms_y, rms_z
    rms = rms_x + rms_y + rms_z
    rms = sqrt(rms / nx / ny / nz / 3)

    test_pass = (rms <= 1.0e-11)
    test_pass = (.not. (rms /= rms)) .and. test_pass
    if (.not. test_pass) then
       print *, "FAIL: RMS dfdx = ", rms, "f=f(y)"
       passing = .false.
    else
       print *, "PASS: RMS dfdx = ", rms, "f=f(y)"
    end if
    passing = .false.
    
  end subroutine check_varying_field_y

  subroutine check_varying_field_z()

    real(pb_dp) :: x, y, z

    integer :: i, j, k

    real(pb_dp) :: expect_x, expect_y, expect_z
    real(pb_dp) :: rms, rms_x, rms_y, rms_z

    logical :: test_pass
    
    z = 0.0_pb_dp
    do k = 1, nz
       y = 0.0_pb_dp
       do j = 1, ny
          x = 0.0_pb_dp
          do i = 1, nx
             f(i, j, k) = sin(z)
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do

    call grad(f, [dx, dy, dz], df)

    rms = 0.0_pb_dp
    rms_x = 0.0_pb_dp; rms_y = 0.0_pb_dp; rms_z = 0.0_pb_dp
    z = 0.5_pb_dp * dz
    do k = 1, nz
       y = 0.5_pb_dp * dy
       expect_z = cos(z)
       do j = 1, ny
          x = 0.5_pb_dp * dx
          expect_y = 0.0_pb_dp
          do i = 1, nx
             expect_x = 0.0_pb_dp

             rms_x = rms_x + (df(i, j, k, 1) - expect_x)**2
             rms_y = rms_y + (df(i, j, k, 2) - expect_y)**2
             rms_z = rms_z + (df(i, j, k, 3) - expect_z)**2
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    ! print *, rms_x, rms_y, rms_z
    rms = rms_x + rms_y + rms_z
    rms = sqrt(rms / nx / ny / nz / 3)

    test_pass = (rms <= 1.0e-11)
    test_pass = (.not. (rms /= rms)) .and. test_pass
    if (.not. test_pass) then
       print *, "FAIL: RMS dfdx = ", rms, "f=f(z)"
       passing = .false.
    else
       print *, "PASS: RMS dfdx = ", rms, "f=f(z)"
    end if
    passing = .false.
    
  end subroutine check_varying_field_z
  
end program test_grad_3d
