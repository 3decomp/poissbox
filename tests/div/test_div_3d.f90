program test_div_3d

  use constants
  use compact_schemes
  
  implicit none

  real(pb_dp), parameter :: pi = 4 * atan(1.0_pb_dp)

  integer, parameter :: nx = 64    ! Problem size (nodes)
  integer, parameter :: ny = 64    ! Problem size (nodes)
  integer, parameter :: nz = 64    ! Problem size (nodes)
  real(pb_dp), parameter :: L = 2 * pi  ! Domain size
  real(pb_dp), parameter :: dx = L / nx ! Grid spacing
  real(pb_dp), parameter :: dy = L / ny ! Grid spacing
  real(pb_dp), parameter :: dz = L / nz ! Grid spacing
  
  real(pb_dp), dimension(:, :, :, :), allocatable :: f ! The (vector) function
  real(pb_dp), dimension(:, :, :), allocatable :: df   ! The divergence

  logical :: passing

  passing = .true.
  
  call init()

  call check_constant_field()
  call check_varying_field()
  
  call fin()

  if (.not. passing) then
     print *, "FAIL"
     stop 1
  end if
  
contains

  subroutine init()

    allocate(f(nx, ny, nz, 3))
    allocate(df(nx, ny, nz))
    
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
    
    f(:, :, :, :) = 2.8170923_pb_dp ! Arbitrary constant field
    df(:, :, :) = 73.29_pb_dp       ! Non-zero (i.e. wrong) divergence
    
    call div(f, [dx, dy, dz], df)

    rms = sqrt(sum(df**2) / nx / ny / nz)
    if (rms > (100 * epsilon(rms))) then
       print *, "FAIL: RMS div(f) = ", rms, "(f = const)"
       passing = .false.
    else
       print *, "PASS: div(f) (f = const)"
    end if

    call interp_div(f(:, :, :, 1), df)

    rms = sqrt(sum((df - f(:, :, :, 1))**2) / nx / ny / nz)
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

    real(pb_dp) :: expect
    real(pb_dp) :: rms

    logical :: test_pass
    
    z = 0.0_pb_dp * dz
    do k = 1, nz
       y = 0.0_pb_dp * dy
       do j = 1, ny
          x = 0.0_pb_dp * dx
          do i = 1, nx
             f(i, j, k, 1) = sin(x)
             f(i, j, k, 2) = sin(y)
             f(i, j, k, 3) = sin(z)
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do

    call div(f, [dx, dy, dz], df)

    rms = 0.0_pb_dp
    z = 0.5_pb_dp * dz
    do k = 1, nz
       y = 0.5_pb_dp * dy
       do j = 1, ny
          x = 0.5_pb_dp * dx
          do i = 1, nx
             expect = cos(x) + cos(y) + cos(z)

             rms = rms + (df(i, j, k) - expect)**2
             
             x = x + dx
          end do
          y = y + dy
       end do
       z = z + dz
    end do
    rms = sqrt(rms / nx / ny / nz)
    
    test_pass = (rms <= 1.0e-9)
    test_pass = (.not. (rms /= rms)) .and. test_pass
    if (.not. test_pass) then
       print *, "FAIL: RMS div(dx) = ", rms, "variable f"
       passing = .false.
    else
       print *, "PASS: RMS div(dx) = ", rms, "variable f"
    end if

  end subroutine check_varying_field
  
end program test_div_3d
