program test_grad_1d

  use constants
  use compact_schemes
  
  implicit none

  integer, parameter :: n = 128            ! Problem size (nodes)
  real(pb_dp), parameter :: L = 3.14_pb_dp ! Domain size
  real(pb_dp), parameter :: dx = L / n     ! Grid spacing
  
  real(pb_dp), dimension(:), allocatable :: f  ! The function
  real(pb_dp), dimension(:), allocatable :: df ! The gradient

  logical :: passing

  passing = .true.
  
  call init()

  call check_constant_field()
  
  call fin()

  if (.not. passing) then
     print *, "FAIL"
     stop 1
  end if
  
contains

  subroutine init()

    allocate(f(n))
    allocate(df(n))
    
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
    
    f(:) = 2.8170923_pb_dp ! Arbitrary constant field
    df(:) = 73.29_pb_dp    ! Non-zero (i.e. wrong) gradient
    
    call grad_1d(f, dx, df)

    rms = sqrt(sum(df**2) / n)
    if (rms > (100 * epsilon(rms))) then
       print *, "FAIL: RMS dfdx = ", rms, "(f = const)"
       passing = .false.
    else
       print *, "PASS: dfdx (f = const)"
    end if

    call interp_1d(f, df)

    rms = sqrt(sum((df - f)**2) / n)
    if (rms > (100 * epsilon(rms))) then
       print *, "FAIL: RMS [f] = ", rms, "(f = const)"
       passing = .false.
    else
       print *, "PASS: [f] (f = const)"
    end if
    
  end subroutine check_constant_field
  
end program test_grad_1d
