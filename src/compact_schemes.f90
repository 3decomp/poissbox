module compact_schemes

  use constants
  use tridsol, only: tdma_periodic
  
  implicit none

  private
  public :: grad, grad_1d
  public :: interp, interp_1d
contains

  ! Compute the 3D gradient tensor (staggered) of a field
  !
  ! XXX: Implemented as Z->Y->X (cell->face->edge->vert)
  subroutine grad(f, dx, df)

    real(pb_dp), dimension(:, :, :), intent(in) :: f      ! Field
    real(pb_dp), dimension(3), intent(in) :: dx           ! Grid spacing
    real(pb_dp), dimension(:, :, :, :), intent(out) :: df ! Gradient tensor

    real(pb_dp), dimension(:, :, :, :), allocatable :: dff ! Face-interpolated field / Z gradient
    real(pb_dp), dimension(:, :, :, :), allocatable :: dfe ! Edge-interpolated field / Z + Y gradient

    integer :: i, j, k
    integer :: nx, ny, nz
    
    nx = size(f, 1)
    ny = size(f, 2)
    nz = size(f, 3)
    
    ! Z gradient/cell->face interpolation
    allocate(dff, mold=df)
    do j = 1, ny
       do i = 1, nx
          call interp_1d(f(i, j, :), dff(i, j, :, 1))
          dff(i, j, :, 2) = dff(i, j, :, 1)
          call grad_1d(f(i, j, :), dx(3), dff(i, j, :, 3))
       end do
    end do

    ! Y gradient/face->edge interpolation
    allocate(dfe, mold=dff)
    do k = 1, nz
       do i = 1, nx
          call interp_1d(dff(i, :, k, 1), dfe(i, :, k, 1))
          call grad_1d(dff(i, :, k, 2), dx(2), dfe(i, :, k, 2))
          call interp_1d(dff(i, :, k, 3), dfe(i, :, k, 3))
       end do
    end do
    deallocate(dff)

    ! X gradient/edge->vert interpolation
    do k = 1, nz
       do j = 1, ny
          call grad_1d(dfe(:, j, k, 1), dx(1), df(:, j, k, 1))
          call interp_1d(dfe(:, j, k, 2), df(:, j, k, 2))
          call interp_1d(dfe(:, j, k, 3), df(:, j, k, 3))
       end do
    end do
    
  end subroutine grad

  ! Compute the 3D interpolation of a field
  !
  ! XXX: Implemented as Z->Y->X (cell->face->edge->vert)
  subroutine interp(f, fi)

    real(pb_dp), dimension(:, :, :), intent(in) :: f   ! Field
    real(pb_dp), dimension(:, :, :), intent(out) :: fi ! Interpolated field

    real(pb_dp), dimension(:, :, :), allocatable :: ff ! Face-interpolated field
    real(pb_dp), dimension(:, :, :), allocatable :: fe ! Edge-interpolated field

    integer :: i, j, k
    integer :: nx, ny, nz
    
    nx = size(f, 1)
    ny = size(f, 2)
    nz = size(f, 3)
    
    ! Cell->face interpolation
    allocate(ff, mold=f)
    do j = 1, ny
       do i = 1, nx
          call interp_1d(f(i, j, :), ff(i, j, :))
       end do
    end do

    ! Face->edge interpolation
    allocate(fe, mold=ff)
    do k = 1, nz
       do i = 1, nx
          call interp_1d(ff(i, :, k), fe(i, :, k))
       end do
    end do
    deallocate(ff)

    ! Edge->vert interpolation
    do k = 1, nz
       do j = 1, ny
          call interp_1d(fe(:, j, k), fi(:, j, k))
       end do
    end do
    deallocate(fe)
    
  end subroutine interp
  
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
