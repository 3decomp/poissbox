module compact_schemes

  use constants
  use tridsol, only: tdma_periodic
  
  implicit none

  private
  public :: grad, grad_1d
  public :: interp, interp_1d
  public :: div, div_1d
  public :: interp_div, interp_1d_div
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
  subroutine interp(f, fi, opt_stagger)

    real(pb_dp), dimension(:, :, :), intent(in) :: f   ! Field
    real(pb_dp), dimension(:, :, :), intent(out) :: fi ! Interpolated field
    integer, intent(in), optional :: opt_stagger
    
    real(pb_dp), dimension(:, :, :), allocatable :: ff ! Face-interpolated field
    real(pb_dp), dimension(:, :, :), allocatable :: fe ! Edge-interpolated field

    integer :: i, j, k
    integer :: nx, ny, nz

    integer :: stagger

    if (present(opt_stagger)) then
       stagger = opt_stagger
    else
       stagger = -1
    end if
    
    nx = size(f, 1)
    ny = size(f, 2)
    nz = size(f, 3)
    
    ! Cell->face interpolation
    allocate(ff, mold=f)
    do j = 1, ny
       do i = 1, nx
          call interp_1d(f(i, j, :), ff(i, j, :), opt_stagger)
       end do
    end do

    ! Face->edge interpolation
    allocate(fe, mold=ff)
    do k = 1, nz
       do i = 1, nx
          call interp_1d(ff(i, :, k), fe(i, :, k), opt_stagger)
       end do
    end do
    deallocate(ff)

    ! Edge->vert interpolation
    do k = 1, nz
       do j = 1, ny
          call interp_1d(fe(:, j, k), fi(:, j, k), opt_stagger)
       end do
    end do
    deallocate(fe)
    
  end subroutine interp

  subroutine interp_div(f, fi)

    real(pb_dp), dimension(:, :, :), intent(in) :: f   ! Field
    real(pb_dp), dimension(:, :, :), intent(out) :: fi ! Interpolated field

    ! Interpolate using forward (vertex->cell) staggering
    call interp(f, fi, +1)
    
  end subroutine interp_div
  
  ! Compute the 1D gradient (staggered) of a field, default backward stagger (cell->vertex)
  subroutine grad_1d(f, dx, df, opt_stagger)

    real(pb_dp), dimension(:), intent(in) :: f   ! Field
    real(pb_dp), intent(in) :: dx                ! Grid spacing
    real(pb_dp), dimension(:), intent(out) :: df ! Gradient
    integer, intent(in), optional :: opt_stagger
    
    real(pb_dp) :: a, b, alpha ! Scheme parameters

    real(pb_dp), dimension(:), allocatable :: ld, d, ud ! Tridiagonal coefficients

    integer :: n

    integer :: stagger
    
    if (present(opt_stagger)) then
       stagger = opt_stagger
    else
       stagger = -1
    end if
    
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
    call eval_1d_rhs(a, b, -1, stagger, f, df)
    
    !! Solve system
    call tdma_periodic(ld, d, ud, df)
    
    !! Cleanup
    deallocate(ld)
    deallocate(d)
    deallocate(ud)
    
  end subroutine grad_1d

  ! Compute the divergence of a vector field
  subroutine div(f, dx, df)

    real(pb_dp), dimension(:, :, :, :), intent(in) :: f ! Vector field
    real(pb_dp), dimension(3), intent(in) :: dx
    real(pb_dp), dimension(:, :, :), intent(out) :: df ! Divergence

    real(pb_dp), dimension(:, :, :, :), allocatable :: dfe
    real(pb_dp), dimension(:, :, :, :), allocatable :: dff
    real(pb_dp), dimension(:), allocatable :: dfc

    integer :: i, j, k
    integer :: nx, ny, nz
    
    nx = size(f, 1)
    ny = size(f, 2)
    nz = size(f, 3)

    ! DDX (Vert->Edge)
    allocate(dfe, mold=f)
    do k = 1, nz
       do j = 1, ny
          call div_1d(f(:, j, k, 1), dx(1), dfe(:, j, k, 1))
          call interp_1d_div(f(:, j, k, 2), dfe(:, j, k, 2))
          call interp_1d_div(f(:, j, k, 3), dfe(:, j, k, 3))
       end do
    end do
    
    ! DDY (Edge->Face)
    allocate(dff, mold=dfe)
    do k = 1, nz
       do i = 1, nx
          call interp_1d_div(dfe(i, :, k, 1), dff(i, :, k, 1))
          call div_1d(dfe(i, :, k, 2), dx(2), dff(i, :, k, 2))
          call interp_1d_div(dfe(i, :, k, 3), dff(i, :, k, 3))
       end do
    end do
    deallocate(dfe)
    
    ! DDZ (Face->Cell)
    allocate(dfc(nz))
    do j = 1, ny
       do i = 1, nx
          call interp_1d_div(dff(i, j, :, 1) + dff(i, j, :, 2), dfc)
          call div_1d(dff(i, j, :, 3), dx(3), df(i, j, :))
          df(i, j, :) = df(i, j, :) + dfc(:)
       end do
    end do
    
    deallocate(dff)
    deallocate(dfc)
  end subroutine div

  ! Compute the divergence of a 1-D field (forward stagger: vertex->cell)
  subroutine div_1d(f, dx, df)

    real(pb_dp), dimension(:), intent(in) :: f   ! Field
    real(pb_dp), intent(in) :: dx                ! Grid spacing
    real(pb_dp), dimension(:), intent(out) :: df ! Gradient

    call grad_1d(f, dx, df, +1)
    
  end subroutine div_1d

  ! Interpolate a 1-D field, default backwards stagger (cell->vertex)
  subroutine interp_1d(f, fi, opt_stagger)

    real(pb_dp), dimension(:), intent(in) :: f   ! Field
    real(pb_dp), dimension(:), intent(out) :: fi ! Interpolated field
    integer, intent(in), optional :: opt_stagger
    
    real(pb_dp) :: a, b, alpha ! Scheme parameters

    real(pb_dp), dimension(:), allocatable :: ld, d, ud ! Tridiagonal coefficients

    integer :: n

    integer :: stagger

    if (present(opt_stagger)) then
       stagger = opt_stagger
    else
       stagger = -1
    end if
    
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
    call eval_1d_rhs(a, b, +1, stagger, f, fi)
    
    !! Solve system
    call tdma_periodic(ld, d, ud, fi)
    
    !! Cleanup
    deallocate(ld)
    deallocate(d)
    deallocate(ud)
    
  end subroutine interp_1d

  ! Interpolate a 1-D field (forward stagger: vertex->cell)
  subroutine interp_1d_div(f, fi)

    real(pb_dp), dimension(:), intent(in) :: f   ! Field
    real(pb_dp), dimension(:), intent(out) :: fi ! Interpolated field

    call interp_1d(f, fi, +1)
    
  end subroutine interp_1d_div

  ! General function for evaluating RHS of staggered compact scheme.
  pure subroutine eval_1d_rhs(a, b, opsign, stagger, f, rhs)

    real(pb_dp), intent(in) :: a, b ! Scheme parameters
    integer, intent(in) :: opsign   ! Set the sign of the finite difference scheme
                                    ! (-1 difference, +1 interpolation)
    integer, intent(in) :: stagger  ! Set the direction of the staggering operation 
                                    ! (-1 from cells->vertices, +1 from vertices->cells)
    real(pb_dp), dimension(:), intent(in) :: f    ! Field
    real(pb_dp), dimension(:), intent(out) :: rhs ! RHS

    integer :: i
    integer :: n

    integer :: shift

    n = size(f)

    ! Set index shift
    if (stagger == -1) then
       shift = 0
    else
       shift = 1
    end if
    
    if (stagger == -1) then
       rhs(1) = a * (f(1) + opsign * f(n)) + b * (f(2) + opsign * f(n - 1))
       rhs(2) = a * (f(2) + opsign * f(1)) + b * (f(3) + opsign * f(n))
    else
       rhs(1) = a * (f(2) + opsign * f(1)) + b * (f(3) + opsign * f(n))
    end if
    do i = 3 - shift, n - 1 - shift
       rhs(i) = a * (f(i + shift) + opsign * f(i - 1 + shift)) + b * (f(i + 1 + shift) + opsign * f(i - 2 + shift))
    end do
    if (stagger == -1) then
       rhs(n) = a * (f(n) + opsign * f(n - 1)) + b * (f(1) + opsign * f(n - 2))
    else
       rhs(n - 1) = a * (f(n) + opsign * f(n - 1)) + b * (f(1) + opsign * f(n - 2))
       rhs(n) = a * (f(1) + opsign * f(n)) + b * (f(2) + opsign * f(n - 1))
    end if
    
  end subroutine eval_1d_rhs
  
end module compact_schemes
