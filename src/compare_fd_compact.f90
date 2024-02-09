program compare_fd_compact

  use constants
  use compute_lapl, only : evaluate_laplacian_pointwise
  use compact_schemes, only : lapl

  implicit none
  
  !! Grid dimensions - hardcoded for simplicity
  integer, parameter :: nx = 64
  integer, parameter :: ny = 64
  integer, parameter :: nz = 64
  integer, dimension(3), parameter :: n = [nx, ny, nz]

  !! Problem dimension - hardcoded for simplicity
  real(pb_dp), parameter :: pi = 4 * atan(1.0_pb_dp)
  real(pb_dp), parameter :: Lx = 2 * pi
  real(pb_dp), parameter :: Ly = 2 * pi
  real(pb_dp), parameter :: Lz = 2 * pi
  real(pb_dp), parameter :: dx = Lx / nx
  real(pb_dp), parameter :: dy = Ly / ny
  real(pb_dp), parameter :: dz = Lz / nz

  real(pb_dp), dimension(nx, ny, nz) :: f
  real(pb_dp), dimension(nx, ny, nz) :: d2fdx2
  real(pb_dp) :: f_p
  
  integer :: i, j, k
  real(pb_dp) :: x, y, z

  integer :: il, ir, jl, jr, kl, kr
  integer, dimension(3) :: iidx, jidx, kidx

  real(pb_dp) :: rms, linf, l2, delta
  real(pb_dp) :: d2fdx2_pw
  
  z = 0.5_pb_dp * dz
  do k = 1, nz
     y = 0.5_pb_dp * dy
     do j = 1, ny
        x = 0.5_pb_dp * dx
        do i = 1, nx
           f_p = sin(x) + sin(y) + sin(z)
           f(i, j, k) = f_p
           
           x = x + dx
        end do
        y = y + dy
     end do
     z = z + dz
  end do

  call lapl(f, [dx, dy, dz], d2fdx2)

  linf = 0.0_pb_dp
  l2 = 0.0_pb_dp
  do k = 1, nz
     if (k == 1) then
        kl = nz
     else
        kl = k - 1
     end if
     if (k == nz) then
        kr = 1
     else
        kr = k + 1
     end if
     kidx = [kl, k, kr]
     
     do j = 1, ny
        if (j == 1) then
           jl = ny
        else
           jl = j - 1
        end if
        if (j == ny) then
           jr = 1
        else
           jr = j + 1
        end if
        jidx = [jl, j, jr]
        
        do i = 1, nx
           if (i == 1) then
              il = nx
           else
              il = i - 1
           end if
           if (i == nx) then
              ir = 1
           else
              ir = i + 1
           end if
           iidx = [il, i, ir]

           d2fdx2_pw = evaluate_laplacian_pointwise(f(iidx, jidx, kidx), [dx, dy, dz])

           delta = d2fdx2(i, j, k) - d2fdx2_pw
           linf = max(linf, abs(delta))
           l2 = l2 + delta**2
        end do
     end do
  end do

  rms = sqrt(l2 / nx / ny / nz)
  l2 = sqrt(l2)
  print *, linf, l2, rms
  
end program compare_fd_compact
