!!! src/constants.f90
!
! Constants module of Poissbox, defines constant values used by Poissbox.
!
!! License
!
! SPDX-License-Identifier: BSD-3-Clause

module constants

  implicit none

  private

  integer, parameter, public :: pb_dp = kind(0.0d0) ! Poissbox double precision (must match PETSc REAL)
  
end module constants
