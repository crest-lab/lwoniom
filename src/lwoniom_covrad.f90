!================================================================================!
! This file is part of lwoniom.
!
! Copyright (C) 2023 Patryk Wesolowski, Philipp Pracht
!
! lwoniom is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! lwoniom is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with lwoniom. If not, see <https://www.gnu.org/licenses/>.
!================================================================================!
module lwoniom_covrad
!************************************************************
!* This module implements the ONIOM energy and
!* gradient reconstruction.
!************************************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  implicit none
  private

  real(wp),parameter :: autoaa = 0.52917726_wp
  real(wp),parameter :: aatoau = 1.0_wp/autoaa

  !> covalent radii
  public :: covalent_radius

!&<
  !> covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  !  188-197), values for metals decreased by 10 %
  real(wp),parameter :: covalent_radius(118) = aatoau * [ &
  & 0.32_wp,                                                   0.46_wp, & ! H,He
  & 1.20_wp,0.94_wp,   0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
  & 1.40_wp,1.25_wp,   1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
  & 1.76_wp,1.54_wp,                                                    & ! K,Ca
  &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp,            & ! Sc-
  &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp,            & ! -Zn
  &                    1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
  & 1.89_wp,1.67_wp,                                                    & ! Rb,Sr
  &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp,            & ! Y-
  &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp,            & ! -Cd
  &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp,    & ! In-Xe
  & 2.09_wp,1.76_wp,                                                    & ! Cs,Ba
  &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp,    & ! La-Eu
  &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp,    & ! Gd-Yb
  &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp,            & ! Lu-
  &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp,            & ! -Hg
  &                    1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
  & 2.01_wp,1.81_wp,                                                    & ! Fr,Ra
  &      1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp,       & ! Ac-Am
  &      1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp,       & ! Cm-No
  &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp,            & ! Lr-
  &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp,            & ! -Cn
  &                    1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og
!&>

  public :: lwoniom_rcov_bonds
  public :: link_ratio_k

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_rcov_bonds(nat,at,xyz,factor,bond_tmp)
!*****************************************************
!* Determine a "mock-up" connectivity information
!* based on the sum of covalent radii of two atoms:
!* if R_ab <= (rcov_a + rcov_b)*factor, assume a bond
!******************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: factor
    !> OUTPUT
    integer,intent(out) :: bond_tmp(nat,nat)
    !> LOCAL
    integer :: i,j,k
    real(wp) :: rab,rcovsum

    bond_tmp(:,:) = 0
    ! Loop over all atom pairs (a, b), calculate their distance, and check against
    ! the sum of their covalent radii * factor
    do i = 1,nat
      do j = i+1,nat
        ! Calculate the distance between atoms 'i' and 'j'
        rab = sqrt(sum((xyz(:,i)-xyz(:,j))**2)) ! Loop over all atom pairs (a,b), calculate their distance, and check against the sum of their covalent radii* factor

        ! Determine the sum of covalent radii for atoms 'i' and 'j'
        rcovsum = covalent_radius(at(i))+covalent_radius(at(j))

        ! Check if the distance is smaller than the sum of covalent radii times
        ! 'factor'
        if (rab <= rcovsum*factor) then
          bond_tmp(i,j) = 1  ! Set a bond between atoms 'i' and 'j'
          bond_tmp(j,i) = 1  ! (assuming bond connectivity is symmetric)
          ! write (*,*) i,j,1
        end if
      end do
    end do

  end subroutine lwoniom_rcov_bonds

!========================================================================================!

  function link_ratio_k(at_a,at_b,at_l) result(k)
!**************************************************
!* Calculates (rcov_b + rcov_l) / (rcov_a + rcov_b)
!***************************************************
    implicit none
    integer,intent(in) :: at_a,at_b,at_l  !> atomic number of atoms a,b, and l
    real(wp) :: k

    ! Obtain the covalent radii of atoms a, b, and l
    real(wp) :: rcov_a,rcov_b,rcov_l
    rcov_a = covalent_radius(at_a)
    rcov_b = covalent_radius(at_b)
    rcov_l = covalent_radius(at_l)

    ! Calculate the ratio from the covalent radii of the atoms
    k = (rcov_b+rcov_l)/(rcov_a+rcov_b)

  end function link_ratio_k
!========================================================================================!
!========================================================================================!
end module lwoniom_covrad

