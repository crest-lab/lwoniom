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
module lwoniom_engrad
!************************************************************
!* This module implements the ONIOM energy and
!* gradient reconstruction.
!************************************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use lwoniom_setup
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: lwoniom_singlepoint

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_singlepoint(nat,dat,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    !integer,intent(in)  :: at(nat)    !> atom types
    !real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    logical,intent(in),optional    :: verbose  !> printout activation
    type(lwoniom_data),intent(inout) :: dat  !> collection of lwoniom datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io,i,j,k,kk
    integer :: nchilds
    logical :: pr

    !> printout activation via verbosity
    if (present(verbose)) then
      pr = verbose
    else
      pr = .false. !> (there is close to no printout anyways)
    end if

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    io = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TODO do the ONIOM energy and gradient reconstruction here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Loop through layers, starting from the 2nd lowest (highest number)
    do i = dat%nlayer-1,1,-1

      !> check for fragments belonging to this layer
      do j = 1,dat%nfrag
        if (dat%fragment(j)%layer == i) then

          if (.not.allocated(dat%fragment(j)%gradient_qq)) then
            allocate (dat%fragment(j)%gradient_qq(3,nat))
          end if
          dat%fragment(j)%energy_qq = 0.0d0
          dat%fragment(j)%gradient_qq = 0.0d0

          !> if so, check for children fragments
          if (allocated(dat%fragment(j)%child)) then
            nchilds = size(dat%fragment(j)%child(:),1)
            !> loop over all children nodes
            do k = 1,nchilds
              kk = dat%fragment(j)%child(k)

              if (allocated(dat%fragment(kk)%gradient_qq)) then
                dat%fragment(j)%energy_qq = dat%fragment(kk)%energy_qq &
                &                           -dat%fragment(kk)%energy_low
                dat%fragment(j)%gradient_qq(:,:) = dat%fragment(kk)%gradient_qq(:,:) &
                &                                  -dat%fragment(kk)%gradient_low(:,:)
              else
                dat%fragment(j)%energy_qq = dat%fragment(kk)%energy_high &
                &                           -dat%fragment(kk)%energy_low
                dat%fragment(j)%gradient_qq(:,:) = dat%fragment(kk)%gradient_high(:,:) &
                &                                  -dat%fragment(kk)%gradient_low(:,:)
              end if

            end do

            !> add "high" level of parent node (which must be the same level of theory
            !> as the children node's low level)
            dat%fragment(j)%energy_qq = dat%fragment(j)%energy_qq+dat%fragment(j)%energy_high
            dat%fragment(j)%gradient_qq(:,:) = dat%fragment(j)%gradient_qq(:,:) &
            &                                  +dat%fragment(j)%gradient_high(:,:)
          end if
        end if
      end do

    end do

    if (present(iostat)) then
      iostat = io
    end if

  end subroutine lwoniom_singlepoint

!========================================================================================!

  recursive subroutine engrad_recursion(dat,F,nat)
!*********************************************************************
!* Recursion construction of the MC-ONIOM energy and gradient.
!* The energy (or gradient) of the current node F_i is calculated as:
!*  
!*   F_i = F_i^h + ∑_j(f_j - f_j^l)
!*
!* where F_i^h is the high-level energy/gradient.
!* To this, (corrected) energy/gradient contibutions of all 
!* children nodes f_j, minus their low-level contibution f_j^l.
!* Keep in mind that the high-level of the parent node (F_i^h) refers
!* to the same as the low-level of the child node (f_j^l).
!* f_j is itself calculated by the same equation as F_i, which is where
!* the recursion comes into play. The recusion terminates when
!* there are no child nodes. In this case F_i = F_i^h.
!*********************************************************************
    implicit none
    type(lwoniom_data),intent(inout) :: dat
    integer,intent(in) :: F
    integer,intent(in) :: nat
    integer :: nchilds,k,KK

    if (.not.allocated(dat%fragment(F)%gradient_qq)) then
      allocate (dat%fragment(F)%gradient_qq(3,nat))
    end if
    dat%fragment(F)%energy_qq = dat%fragment(F)%energy_high
    dat%fragment(F)%gradient_qq(:,:) = dat%fragment(F)%gradient_high(:,:)

    !> loop over all children nodes
    if (allocated(dat%fragment(F)%child)) then
      nchilds = size(dat%fragment(F)%child(:),1)
      do k = 1,nchilds
        KK = dat%fragment(F)%child(k) !> the child-node ID
        !==============================================!
        call engrad_recursion(dat,KK,nat) !> RECURSION
        !==============================================!

        dat%fragment(F)%energy_qq = dat%fragment(F)%energy_qq+ &
        &  dat%fragment(KK)%energy_qq-dat%fragment(KK)%energy_low

        dat%fragment(F)%gradient_qq(:,:) = dat%fragment(F)%gradient_qq(:,:)+ &
        & dat%fragment(KK)%gradient_qq(:,:)-dat%fragment(kk)%gradient_low(:,:)

      end do
    else
      !> RECURSION TERMINATION: NO FURTHER CHILDREN NODES
      !> i.e., we take energy_qq/gradient_qq for this node 'as is'
      return
    end if

  end subroutine engrad_recursion

!========================================================================================!
!========================================================================================!
end module lwoniom_engrad

