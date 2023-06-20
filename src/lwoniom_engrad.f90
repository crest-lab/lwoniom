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

  subroutine lwoniom_singlepoint(nat,at,xyz,dat,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    logical,intent(in),optional    :: verbose  !> printout activation
    type(lwoniom_data),intent(inout) :: dat  !> collection of lwoniom datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io
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

    !> singlpoint + gradient call goes here (best would be another module)
    !call lwoniom_eg(  )

    if (present(iostat)) then
      iostat = io
    end if

  end subroutine lwoniom_singlepoint

!========================================================================================!
!========================================================================================!
end module lwoniom_engrad

