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
module lwoniom_structures
!**************************************************************
!* This module implements the structure_data type.
!* Each structure_data tracks coordinate and layer information
!**************************************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  implicit none
  private

  
!> a single structure_data required for a lwONIOM calculation 
  public :: structure_data 
  type :: structure_data

    !> layer tracker
    integer :: layer = 0

    !> system coordinates
    integer  :: nat = 0
    integer,allocatable  :: at(:)
    real(wp),allocatable :: xyz(:,:)

    !> link atom coordinates
    integer :: nlink = 0
    integer,allocatable :: linkat(:)
    real(wp),allocatable :: linkxyz(:,:)

    !> embedding information (TODO)
    !integer :: npoint
    !real(wp),allocatable :: pointc(:)
    !real(wp),allocatable :: pointxyz(:,:)

  contains
    procedure :: deallocate => structure_data_deallocate
  end type structure_data

  

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ROUTINES GO HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!========================================================================================!
  subroutine structure_data_deallocate(self)
    implicit none
    class(structure_data) :: self
    self%layer = 0
    self%nat = 0
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    self%nlink = 0
    if (allocated(self%linkat)) deallocate (self%linkat)
    if (allocated(self%linkxyz)) deallocate (self%linkxyz)
    !self%npoint = 0
    !if (allocated(self%pointc)) deallocate (self%pointc)
    !if (allocated(self%pointxyz)) deallocate (self%pointxyz)
  end subroutine structure_data_deallocate

!========================================================================================!
end module lwoniom_structures

