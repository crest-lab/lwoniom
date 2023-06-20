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
!* Each structure_data tracks coordinate and layer informatio
!**************************************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  implicit none
  private

!> a single structure_data required for a lwONIOM calculation
!> tracking information will be stored in an undirected tree
  public :: structure_data
  type :: structure_data

    !> an identifier (position in the oniom list)
    integer :: id = 0

    !> layer tracker, one level of theory per layer
    integer :: layer = 0
    integer :: parent = 0
    integer,allocatable :: child(:)
    !> energies for level of theory at current layer (high) and parent layer (low)
    real(wp) :: energy_high = 0.0_wp
    real(wp) :: energy_low = 0.0_wp

    !> system coordinates
    integer  :: nat = 0
    integer,allocatable  :: at(:)     !> atomic number
    real(wp),allocatable :: xyz(:,:)  !> also atomic units -> Bohr
    real(wp),allocatable :: grd(:,:,:)
    !> grd should have dimension(3,nat,2)to store two gradients:
    !> one for the layer (high) and one for the parent layer (low)

    !> link atom coordinates
    integer :: nlink = 0
    integer,allocatable :: linkat(:)
    real(wp),allocatable :: linkxyz(:,:)
    real(wp),allocatable :: linkgrd(:,:,:) !> similar to grd, but for link atoms

    !> embedding information (TODO)
    !integer :: npoint
    !real(wp),allocatable :: pointc(:)
    !real(wp),allocatable :: pointxyz(:,:)

  contains
    procedure :: deallocate => structure_data_deallocate
    procedure :: add_child => structure_add_child
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
    self%parent = 0
    if (allocated(self%child)) deallocate (self%child)
    self%nat = 0
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%grd)) deallocate (self%grd)
    self%nlink = 0
    if (allocated(self%linkat)) deallocate (self%linkat)
    if (allocated(self%linkxyz)) deallocate (self%linkxyz)
    if (allocated(self%linkgrd)) deallocate (self%linkgrd)
    !self%npoint = 0
    !if (allocated(self%pointc)) deallocate (self%pointc)
    !if (allocated(self%pointxyz)) deallocate (self%pointxyz)
  end subroutine structure_data_deallocate
!========================================================================================!
  subroutine structure_add_child(self,child)
    implicit none
    class(structure_data) :: self
    type(structure_data)  :: child
    integer :: n,n1
    integer,allocatable :: tmp(:)
    if (child%parent .ne. 0) then
      write (stderr,'(a,i0,a)') 'warning: lwONIOM structure of layer ',child%layer, &
      & ' was already associated with another parent layer!'
    end if
    child%parent = self%id
    if (.not.allocated(self%child)) then  !> first child
      allocate (self%child(1))
      self%child(1) = child%id
    else !> new child
      n = size(self%child,1)
      n1 = n+1
      allocate (tmp(n1))
      tmp(1:n) = self%child
      tmp(n1) = child%id
      call move_alloc(tmp,self%child)
    end if
  end subroutine structure_add_child

!========================================================================================!
end module lwoniom_structures

