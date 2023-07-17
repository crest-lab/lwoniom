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

!**************************************************************************!
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
    integer,allocatable  :: opos(:)   !> mapping of each atom in the original (topmost) layer
    integer,allocatable  :: at(:)     !> atomic number
    real(wp),allocatable :: xyz(:,:)  !> Cartesian coordinates, also atomic units -> Bohr
    real(wp),allocatable :: grd(:,:,:)
    !> grd should have dimension(3,nat,2)to store two gradients:
    !> one for the layer (high) and one for the parent layer (low)

    !> link atom coordinates
    integer :: nlink = 0
    integer,allocatable :: linkopos(:)  !> corresponds to which atom in original structure?
    integer,allocatable :: linksto(:)   !> links to which atom in this fragment?
    integer,allocatable :: linkat(:)    !> atom type (will be mostly H)
    real(wp),allocatable :: linkxyz(:,:)  !> Cartesian coordinates, in Bohr
    real(wp),allocatable :: linkgrd(:,:,:) !> similar to grd, but for link atoms

    !> embedding information (TODO, for the future)
    !integer :: npoint
    !real(wp),allocatable :: pointc(:)
    !real(wp),allocatable :: pointxyz(:,:)

  contains
    procedure :: deallocate => structure_data_deallocate
    procedure :: add_child => structure_add_child
    procedure :: write => structure_data_write_structure
    procedure :: extract => extract_atoms_and_links
  end type structure_data
!**************************************************************************!

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ROUTINES GO HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO: we need a routine that takes a structure_data object
!       and returns  the variables nat, at, xyz for a newly build
!       structure consisting out of the fragments atoms+link atoms

!> Extracts the atoms and link atoms of the structure_data object into a new
!structure
  subroutine extract_atoms_and_links(self, nat_new, at_new, xyz_new)
    implicit none
    class(structure_data) :: self
    integer, intent(out) :: nat_new    !> Number of atoms in the new structure
    integer, dimension(:), intent(out) :: at_new   !> Atomic numbers of atoms in the new structure
    real(wp), dimension(:,:), intent(out) :: xyz_new !> Cartesian coordinates of atoms in the new structure

    !> Combine the original atoms and link atoms
    nat_new = self%nat + self%nlink
    at_new = [self%at, self%linkat]
    xyz_new = reshape([self%xyz, self%linkxyz], [3,nat_new])

  end subroutine extract_atoms_and_links

  subroutine structure_data_write_structure(self)
    implicit none
    class(structure_data) :: self

  !>   fragment(i)%write()

  end subroutine structure_data_write_structure


! TODO: a second routine should recieve energy and gradients and distribute it
!       into grd and linkgrd accordingly Eq.6+8

!========================================================================================!
  subroutine structure_data_deallocate(self)
!**************************************************
!* Deallocates a given structure_data object
!**************************************************
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
!********************************************************
!* links two structure_data objects as parent and child
!*******************************************************
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

