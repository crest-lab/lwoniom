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
    integer :: truenat
    real(wp),allocatable :: gradient(:,:) !> gradient in the original system's dimension
    !> projected via Jacobian

    !> system coordinates
    integer  :: nat = 0
    integer,allocatable  :: opos(:)   !> mapping of each atom in the original (topmost) layer
    integer,allocatable  :: at(:)     !> atomic number
    real(wp),allocatable :: xyz(:,:)  !> Cartesian coordinates, also atomic units -> Bohr
    real(wp),allocatable :: grd(:,:)

    !> link atom coordinates
    integer :: nlink = 0
    integer,allocatable :: linkopos(:)    !> corresponds to which atom in original structure?
    integer,allocatable :: linksto(:)     !> links to which atom in this fragment?
    real(wp),allocatable :: link_g(:)     !> link model scaling parameter g
    integer,allocatable :: linkat(:)      !> atom type (will be mostly H)
    real(wp),allocatable :: linkxyz(:,:)  !> Cartesian coordinates, in Bohr
    real(wp),allocatable :: linkgrd(:,:)  !> similar to grd, but for link atoms

    !> embedding information (TODO, for the future)
    !integer :: npoint
    !real(wp),allocatable :: pointc(:)
    !real(wp),allocatable :: pointxyz(:,:)

  contains
    procedure :: deallocate => structure_data_deallocate
    procedure :: add_child => structure_add_child
    procedure :: extract => extract_atoms_and_links
    procedure :: gradient_distribution
    procedure :: allocate_link => allocating_linking_atoms
    procedure :: set_link => set_linking_atoms

  end type structure_data
!**************************************************************************!

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine extract_atoms_and_links(self,nat_new,at_new,xyz_new)
!*****************************************************************
!* Extracts the atoms and link atoms of the structure_data object
!* into a new structure.
!*****************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(out) :: nat_new    !> Number of atoms in the new structure
    integer,allocatable,dimension(:),intent(out) :: at_new   !> Atomic numbers of atoms in the new structure
    real(wp),allocatable,dimension(:,:),intent(out) :: xyz_new !> Cartesian coordinates of atoms in the new structure

    !> Combine the original atoms and link atoms
    nat_new = self%nat+self%nlink
    allocate (at_new(nat_new))
    allocate (xyz_new(3,nat_new))
    at_new = [self%at,self%linkat]
    xyz_new = reshape([self%xyz,self%linkxyz], [3,nat_new])

  end subroutine extract_atoms_and_links

!========================================================================================!

  subroutine gradient_distribution(self,energy,grd)
!****************************************************
!* Distribute gradient of the model system to the
!* structure_data type.
!****************************************************
    implicit none
    class(structure_data) :: self
    real(wp),intent(in) :: energy
    real(wp),dimension(:,:),intent(in) :: grd
    integer :: nat_new

    self%energy_high = energy
    !> Distribute the gradient to grd and linkgrd
    nat_new = size(grd,2)
    self%grd = grd(:,1:self%nat)
    self%linkgrd = grd(:,self%nat+1:nat_new)

  end subroutine gradient_distribution

!========================================================================================!

  subroutine project_gradient(self,truenat,truegrd)
!*********************************************************************
!* The Jacobian is a 3m x 3n matrix ( Jaco(3*m,3*n) ), where
!* m = nat+nlink of the fragment and
!* n = nat of the full system
!* assuming J_ij is a 3x3 diagonal matrix and o(i) is the mapping
!* of the atom index i to the original system, then
!*
!*                                         ⎧ h   if ...
!*      ⎛ J_11 ... J_1n ⎞                  ⎪ 1-h if ...
!*  J = ⎜  ⋮   ⋱   ⋮    ⎟  with J_ij = E * ⎨ 1   if o(i) == j
!*      ⎝ J_m1 ... J_mn ⎠                  ⎩ 0   else
!*
!*
!* The gradient of the model system g' can then be projected into the
!* basis of the real system  g = g'J
!*********************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: truenat
    real(wp),intent(out),optional :: truegrd(3,truenat)
!&<
    real(wp),parameter :: E(3, 3) = reshape( &
                & [ 1.0_wp, 0.0_wp, 0.0_wp,  &
                &   0.0_wp, 1.0_wp, 0.0_wp,  &
                &   0.0_wp, 0.0_wp, 1.0_wp], &
                & shape(E))
!&>
    real(wp),allocatable :: Jaco(:,:)
    real(wp),allocatable :: g(:)
    integer :: m,n

!>--- if self%gradient wasn't allocated so far, do it now
    if (.not.allocated(self%gradient)) then
      allocate (self%gradient(3,truenat),source=0.0_wp)
      self%truenat = truenat
    end if

!>--- allocate J and g'
    m = self%nat+self%nlink
    n = truenat
    allocate (g(3*m),source=0.0_wp)
    g = reshape([self%grd,self%linkgrd], [3*m])
    allocate (Jaco(3*m,3*n),source=0.0_wp)

!>--- Set up the Jacobian

!TODO, set up the Jacobian (Eq. 6+8), see comment above.
! we might need to test around a bit
! h is saved in self%link_g(:)

!>--- calculate g = g'J  and save to self%gradient
    self%gradient = reshape(matmul(g,Jaco), [3,truenat])
    if (present(truegrd)) truegrd = self%gradient

    deallocate (Jaco,g)
  end subroutine project_gradient

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
  subroutine allocating_linking_atoms(self,m)
!********************************************************
!* linking atoms space allocating
!*******************************************************
    implicit none
    class(structure_data) :: self
    integer :: m

    !> link atom coordinates
    self%nlink = m
    allocate (self%linkopos(m))    !> corresponds to which atom in original structure?
    allocate (self%linksto(m))     !> links to which atom in this fragment?
    allocate (self%link_g(m))     !> link model scaling parameter g
    allocate (self%linkat(m))      !> atom type (will be mostly H)
    allocate (self%linkxyz(3,m))  !> Cartesian coordinates, in Bohr
    allocate (self%linkgrd(3,m))  !> similar to grd, but for link atoms

  end subroutine allocating_linking_atoms

!========================================================================================!
  subroutine set_linking_atoms(self,nat,at,xyz,linking_atoms)
!********************************************************
!* linking atoms space allocating
!*******************************************************
    use lwoniom_covrad
    implicit none
    class(structure_data) :: self
    integer :: m,nat
    integer :: at(nat),linking_atoms(2,nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer :: i,k,j

    !> link atom coordinates
    m = self%nlink
    do i = 1,m
      k = linking_atoms(1,i)
      self%linkopos(i) = k
      j = linking_atoms(2,i)
      self%linksto(i) = j
      self%linkat(i) = 1
      self%link_g(i) = link_ratio_g(at(j),at(k),self%linkat(i))
      self%linkxyz(:,i) = link_position(xyz(:,j),xyz(:,k),self%link_g(i))
    end do

  end subroutine set_linking_atoms

!========================================================================================!
  function link_position(ra,rb,g) result(rl)
!**************************************************
!* Calculates rl = rb + g(ra – rb)
!***************************************************
    implicit none
    real(wp) :: g
    real(wp) :: rl(3)
    real(wp),intent(in) :: ra(3),rb(3)

    rl(:) = rb(:)+g*(ra(:)-rb(:))

  end function link_position

!========================================================================================!

end module lwoniom_structures

