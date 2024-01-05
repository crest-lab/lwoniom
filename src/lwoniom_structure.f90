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
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  implicit none
  private
  logical,parameter,private :: debug = .false.

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
    real(wp) :: energy_high = 0.0_wp  !> high = QM
    real(wp) :: energy_low = 0.0_wp   !> low  = MM/CG
    real(wp) :: energy_qq = 0.0_wp
    integer :: truenat
    !> gradient storage 
    real(wp),allocatable :: gradient_high(:,:) !> gradient in the original system's basis
    real(wp),allocatable :: gradient_low(:,:)  !> gradient in the original system's basis
    real(wp),allocatable :: gradient_qq(:,:)   !> constructed gradient 
    !> Hessian storage. Note: all Hessians are stored in 1D packed form
    real(wp),allocatable :: Hss_high(:) 
    real(wp),allocatable :: Hss_low(:)  
    real(wp),allocatable :: Hss_qq(:)   
    !> projected via Jacobian
    real(wp),allocatable :: Jaco(:,:)
    
    !> system coordinates
    integer  :: nat = 0
    integer,allocatable  :: opos(:)   !> mapping of each atom in the original (topmost) layer
    integer,allocatable  :: at(:)     !> atomic number
    logical :: replace_at = .false.
    integer,allocatable  :: at_high(:)
    integer,allocatable  :: at_low(:)
    real(wp),allocatable :: xyz(:,:)  !> Cartesian coordinates, also atomic units -> Bohr
    real(wp),allocatable :: grd_high(:,:)
    real(wp),allocatable :: grd_low(:,:)

    !> link atom coordinates
    integer :: nlink = 0
    integer,allocatable :: linkopos(:)    !> corresponds to which atom in original structure?
    integer,allocatable :: linksto(:)     !> links to which atom in this fragment?
    real(wp),allocatable :: link_k(:)     !> link model scaling parameter k
    integer,allocatable :: linkat(:)      !> atom type (will be mostly H)
    real(wp),allocatable :: linkxyz(:,:)  !> Cartesian coordinates, in Bohr
    real(wp),allocatable :: linkgrd_high(:,:)  !> similar to grd, but for link atoms
    real(wp),allocatable :: linkgrd_low(:,:)  !> similar to grd, but for link atoms

    !> Additional user-defined information for the subsystem
    integer,allocatable :: chrg
    integer,allocatable :: uhf

    !> embedding information (TODO, for the future)
    integer :: npoint = 0
    real(wp),allocatable :: pointc(:)
    real(wp),allocatable :: pointxyz(:,:)

  contains
    procedure :: deallocate => structure_data_deallocate
    procedure :: add_child => structure_add_child
    procedure :: extract => extract_atoms_and_links
    procedure :: gradient_distribution
    procedure :: project_gradient_init
    procedure :: allocate_link => allocating_linking_atoms
    procedure :: set_link => set_linking_atoms
    procedure :: dump_fragment !> for testing
    procedure :: jacobian => project_gradient_highlow
    procedure :: update => update_fragment
    procedure :: dump_bin => lwoniom_structure_dump_bin
    procedure :: read_bin => lwoniom_structure_read_bin
    procedure :: construct_jacobian
  end type structure_data
!**************************************************************************!

  !>--- Element symbols
!&<
  character(len=2),parameter :: PSE(118) = [ &
 & 'H ',                                                                                'He', &
 & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
 & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
 & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
 & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
 & 'Cs','Ba','La',                                                                            &
 &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
 &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
 & 'Fr','Ra','Ac',                                                                            &
 &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
 &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]

  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: angstrom = 1.0_wp / bohr
  real(wp),parameter :: autoaa = bohr
  real(wp),parameter :: aatoau = angstrom
!&>
  public :: zsym_to_at,lowercase,PSE

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
    !> Number of atoms in the new structure
    integer,intent(out) :: nat_new
    !> Atomic numbers of atoms in the new   structure
    integer,allocatable,dimension(:),intent(out) :: at_new
    !> Cartesian coordinates of atoms in  the new structure
    real(wp),allocatable,dimension(:,:),intent(out) :: xyz_new

    !> Combine the original atoms and link atoms
    nat_new = self%nat+self%nlink
    allocate (at_new(nat_new))
    allocate (xyz_new(3,nat_new))
    at_new = [self%at,self%linkat]
    xyz_new = reshape([self%xyz,self%linkxyz], [3,nat_new])

  end subroutine extract_atoms_and_links

!========================================================================================!

  subroutine gradient_distribution(self,energy_high,grd_high,energy_low,grd_low)
!****************************************************
!* Distribute gradient of the model system to the
!* structure_data type.
!****************************************************
    implicit none
    class(structure_data) :: self
    real(wp),intent(in) :: energy_high,energy_low
    real(wp),dimension(:,:),intent(in) :: grd_high
    real(wp),dimension(:,:),intent(in) :: grd_low
    integer :: nat_new

    self%energy_high = energy_high
    !> Distribute the gradient to grd and linkgrd
    nat_new = size(grd_high,2)
    self%grd_high = grd_high(:,1:self%nat)
    self%linkgrd_high = grd_high(:,self%nat+1:nat_new)

    self%energy_low = energy_low
    !> Distribute the gradient to grd and linkgrd
    nat_new = size(grd_low,2)
    self%grd_low = grd_low(:,1:self%nat)
    self%linkgrd_low = grd_low(:,self%nat+1:nat_new)

  end subroutine gradient_distribution

!========================================================================================!


  subroutine construct_jacobian(self,fragnat,truenat)
!*********************************************************************
!* The Jacobian is a 3m x 3n matrix ( Jaco(3*m,3*n) ), where
!* m = nat+nlink of the fragment and
!* n = nat of the full system
!* assuming J_ij is a 3x3 diagonal matrix and o(i) is the mapping
!* of the atom index i to the original system, then
!*
!*                                         ⎧ k   if i==linkatom and j==linkopos(i)
!*      ⎛ J_11 ... J_1n ⎞                  ⎪ 1-k if i==linkatom and j==linksto(i)
!*  J = ⎜  ⋮   ⋱   ⋮    ⎟  with J_ij = E * ⎨ 1   if o(i) == j
!*      ⎝ J_m1 ... J_mn ⎠                  ⎩ 0   else
!*
!* k is saved in self%link_k(:), which is our link ratio from covrad
!* Since k is constant in this implementation (determined from covalent rad)
!* the corresponding J elements are easy. If k is calculated dynamically,
!* additional derivatives would be introduced.
!*********************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: fragnat
    integer,intent(in) :: truenat
!&<
    real(wp),parameter :: E(3, 3) = reshape( &
                & [ 1.0_wp, 0.0_wp, 0.0_wp,  &
                &   0.0_wp, 1.0_wp, 0.0_wp,  &
                &   0.0_wp, 0.0_wp, 1.0_wp], &
                & shape(E))
!&>
    real(wp),allocatable :: Jaco(:,:)
    real(wp) :: Jij(3,3),x
    real(wp),allocatable :: g(:)
    integer :: m,n,i,j,l,k,i2,i3,i4,ii4

!>--- reset
    self%truenat = truenat

!>--- allocate J and g'
    m = fragnat != self%nat+self%nlink
    n = truenat
    if (.not.allocated(self%Jaco)) then
      allocate (Jaco(3*m,3*n),source=0.0_wp)
!>--- Set up the Jacobian
      do i = 1,m
        do j = 1,n
          x = 0.0d0
          if (i <= self%nat) then !> "true" atoms
            if (self%opos(i) .eq. j) then
              x = 1.0d0
            else
              x = 0.0d0
            end if
          else !> link atoms
            i2 = i-self%nat
            i3 = self%linkopos(i2)
            i4 = self%linksto(i2)
            ii4 = self%opos(i4)
            if (j == i3) then  !> linkopos
              x = self%link_k(i2)
            elseif (j == ii4) then !> linksto
              x = 1.0d0-self%link_k(i2)
            else
              x = 0.0d0
            end if
          end if
          Jij(:,:) = E(:,:)*x

!>--- Putting Jij into Jaco
          l = (i-1)*3+1
          k = (j-1)*3+1
          Jaco(l:l+2,k:k+2) = Jij(:,:)
        end do
      end do
!>--- we save the Jacobian to self so we don't need to set it up again
!>--- this may only be done if the Jacobian is static
      call move_alloc(Jaco,self%Jaco)
    end if
  end subroutine construct_jacobian

!========================================================================================!

  subroutine project_gradient(self,fragnat,fraggrd,truenat,truegrd)
!*********************************************************************
!* The Jacobian is a 3m x 3n matrix ( Jaco(3*m,3*n) ), where
!* m = nat+nlink of the fragment and
!* n = nat of the full system
!*
!* The gradient of the model system g' can then be projected into the
!* basis of the real system  g = g'J
!*********************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: fragnat
    real(wp),intent(in) :: fraggrd(3,fragnat)
    integer,intent(in) :: truenat
    real(wp),intent(out) :: truegrd(3,truenat)
!&<
    real(wp),parameter :: E(3, 3) = reshape( &
                & [ 1.0_wp, 0.0_wp, 0.0_wp,  &
                &   0.0_wp, 1.0_wp, 0.0_wp,  &
                &   0.0_wp, 0.0_wp, 1.0_wp], &
                & shape(E))
!&>
    real(wp),allocatable :: Jaco(:,:)
    real(wp) :: Jij(3,3),x
    real(wp),allocatable :: g(:)
    integer :: m,n,i,j,l,k,i2,i3,i4,ii4

!>--- reset
    truegrd = 0.d0
    self%truenat = truenat

!>--- allocate J and g'
    m = fragnat != self%nat+self%nlink
    n = truenat
    allocate (g(3*m),source=0.0_wp)
    g = reshape(fraggrd, [3*m])

    if (.not.allocated(self%Jaco)) then
       call self%construct_jacobian(m,n)
    end if

!>--- calculate g = g'J  and save to self%gradient
    truegrd(:,:) = reshape(matmul(g,self%Jaco), [3,truenat])

    deallocate (g)
  end subroutine project_gradient

!========================================================================================!

  subroutine project_gradient_highlow(self,truenat)
!*******************************************************
!* Project both the high and low level gradients
!* of the fragment into dimensions of the "true" system
!*******************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: truenat
    integer :: m,n,i,j,l,k,i2,i3
    integer :: fragnat
    real(wp),allocatable :: g_high(:,:),g_low(:,:)

    !> determine the total number of atoms in fragment
    fragnat = self%nat+self%nlink

    if (fragnat <= truenat) then
      if (allocated(self%grd_high)) then
        !> allocate and call project_gradient for high level
        if (.not.allocated(self%gradient_high)) allocate (self%gradient_high(3,truenat))
        call project_gradient(self,fragnat, &
        &  reshape([self%grd_high,self%linkgrd_high], [3,fragnat]), &
        &  truenat,self%gradient_high)
      end if

      if (allocated(self%grd_low)) then
        !> allocate and call project_gradient for low level
        if (.not.allocated(self%gradient_low)) allocate (self%gradient_low(3,truenat))
        call project_gradient(self,fragnat, &
        &  reshape([self%grd_low,self%linkgrd_low], [3,fragnat]), &
        &  truenat,self%gradient_low)
      end if
    end if
  end subroutine project_gradient_highlow

!========================================================================================!

  subroutine project_gradient_init(self,truenat)
!****************************************************************
!* Do a dry run of projecting both the high and low level
!* gradients of the fragment into dimensions of the "true" system
!*****************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: truenat
    integer :: m,n,i,j,l,k,i2,i3
    integer :: fragnat
    real(wp),allocatable :: g_high(:,:),g_low(:,:)

    !> determine the total number of atoms in fragment
    fragnat = self%nat+self%nlink

    if (fragnat <= truenat) then
      if (.not.allocated(self%grd_high)) then
        allocate (self%grd_high(3,self%nat),source=0.0_wp)
      end if
      if (.not.allocated(self%linkgrd_high).and.self%nlink > 0) then
        allocate (self%linkgrd_high(3,self%nlink),source=0.0_wp)
      end if
      !> allocate and call project_gradient for high level
      if (.not.allocated(self%gradient_high)) allocate (self%gradient_high(3,truenat))
      call project_gradient(self,fragnat, &
      &  reshape([self%grd_high,self%linkgrd_high], [3,fragnat]), &
      &  truenat,self%gradient_high)

      if (.not.allocated(self%grd_low)) then
        allocate (self%grd_low(3,self%nat),source=0.0_wp)
      end if
      if (.not.allocated(self%linkgrd_low).and.self%nlink > 0) then
        allocate (self%linkgrd_low(3,self%nlink),source=0.0_wp)
      end if

      !> allocate and call project_gradient for low level
      if (.not.allocated(self%gradient_low)) allocate (self%gradient_low(3,truenat))
      call project_gradient(self,fragnat, &
      &  reshape([self%grd_low,self%linkgrd_low], [3,fragnat]), &
      &  truenat,self%gradient_low)

      if (.not.allocated(self%gradient_qq)) then
        allocate (self%gradient_qq(3,truenat),source=0.0_wp)
      end if

    end if
  end subroutine project_gradient_init

!========================================================================================!

  subroutine update_fragment(self,truenat,xyz,pc)
!**************************************************
!* Updates the fragment geometry (and point charges)
!* upon change of the original system geometry
!**************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: truenat
    real(wp),intent(in) :: xyz(3,truenat)
    real(wp),intent(in),optional :: pc(truenat)
    integer :: i,j,k,l,a,aa

    k = 0
    do i = 1,truenat
      if (any(self%opos(:) .eq. i)) then
!>--- update atoms
        do j = 1,self%nat
          if (self%opos(j) == i) then
            self%xyz(:,j) = xyz(:,i)
          end if
        end do
      else if (any(self%linkopos(:) .eq. i)) then
!>--- update link atoms
        do j = 1,self%nlink
          if (self%linkopos(j) == i) then
            a = self%linksto(j)
            aa = self%opos(a)
            self%linkxyz(:,j) = link_position(xyz(:,i),xyz(:,aa),self%link_k(j))
          end if
        end do
      else if (present(pc)) then
!>--- update point charges (optional)
        k = k+1
        self%pointc(k) = pc(i)
        self%pointxyz(:,:) = xyz(:,:)
      end if
    end do

  end subroutine update_fragment

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
    if (allocated(self%gradient_high)) deallocate (self%gradient_high)
    if (allocated(self%gradient_low)) deallocate (self%gradient_low)
    if (allocated(self%gradient_qq)) deallocate (self%gradient_qq)
    self%nat = 0
    if (allocated(self%at)) deallocate (self%at)
    self%replace_at = .false.
    if (allocated(self%at_low)) deallocate (self%at_low)
    if (allocated(self%at_high)) deallocate (self%at_high)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%grd_high)) deallocate (self%grd_high)
    if (allocated(self%grd_low)) deallocate (self%grd_low)
    if (allocated(self%Hss_high)) deallocate (self%Hss_high)
    if (allocated(self%Hss_low)) deallocate (self%Hss_low)
    if (allocated(self%Hss_qq)) deallocate (self%Hss_qq)
    if (allocated(self%Jaco)) deallocate (self%Jaco)
    self%nlink = 0
    if (allocated(self%linkat)) deallocate (self%linkat)
    if (allocated(self%linkxyz)) deallocate (self%linkxyz)
    if (allocated(self%linkgrd_high)) deallocate (self%linkgrd_high)
    if (allocated(self%linkgrd_low)) deallocate (self%linkgrd_low)
    self%npoint = 0
    if (allocated(self%pointc)) deallocate (self%pointc)
    if (allocated(self%pointxyz)) deallocate (self%pointxyz)
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
      write (stderr,'(a,i0,a)') '**WARNING** lwONIOM structure of layer ',child%layer, &
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
    allocate (self%link_k(m))      !> link model scaling parameter g
    allocate (self%linkat(m))      !> atom type (will be mostly H)
    allocate (self%linkxyz(3,m))   !> Cartesian coordinates, in Bohr
    allocate (self%linkgrd_high(3,m))  !> similar to grd, but for link atoms (high level)
    allocate (self%linkgrd_low(3,m))   !> similar to grd, but for link atoms (low level)

  end subroutine allocating_linking_atoms

!========================================================================================!
  subroutine set_linking_atoms(self,m,nat,at,xyz,linking_atoms)
!********************************************************
!* linking atoms space allocating
!*******************************************************
    use lwoniom_covrad
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: nat,m
    integer,intent(in) :: at(nat)
    integer,intent(in) :: linking_atoms(3,nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer :: i,k,j,jj

    !> link atom coordinates
    if (debug) write (*,*) 'atomtypes and link factors'
    self%nlink = m
    do i = 1,m
      k = linking_atoms(1,i)
      self%linkopos(i) = k
      j = linking_atoms(2,i)
      jj = self%opos(j)
      self%linksto(i) = j
      self%linkat(i) = 1 !> Hydrogen for cuts through single bonds
      if (debug) self%linkat(i) = 2
      if (linking_atoms(3,i) > 1) then
        !> the linking atom is bound to multiple atoms in the fragment
        !> which means this is a bad setup. We set g to one
        self%link_k(i) = 1.0_wp
      else
        !> the regular case, cuts through single bonds
        self%link_k(i) = link_ratio_k(at(k),at(jj),self%linkat(i))
        if (debug) then
          write (*,*) at(k),at(jj),self%linkat(i),self%link_k(i)
        end if
      end if

      self%linkxyz(:,i) = link_position(xyz(:,k),xyz(:,jj),self%link_k(i))
    end do

  end subroutine set_linking_atoms

!========================================================================================!
  function link_position(ra,rb,k) result(rl)
!************************************************************
!* Calculates rl = rb + k(ra – rb)
!* a is the atom that is being replaced by the linking atom l
!* b is the atom l is attached to.
!************************************************************
    implicit none
    real(wp) :: k
    real(wp) :: rl(3)
    real(wp),intent(in) :: ra(3),rb(3)

    rl(:) = rb(:)+k*(ra(:)-rb(:))

  end function link_position

!========================================================================================!
  subroutine dump_fragment(self,xyz)
!**********************************
!* linking atoms space allocating
!**********************************
    implicit none
    class(structure_data) :: self
    real(wp),intent(in),optional :: xyz(:,:)
    character(len=100) :: fname
    integer :: n_tot,i,j,k,l,ich

    write (fname,'(a,i0,a)') 'fragment.',self%id,'.xyz'

    n_tot = self%nat+self%nlink
    open (newunit=ich,file=trim(fname))
    if (present(xyz)) then
      write (ich,*) size(xyz,2)
    else
      write (ich,*) n_tot
    end if
    write (ich,*)

    do i = 1,self%nat
      write (ich,'(a2,3f16.8)') PSE(self%at(i)),self%xyz(:,i)*autoaa
    end do

    do i = 1,self%nlink
      write (ich,'(a2,3f16.8)') PSE(self%linkat(i)),self%linkxyz(:,i)*autoaa
    end do

    if (present(xyz)) then
      do i = 1,size(xyz,2)
        if (.not.any(self%opos(:) .eq. i).and. &
        &  .not.any(self%linkopos(:) .eq. i)) then
          write (ich,'(a2,3f16.8)') 'He',xyz(:,i)*autoaa
        end if
      end do
    end if

    close (ich)
  end subroutine dump_fragment

!========================================================================================!

  function zsym_to_at(zsym) result(at)
    implicit none
    character(len=*) :: zsym
    integer :: at
    integer :: i,j
    at = 0
    do i = 1,size(PSE,1)
      if (trim(lowercase(zsym)) .eq. trim(lowercase(PSE(i)))) at = i
    end do
  end function zsym_to_at

  function lowerCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: lowerCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(high,s(i:i))
      if (ic > 0) sout(i:i) = low(ic:ic)
    end do
    call move_alloc(sout,lowerCase)
  end function lowerCase

!========================================================================================!

  subroutine lwoniom_structure_dump_bin(self,bin)
!***********************************************************************
!* handling of unformatted write for an entire lwoniom_structure object
!* Enables easy acces I/O to the model setup
!***********************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: bin  !> output channel of the binary file
    logical :: bdum
    integer :: idum,i,j,k,l,m,n
    integer :: d1,d2,d3

    !> an identifier (position in the oniom list)
    write (bin) self%id

    !> layer tracker, one level of theory per layer
    write (bin) self%layer
    write (bin) self%parent

    bdum = allocated(self%child)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%child,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%child(i)
      end do
    end if
    !> energies for level of theory at current layer (high) and parent layer (low)
    write (bin) self%energy_high  !> high = QM
    write (bin) self%energy_low   !> low  = MM/CC
    write (bin) self%energy_qq
    write (bin) self%truenat

    !> high level gradient
    bdum = allocated(self%gradient_high)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%gradient_high,1)
      d2 = size(self%gradient_high,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%gradient_high(i,j)
        end do
      end do
    end if
    !> low level gradient
    bdum = allocated(self%gradient_low)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%gradient_low,1)
      d2 = size(self%gradient_low,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%gradient_low(i,j)
        end do
      end do
    end if
    !> reconstruction gradient
    bdum = allocated(self%gradient_qq)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%gradient_qq,1)
      d2 = size(self%gradient_qq,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%gradient_qq(i,j)
        end do
      end do
    end if

    !> projected via Jacobian
    bdum = allocated(self%Jaco)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%Jaco,1)
      d2 = size(self%Jaco,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%Jaco(i,j)
        end do
      end do
    end if

    !> system coordinates
    write (bin) self%nat

    !integer,allocatable  :: opos(:)   !> mapping of each atom in the original (topmost) layer
    bdum = allocated(self%opos)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%opos,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%opos(i)
      end do
    end if

    !integer,allocatable  :: at(:)     !> atomic number
    bdum = allocated(self%at)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%at,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%at(i)
      end do
    end if

    !logical :: replace_at = .false.
    write (bin) self%replace_at

    !integer,allocatable  :: at_high(:)
    bdum = allocated(self%at_high)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%at_high,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%at_high(i)
      end do
    end if

    !integer,allocatable  :: at_low(:)
    bdum = allocated(self%at_low)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%at_low,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%at_low(i)
      end do
    end if

    !real(wp),allocatable :: xyz(:,:)  !> Cartesian coordinates, also atomic units -> Bohr
    bdum = allocated(self%xyz)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%xyz,1)
      d2 = size(self%xyz,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%xyz(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: grd_high(:,:)
    bdum = allocated(self%grd_high)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%grd_high,1)
      d2 = size(self%grd_high,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%grd_high(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: grd_low(:,:)
    bdum = allocated(self%grd_low)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%grd_low,1)
      d2 = size(self%grd_low,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%grd_low(i,j)
        end do
      end do
    end if

    !> link atom coordinates
    !integer :: nlink
    write (bin) self%nlink

    !integer,allocatable :: linkopos(:)    !> corresponds to which atom in original structure?
    bdum = allocated(self%linkopos)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%linkopos,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%linkopos(i)
      end do
    end if

    !integer,allocatable :: linksto(:)     !> links to which atom in this fragment?
    bdum = allocated(self%linksto)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%linksto,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%linksto(i)
      end do
    end if

    !real(wp),allocatable :: link_k(:)     !> link model scaling parameter k
    bdum = allocated(self%link_k)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%link_k,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%link_k(i)
      end do
    end if

    !integer,allocatable :: linkat(:)      !> atom type (will be mostly H)
    bdum = allocated(self%linkat)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%linkat,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%linkat(i)
      end do
    end if

    !real(wp),allocatable :: linkxyz(:,:)  !> Cartesian coordinates, in Bohr
    bdum = allocated(self%linkxyz)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%linkxyz,1)
      d2 = size(self%linkxyz,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%linkxyz(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: linkgrd_high(:,:)  !> similar to grd, but for link atoms
    bdum = allocated(self%linkgrd_high)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%linkgrd_high,1)
      d2 = size(self%linkgrd_high,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%linkgrd_high(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: linkgrd_low(:,:)  !> similar to grd, but for link atoms
    bdum = allocated(self%linkgrd_low)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%linkgrd_low,1)
      d2 = size(self%linkgrd_low,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%linkgrd_low(i,j)
        end do
      end do
    end if

    !> embedding information (TODO, for the future)
    !integer :: npoint = 0
    write (bin) self%npoint

    !real(wp),allocatable :: pointc(:)
    bdum = allocated(self%pointc)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%pointc,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%pointc(i)
      end do
    end if

    !real(wp),allocatable :: pointxyz(:,:)
    bdum = allocated(self%pointxyz)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%pointxyz,1)
      d2 = size(self%pointxyz,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%pointxyz(i,j)
        end do
      end do
    end if

  end subroutine lwoniom_structure_dump_bin

  subroutine lwoniom_structure_read_bin(self,bin)
!***********************************************************************
!* handling of unformatted read for an entire lwoniom_structure object
!* Enables easy acces I/O to the model setup
!***********************************************************************
    implicit none
    class(structure_data) :: self
    integer,intent(in) :: bin  !> output channel of the binary file
    logical :: bdum
    integer :: idum,i,j,k,l,m,n
    integer :: d1,d2,d3
    logical,parameter :: debug = .false.

    call self%deallocate()

    !> an identifier (position in the oniom list)
    read (bin) self%id
    if (debug) write (*,*) self%id

    !> layer tracker, one level of theory per layer
    read (bin) self%layer
    read (bin) self%parent

    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%child'
      read (bin) d1
      allocate (self%child(d1))
      do i = 1,d1
        read (bin) self%child(i)
      end do
    end if
    !> energies for level of theory at current layer (high) and parent layer (low)
    read (bin) self%energy_high  !> high = QM
    read (bin) self%energy_low   !> low  = MM/CC
    read (bin) self%energy_qq
    read (bin) self%truenat

    !> high level gradient
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%gradient_high'
      read (bin) d1
      read (bin) d2
      allocate (self%gradient_high(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%gradient_high(i,j)
        end do
      end do
    end if
    !> low level gradient
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%gradient_low'
      read (bin) d1
      read (bin) d2
      allocate (self%gradient_low(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%gradient_low(i,j)
        end do
      end do
    end if
    !> reconstruction gradient
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%gradient_qq'
      read (bin) d1
      read (bin) d2
      allocate (self%gradient_qq(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%gradient_qq(i,j)
        end do
      end do
    end if

    !> projected via Jacobian
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%Jaco'
      read (bin) d1
      read (bin) d2
      allocate (self%Jaco(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%Jaco(i,j)
        end do
      end do
    end if

    !> system coordinates
    read (bin) self%nat

    !integer,allocatable  :: opos(:)   !> mapping of each atom in the original (topmost) layer
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%opos(d1))
      do i = 1,d1
        read (bin) self%opos(i)
      end do
    end if

    !integer,allocatable  :: at(:)     !> atomic number
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%at(d1))
      do i = 1,d1
        read (bin) self%at(i)
      end do
    end if

    !logical :: replace_at = .false.
    read (bin) self%replace_at

    !integer,allocatable  :: at_high(:)
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%at_high'
      read (bin) d1
      allocate (self%at_high(d1))
      do i = 1,d1
        read (bin) self%at_high(i)
      end do
    end if

    !integer,allocatable  :: at_low(:)
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%at_low'
      read (bin) d1
      allocate (self%at_low(d1))
      do i = 1,d1
        read (bin) self%at_low(i)
      end do
    end if

    !real(wp),allocatable :: xyz(:,:)  !> Cartesian coordinates, also atomic units -> Bohr
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%xyz'
      read (bin) d1
      read (bin) d2
      allocate (self%xyz(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%xyz(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: grd_high(:,:)
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%grd_high'
      read (bin) d1
      read (bin) d2
      allocate (self%grd_high(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%grd_high(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: grd_low(:,:)
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%grd_low'
      read (bin) d1
      read (bin) d2
      allocate (self%grd_low(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%grd_low(i,j)
        end do
      end do
    end if

    !> link atom coordinates
    !integer :: nlink
    read (bin) self%nlink

    !integer,allocatable :: linkopos(:)    !> corresponds to which atom in original structure?
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%linkopos(d1))
      do i = 1,d1
        read (bin) self%linkopos(i)
      end do
    end if

    !integer,allocatable :: linksto(:)     !> links to which atom in this fragment?
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%linksto(d1))
      do i = 1,d1
        read (bin) self%linksto(i)
      end do
    end if

    !real(wp),allocatable :: link_k(:)     !> link model scaling parameter k
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%link_k(d1))
      do i = 1,d1
        read (bin) self%link_k(i)
      end do
    end if

    !integer,allocatable :: linkat(:)      !> atom type (will be mostly H)
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%linkat(d1))
      do i = 1,d1
        read (bin) self%linkat(i)
      end do
    end if

    !real(wp),allocatable :: linkxyz(:,:)  !> Cartesian coordinates, in Bohr
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      read (bin) d2
      allocate (self%linkxyz(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%linkxyz(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: linkgrd_high(:,:)  !> similar to grd, but for link atoms
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%linkgrd_high'
      read (bin) d1
      read (bin) d2
      allocate (self%linkgrd_high(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%linkgrd_high(i,j)
        end do
      end do
    end if

    !real(wp),allocatable :: linkgrd_low(:,:)  !> similar to grd, but for link atoms
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%linkgrd_low'
      read (bin) d1
      read (bin) d2
      allocate (self%linkgrd_low(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%linkgrd_low(i,j)
        end do
      end do
    end if

    !> embedding information (TODO, for the future)
    !integer :: npoint = 0
    read (bin) self%npoint

    !real(wp),allocatable :: pointc(:)
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%pointc'
      read (bin) d1
      allocate (self%pointc(d1))
      do i = 1,d1
        read (bin) self%pointc(i)
      end do
    end if

    !real(wp),allocatable :: pointxyz(:,:)
    read (bin) bdum
    if (bdum) then
      if (debug) write (*,*) 'reading self%pointxyz'
      read (bin) d1
      read (bin) d2
      allocate (self%pointxyz(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%pointxyz(i,j)
        end do
      end do
    end if

  end subroutine lwoniom_structure_read_bin

!========================================================================================!
!========================================================================================!
end module lwoniom_structures

