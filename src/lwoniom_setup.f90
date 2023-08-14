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

module lwoniom_setup
!***************************************************************
!* This module implements the lwoniom_data type,
!* the main object used to track fragments and layers
!**************************************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  use lwoniom_structures
  use lwoniom_covrad
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: lwoniom_data
  public :: lwoniom_initialize

!**************************************************************!
!> this type bundles together most of the
!> data required for a lwONIOM calculation
  type :: lwoniom_data

    real(wp) :: lwoniom_energy

    !> number of layers
    integer :: nlayer = 0
    integer,allocatable :: layer(:)

    !> number of fragments
    integer :: nfrag = 0
    !> list of fragments (subsystems)
    type(structure_data),allocatable :: fragment(:)

    !> further system information
    integer,allocatable :: bond(:,:)
    integer,allocatable :: indexf(:)

  contains
    procedure :: info => print_lwoniom_info
    procedure :: deallocate => lwoniom_data_deallocate
    procedure :: add_fragment => lwoniom_add_fragment
  end type lwoniom_data
!************************************************************!

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!========================================================================================!
  subroutine print_lwoniom_info(self,iunit)
    class(lwoniom_data) :: self
    integer,intent(in),optional :: iunit ! file handle (usually output_unit=6)
    character(len=*),parameter :: outfmt = '(2x,a,f16.8,1x,a)'
    character(len=*),parameter :: outfmt2 = '(2x,a,i16,1x,a)'
    integer :: myunit,i,j,k
    character(len=10) :: atmp
    character(len=200) :: btmp
    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if
    write (myunit,'(2x,13x,a)') 'lwONIOM summary'
    write (myunit,'(2x,a40)') repeat('-',40)
    write (myunit,outfmt2) 'Layers        :',self%nlayer,''
    write (myunit,outfmt2) 'Substructures :',self%nfrag,''

    if (allocated(self%fragment)) then
      write (myunit,*)
      do i = 1,self%nlayer
        write (myunit,'(2x,a,i0)') 'Layer ',i
        do j = 1,self%nfrag
          if (self%fragment(j)%layer == i) then
            write (myunit,'(4x,a,i0)',advance='no') '-> fragment ',self%fragment(j)%id
            if (self%fragment(j)%parent .ne. 0) then
              write (myunit,'(1x,a,i0)') ', substructure of fragment ',self%fragment(j)%parent
            else
              write (myunit,*)
            end if
            flush (myunit)
          end if
        end do
      end do
    end if
  end subroutine print_lwoniom_info

!========================================================================================!
  subroutine lwoniom_initialize(nat,at,xyz,dat, &
  &                 layer,subsystem,bond,       &
  &                 print,verbose,iunit,iostat)
!******************************************************************
!* This routine sets up the lwoniom_data object and handles
!* the partitioning of the whole system into different layers.
!*
!* It requires input coordinates for the entire system, as
!* well as "layering" information.
!* The layering instruction can be provided by one or two arrays:
!*  A) If only the "layer" array is provided, layers will contain
!*     all the atoms corresponting to the respective layer
!*  B) If the "layer" and the "subsystem" arrays are provided,
!*     one layer can have mutliple fragments (MC-ONIOM).
!*     However, the array "bond" must also be present in this case
!*
!******************************************************************
    character(len=*),parameter :: source = 'lwoniom_initialize'
    !> INPUT
    integer,intent(in) :: nat          !> number of atoms
    integer,intent(in) :: at(nat)      !> atom number for each atom
    real(wp),intent(in) :: xyz(3,nat)  !> Cartesian coordinates for each atom
    integer,intent(in) :: layer(nat)   !> layer instruction for each atom
    integer,intent(in),optional :: subsystem(nat) !> subsystem instruction (for MC-ONIOM)
    integer,intent(in),optional :: bond(nat,nat) !> bonding/interaction information
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    !> OUTPUT
    integer,intent(out),optional :: iostat
    type(lwoniom_data),intent(inout) :: dat
    !> LOCAL
    type(structure_data) :: tmp
    integer,allocatable :: subsystem_tmp(:)
    integer,allocatable :: layer_tmp(:)
    integer,allocatable :: bond_tmp(:,:)
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun
    integer :: maxf
    integer,allocatable ::   indexf(:),parent(:)
    integer :: i,j,k

    io = 0

!> mapping of optional instuctions
    if (present(print)) then
      pr = print
    else
      pr = .false.
    end if
    if (present(verbose)) then
      pr2 = verbose
    else
      pr2 = .false.
    end if
    if (pr2) pr = pr2
    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ROUTINE CONTENT GOES FROM HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate (subsystem_tmp(nat),layer_tmp(nat),source=0)
    if (present(subsystem)) then
      subsystem_tmp(:) = subsystem(:)
      if (.not.present(bond)) then
        write (stderr,'(a)') "**ERROR** 'bond' array not provided in call to "//source
        error stop
      end if
    else
      subsystem_tmp(:) = 1
    end if
    layer_tmp(:) = layer(:)
    allocate (indexf(nat))
    !> determine maxf and indexf, i.e. the number of "nodes"
    !> and the mapping which atom belongs to which node
    call count_fragments(nat,layer_tmp,subsystem_tmp,maxf,indexf)
    dat%nlayer = maxval(layer_tmp,1)

    !> from the connectivity determine the tree structure
    allocate (parent(maxf))
    if (maxf .eq. maxval(layer_tmp,1)) then
      call construct_tree_ONIOM_classic(maxf,parent)
    else if (present(bond)) then
      call construct_tree_ONIOM_multicenter(nat,layer_tmp,indexf,bond,maxf,parent)
    else
      write (stderr,'(a)') "**ERROR** 'bond' array not provided in call to "//source
      error stop
    end if

    !> go through all the fragments, create structure_data and put
    !> them into the fragment list
    do i = 1,maxf
      call tmp%deallocate()
      call fragment_set_atoms_ONIOM(nat,at,xyz,layer_tmp,maxf,indexf,parent,i,tmp)
      call dat%add_fragment(tmp)
    end do
    !> and connect them in the parent-child relations
    do i = 1,maxf
      if (parent(i) .ne. 0) then
        j = parent(i)
        call dat%fragment(j)%add_child(dat%fragment(i))
      end if
    end do

    !> For the next steps we need some kind of connectivity information
    allocate (bond_tmp(nat,nat),source=0)
    if (present(bond)) then
      !> if bond was provided, use those
      bond_tmp = bond
    else
      call lwoniom_rcov_bonds(nat,at,xyz,1.1_wp,bond_tmp)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  linking atoms between layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> at this point, dat%fragment should be set up
!> what's missing are the linking atoms between the layers
! Determine the linking atoms between layers
    call determine_linking_atoms(dat,at,xyz,bond_tmp)

    if ((io /= 0).and.pr) then
      write (myunit,'("Could not create lwONIOM object ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if

    if (allocated(subsystem_tmp)) deallocate (subsystem_tmp)
    if (allocated(indexf)) deallocate (indexf)
  end subroutine lwoniom_initialize

!========================================================================================!
  subroutine determine_linking_atoms(self,at,xyz,bond_tmp)
!****************************************
!* Set up link atoms for each fragment
!****************************************
    implicit none
    type(lwoniom_data) :: self
    integer,intent(in) :: at(:) !> atomic numbers of the original system    
    real(wp),intent(in) :: xyz(:,:)  !> Cartesian coordinates of the original system
    integer,intent(in) :: bond_tmp(:,:)  !> connectivity info of the original system
    integer :: i,j,k,l,m,nat,nat_fragment

    integer,allocatable :: linking_atoms(:,:)

    !> nat is the original system's number of atoms
    nat = size(bond_tmp,1) 
    !> linking_atoms is our "work" array. it's second dimension is <nat for each fragment
    allocate( linking_atoms(2,nat) )

    !> loop over all fragments
    do i = 1,self%nfrag

      !> loop over all atoms in the fragments i
      !> number of atoms in fragment i inherited from the original structure
      nat_fragment = self%fragment(i)%nat
      !> m will count the total number of link atoms in fragment i
      m = 0
      !> linking_atoms and m must be reset for each fragment
      linking_atoms = 0
      do j = 1,nat_fragment
        
        !> check all bonds of the atom, l is the atom's index in the original structure
        l = self%fragment(i)%opos(j)
        do k = 1,nat

          !> if there is a bond between atoms l and k, ...
          if (bond_tmp(l,k) .ne. 0) then
 
            !> ... check if k is any of the atoms in fragment i
            if (.not.any(self%fragment(i)%opos(:) .eq. k)) then

            !> if k is NOT a member of the fragment, it must be a new link atom   
            !> increment the number of link atoms
              m = m+1
            !> document the reference atom k for the link
              linking_atoms(1,m) = k
            !> document to which atom it is linked to, using the fragment atom index
              linking_atoms(2,m) = j  !> i.e., j not l (!)
            end if
          end if
        end do
      end do
      !>allocating m linking atoms for fragment i 
      call self%fragment(i)%allocate_link(m)
      call self%fragment(i)%set_link(nat,at,xyz,linking_atoms)
    end do


 



    deallocate( linking_atoms )
  end subroutine determine_linking_atoms
!========================================================================================!
  subroutine lwoniom_data_deallocate(self)
!****************************************************
!* This subroutine deallocates a lwoniom_data object
!*****************************************************
    implicit none
    class(lwoniom_data) :: self
    self%lwoniom_energy = 0.0_Wp
    self%nlayer = 0
    self%nfrag = 0
    if (allocated(self%fragment)) deallocate (self%fragment)
  end subroutine lwoniom_data_deallocate

!========================================================================================!
  subroutine lwoniom_add_fragment(self,frag)
!*****************************************************
!* This subroutine adds a single structure_data object
!* to the fragment list of a lwoniom_data object
!******************************************************
    implicit none
    class(lwoniom_data) :: self
    type(structure_data) :: frag
    type(structure_data),allocatable :: tmp(:)
    integer :: n,n1
    if (frag%layer > self%nlayer) then
      self%nlayer = frag%layer
    end if
    if (self%nfrag == 0) then
      self%nfrag = 1
      allocate (self%fragment(1))
      self%fragment(1) = frag
    else
      n = self%nfrag
      n1 = n+1
      allocate (tmp(n1))
      tmp(1:n) = self%fragment(1:n)
      tmp(n1) = frag
      call move_alloc(tmp,self%fragment)
      self%nfrag = n1
    end if
  end subroutine lwoniom_add_fragment

!========================================================================================!

  subroutine count_fragments(nat,layer,subsystem,maxf,indexf)
!**********************************************************************
!* Given the two tracking arrays "layer" and "subsystem",
!* this routine computes the maximum number of fragments (subsets)
!* and enumerates the atoms accordingly.
!**********************************************************************
    implicit none
    character(len=*),parameter :: source = 'lwoniom_count_fragments'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(inout) :: layer(nat)
    integer,intent(in) :: subsystem(nat)
    !> OUTPUT
    integer,intent(out) :: maxf
    integer,intent(out) :: indexf(nat)
    !> LOCAL
    integer :: i,j,k,l
    integer :: maxl,ninlayer
    integer,allocatable :: tmp(:)
    integer :: ntmp

    !> initialize
    maxf = 0
    indexf(:) = 0
    allocate (tmp(nat))

    !> iterate through the layers
    maxl = maxval(layer,1)
    do i = 1,maxval(layer,1)
      if (i > maxl) cycle !> safety, since maxl could decrease
      !> check if there are objects in the layer
      !> if the layer is empty, modify the layer tracking accordingly
      do while (count(layer(:) .eq. i) <= 0.and.i <= maxl)
        ninlayer = count(layer(:) .eq. i)
        if (ninlayer <= 0.and.i <= maxl) then
          maxl = maxl-1
          do j = 1,nat
            if (layer(j) >= i) layer(j) = layer(j)-1
          end do
          cycle
        end if
      end do

      !> iterate through all members of the layer and count subsystem types
      ntmp = 0
      tmp(:) = 0
      do j = 1,nat
        if (layer(j) == i) then
          if (.not.any(subsystem(j) .eq. tmp)) then
            ntmp = ntmp+1
            tmp(ntmp) = subsystem(j)
          end if
        end if
      end do
      !> and then iterate through the subsystem types for this layer
      !> and track all atoms to belonging to both
      do j = 1,ntmp
        maxf = maxf+1
        do k = 1,nat
          if (layer(k) == i.and.subsystem(k) == tmp(j)) then
            indexf(k) = maxf
          end if
        end do
      end do
    end do

    deallocate (tmp)
  end subroutine count_fragments

!========================================================================================!
  subroutine construct_tree_ONIOM_classic(maxf,parent)
!************************************************************************
!* Set layer hierarchy for classic ONIOMn schemes, which is trivial
!* since there are no overlapping fragments and just one center per layer
!************************************************************************
    integer,intent(in) :: maxf
    integer,intent(out) :: parent(maxf)
    integer :: i,j,k
    do i = 1,maxf
      parent(i) = i-1
    end do
  end subroutine construct_tree_ONIOM_classic

!========================================================================================!
  subroutine construct_tree_ONIOM_multicenter(nat,layer,indexf,bond,maxf,parent)
!****************************************************************
!* Set layer hierarchy for MC-ONIOMn schemes, which needs to be
!* done recursively.
!****************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: layer(nat)
    integer,intent(in) :: indexf(nat)
    integer,intent(in) :: bond(nat,nat)
    integer,intent(in) :: maxf
    !> OUTPUT
    integer,intent(out) :: parent(maxf)
    !> LOCAL
    integer :: maxl
    integer :: i,j,k,l
    parent(:) = 0
    maxl = maxval(layer,1)
    !> go through the layers reversely (i.e. start with the innermost layer)
    do i = maxl,1,-1
      !> go through substructures
      do j = 1,maxf
        do k = 1,nat
          if (indexf(k) == j) then
            !> found a matching atom, now check to what
            !> upper layer/fragment it is connected
            do l = 1,nat
              if (abs(bond(l,k)) > 0.and.layer(l) < i) then
                parent(j) = indexf(l)
              end if
            end do
          else
            cycle
          end if
        end do
      end do
    end do

  end subroutine construct_tree_ONIOM_multicenter

!========================================================================================!
  subroutine fragment_set_atoms_ONIOM(nat,at,xyz,layer,maxf,indexf,parent,k,frag)
!**************************************************************************
!* This subroutine creates a structure_data object frag for the
!* k-th fragment index in an ONIOM setup, based on the present information
!***************************************************************************
    implicit none
    character(len=*),parameter :: source = 'fragment_set_atoms'
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: layer(nat)
    integer,intent(in) :: maxf
    integer,intent(in) :: indexf(nat)
    integer,intent(in) :: parent(maxf)
    integer,intent(in) :: k
    type(structure_data),intent(out) :: frag
    logical,allocatable :: taken(:)
    integer :: i,natf,n,layerk

    call frag%deallocate()
    allocate (taken(nat),source=.false.)
    taken(:) = .false.

    if (k > maxf.or.k < 1) then
      write (stderr,'(a,a)') "**ERROR** invalid fragment index in ",source
      error stop
    end if

    !> add all atoms that belong to k anyways
    do i = 1,nat
      if (indexf(i) == k) then
        taken(i) = .true.
        layerk = layer(i)
      end if
    end do
    !> and then go through tree recursively
    call recursive_take(nat,maxf,parent,indexf,k,taken)

    !> determine how many atoms there are in the fragment total
    natf = count(taken(:),1)
    !> and finally, put the information into frag
    frag%id = k
    frag%layer = layerk
    frag%nat = natf
    allocate (frag%opos(natf),source=0)
    allocate (frag%at(natf),source=0)
    allocate (frag%xyz(3,natf),source=0.0_wp)
    allocate (frag%grd(3,natf),source=0.0_wp)
    n = 0
    do i = 1,nat
      if (taken(i)) then
        n = n+1
        frag%opos(n) = i
        frag%at(n) = at(i)
        frag%xyz(1:3,n) = xyz(1:3,i)
      end if
    end do

    deallocate (taken)
  contains
    recursive subroutine recursive_take(nat,maxf,parent,indexf,k,taken)
      implicit none
      integer,intent(in) :: k
      integer,intent(in) :: nat
      integer,intent(in) :: maxf
      integer,intent(in) :: indexf(nat)
      integer,intent(in) :: parent(maxf)
      logical,intent(inout) :: taken(nat)
      integer :: i,j
      do i = 1,maxf
        if (parent(i) .eq. k) then
          do j = 1,nat
            if (indexf(j) == i) then
              taken(j) = .true.
            end if
          end do
          call recursive_take(nat,maxf,parent,indexf,i,taken)
        end if
      end do
    end subroutine recursive_take
  end subroutine fragment_set_atoms_ONIOM

!========================================================================================!
  subroutine fragment_set_atoms_QMMM(nat,at,xyz,layer,maxf,indexf,parent,k,frag)
!**************************************************************************
!* This subroutine creates a structure_data object frag for the
!* k-th fragment index in an QMMM setup, based on the present information
!***************************************************************************
    implicit none
    character(len=*),parameter :: source = 'fragment_set_atoms_QMMM'
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: layer(nat)
    integer,intent(in) :: maxf
    integer,intent(in) :: indexf(nat)
    integer,intent(in) :: parent(maxf)
    integer,intent(in) :: k
    type(structure_data),intent(out) :: frag
    logical,allocatable :: taken(:)
    integer :: i,natf,n,layerk

    call frag%deallocate()
    allocate (taken(nat),source=.false.)
    taken(:) = .false.

    if (k > maxf.or.k < 1) then
      write (stderr,'(a,a)') "**ERROR** invalid fragment index in ",source
      error stop
    end if

    !> add all atoms that belong to k
    do i = 1,nat
      if (indexf(i) == k) then
        taken(i) = .true.
        layerk = layer(i)
      end if
    end do

    !> determine how many atoms there are in the fragment total
    natf = count(taken(:),1)
    !> and finally, put the information into frag
    frag%id = k
    frag%layer = layerk
    frag%nat = natf
    allocate (frag%opos(natf),source=0)
    allocate (frag%at(natf),source=0)
    allocate (frag%xyz(3,natf),source=0.0_wp)
    allocate (frag%grd(3,natf),source=0.0_wp)
    n = 0
    do i = 1,nat
      if (taken(i)) then
        n = n+1
        frag%opos(n) = i
        frag%at(n) = at(i)
        frag%xyz(1:3,n) = xyz(1:3,i)
      end if
    end do

    deallocate (taken)
  end subroutine fragment_set_atoms_QMMM

!========================================================================================!
end module lwoniom_setup

