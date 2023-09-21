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
  logical,parameter,private :: debug = .false.

!> routines/datatypes that can be seen outside the module
  public :: lwoniom_data
  public :: lwoniom_initialize

!**************************************************************!
!> this type bundles together most of the
!> data required for a lwONIOM calculation
  type :: lwoniom_data

    !> number of layers
    integer :: nlayer = 0
    logical,allocatable :: layer(:,:)

    !> number of fragments
    integer :: nfrag = 0
    !> list of fragments (subsystems)
    type(structure_data),allocatable :: fragment(:)
    !> root fragment id (the original system)
    integer :: root_id = 0
    !> replace atom-types in any of the layers?
    logical :: replace_at = .false.

    !> further system information
    integer,allocatable :: bond(:,:)
    integer,allocatable :: indexf(:)

    !> calculation mapping
    integer :: ncalcs !> total number of required calculations
    integer,allocatable :: calcids(:,:)  !> high (1,:) and low (2,:) level IDs for each fragment

  contains
    procedure :: info => print_lwoniom_info
    procedure :: deallocate => lwoniom_data_deallocate
    procedure :: add_fragment => lwoniom_add_fragment
    procedure :: dump_fragments => lwoniom_dump_fragments
    procedure :: dump_layers => lwoniom_dump_layers
    procedure :: update => lwoniom_update_fragments
    procedure :: dump_bin => lwoniom_dump_bin
    procedure :: read_bin => lwoniom_read_bin
  end type lwoniom_data
!************************************************************!

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
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
            else if (self%root_id == j) then
              write (myunit,*) '(root)'
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
  &                 layer,fragment,bond,       &
  &                 print,verbose,iunit,iostat)
!******************************************************************
!* This routine SETS UP the lwoniom_data object "dat" and handles
!* the partitioning of the whole system into different layers.
!*
!* It requires input coordinates for the entire system, as
!* well as "layering" information.
!* The layering instruction can be provided by one or two arrays:
!*
!*  A) If only the "layer" array is provided, layers will contain
!*     all the atoms corresponting to the respective layer
!*
!*  B) If the "layer" and the "fragment" arrays are provided,
!*     one layer can have mutliple fragments (MC-ONIOM).
!*
!* The array "bond" can be present to determine connectivity,
!* otherwise this is calculated from covalent radii
!* This routine needs calling only once, at the beginning
!******************************************************************
    character(len=*),parameter :: source = 'lwoniom_initialize'
    !> INPUT
    integer,intent(in) :: nat          !> number of atoms
    integer,intent(in) :: at(nat)      !> atom number for each atom
    real(wp),intent(in) :: xyz(3,nat)  !> Cartesian coordinates for each atom
    !integer,intent(in) :: layer(nat)   !> layer instruction for each atom
    logical,intent(in) :: layer(:,:)   !> layer instruction
    logical,intent(in),optional :: fragment(:,:) !> subsystem instruction (for MC-ONIOM)
    integer,intent(in),optional :: bond(nat,nat) !> bonding/interaction information
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    !> OUTPUT
    integer,intent(out),optional :: iostat
    type(lwoniom_data),intent(inout) :: dat
    !> LOCAL
    type(structure_data) :: tmp
    logical,allocatable :: fragment_tmp(:,:)
    !integer,allocatable :: layer_tmp(:)
    logical,allocatable :: layer_tmp(:,:)
    integer,allocatable :: bond_tmp(:,:)
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun
    integer :: maxf,nroot,nlayer,nfrag
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
! SUBSYTEM SETUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nlayer = size(layer,2)
    if (present(fragment)) then
      nfrag = size(fragment,2)
      if (debug) write (*,*) 'fragment present',nfrag
      allocate (fragment_tmp(nat,nfrag),source=.false.)
      fragment_tmp(:,:) = fragment(:,:)
    else
      nfrag = nlayer
      allocate (fragment_tmp(nat,nfrag),source=.false.)
      fragment_tmp(:,:) = layer(:,:)
    end if
    allocate (layer_tmp(nat,nlayer),source=.false.)
    layer_tmp(:,:) = layer(:,:)
    allocate (indexf(nat))
    do i = 1,nfrag
      if (debug) write (*,*) 'atoms in fragment',i
      if (debug) then
        do j = 1,nat
          if (fragment_tmp(j,i)) write (*,'(1x,i0)',advance='no') j
        end do
        write (*,*)
      end if
    end do

    !> determine maxf and indexf, i.e. the number of "nodes"
    !> and the mapping which atom belongs to which node
    call count_fragments(nat,layer_tmp,fragment_tmp,maxf,indexf)
    dat%nlayer = nlayer
    if (debug) then
      write (*,*) 'maxf',maxf
      write (*,*) 'nlayer',dat%nlayer
      write (*,*) 'indexf',indexf
    end if

    !> For the next steps we need some kind of connectivity information
    !> either provided by the user or determined from covalent radii
    allocate (bond_tmp(nat,nat),source=0)
    if (present(bond)) then
      !> if bond was provided, use those
      bond_tmp = bond
    else
      call lwoniom_rcov_bonds(nat,at,xyz,1.15_wp,bond_tmp)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! from the connectivity determine the tree structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate (parent(maxf))
    if (maxf .eq. nlayer) then
      !> one fragment per layer is easy
      call construct_tree_ONIOM_classic(maxf,parent)

    else if (allocated(bond_tmp)) then
      if (debug) write (*,*) 'creating MC-ONIOMn tree ...'
      call construct_tree_ONIOM_multicenter(nat,layer_tmp,fragment_tmp,bond_tmp,maxf,parent)
      if (debug) write (*,*) 'parent',parent
    else
      write (stderr,'(a)') "**ERROR** 'bond' array not provided in call to "//source
      error stop
    end if
    nroot = count(parent(:) .eq. 0,1)
    if (nroot .ne. 1) then
      write (stderr,'(a,i0)') "**ERROR** incorrect number of parent fragments ",nroot
      error stop
    end if

    !> go through all the fragments, create structure_data and put
    !> them into the fragment list
    do i = 1,maxf
      call tmp%deallocate()
      if (parent(i) < 0) cycle !> safety
      call fragment_set_atoms_ONIOM(nat,at,xyz,layer_tmp,maxf,fragment_tmp,parent,i,tmp)
      call dat%add_fragment(tmp)
    end do
    !> and connect them in the parent-child relations
    do i = 1,maxf
      if (parent(i) .ne. 0) then
        j = parent(i)
        call dat%fragment(j)%add_child(dat%fragment(i))
      else
        dat%root_id = i
      end if
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  linking atoms between layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> at this point, dat%fragment should be set up
!> what's missing are the linking atoms between the layers
! Determine the linking atoms between layers
    call determine_linking_atoms(dat,at,xyz,bond_tmp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally, create storage for true gradients and Jacobians
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1,dat%nfrag
      call dat%fragment(i)%project_gradient_init(nat)
    end do

    if ((io /= 0).and.pr) then
      write (myunit,'("Could not create lwONIOM object ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if

    if (allocated(layer_tmp)) deallocate (layer_tmp)
    if (allocated(fragment_tmp)) deallocate (fragment_tmp)
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
    integer :: k2,m2
    integer,allocatable :: linking_atoms(:,:)

    !> nat is the original system's number of atoms
    nat = size(bond_tmp,1)
    !> linking_atoms is our "work" array. it's second dimension is <nat for each fragment
    allocate (linking_atoms(3,nat))

    !> loop over all fragments
    do i = 1,self%nfrag
      !> loop over all atoms in the fragments i
      !> number of atoms in fragment i inherited from the original structure
      nat_fragment = self%fragment(i)%nat
      if (debug) then
        write (*,*) 'fragment',i,'atoms'
        do j = 1,nat_fragment
          write (*,*) j,'-->',self%fragment(i)%opos(j)
        end do
      end if

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
            !> if k is NOT a member of the fragment, it must be a new link atom
            if (.not.any(self%fragment(i)%opos(:) .eq. k)) then
              if (debug) write (*,'(1x,a,i0,a,i0)') 'link bond   ',l,' - ',k

              !> The only exception is when k is already a link atom for another atom
              if (any(linking_atoms(1,:) .eq. k)) then
                m2 = minloc(abs(linking_atoms(1,:)-k),1)
                k2 = linking_atoms(1,m2)
                linking_atoms(3,m2) = linking_atoms(3,m2)+1
                !write(*,*) k,k2,linking_atoms(1,k2)-k
                write (stdout,'(a,2(i0,a),i0)') '**WARNING** Bad setup - linking atom ',m2, &
                &  ' (atom ',k2,') is connected to multiple atoms in fragment ',i
                cycle
              end if

              !> increment the number of link atoms
              m = m+1
              !> document the reference atom k for the link
              linking_atoms(1,m) = k
              !> document to which atom it is linked to, using the fragment atom index
              linking_atoms(2,m) = j  !> i.e., j not l (!)
              !> count the number of connections for the link atom
              linking_atoms(3,m) = linking_atoms(3,m)+1
            end if
          end if
        end do
      end do
      if (debug) then
        do j = 1,m
          write (*,*) linking_atoms(:,j)
        end do
      end if
      !>allocating m linking atoms for fragment i
      call self%fragment(i)%allocate_link(m)
      call self%fragment(i)%set_link(m,nat,at,xyz,linking_atoms)
    end do

    deallocate (linking_atoms)
  end subroutine determine_linking_atoms

!========================================================================================!

  subroutine lwoniom_data_deallocate(self)
!****************************************************
!* This subroutine deallocates a lwoniom_data object
!*****************************************************
    implicit none
    class(lwoniom_data) :: self
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

  subroutine count_fragments(nat,layer,fragment,maxf,indexf)
!**********************************************************************
!* Given the two tracking arrays "layer" and "fragment",
!* this routine returns the maximum number of fragments (subsets)
!* and enumerates the atoms accordingly.
!**********************************************************************
    implicit none
    character(len=*),parameter :: source = 'lwoniom_count_fragments'
    !> INPUT
    integer,intent(in) :: nat
    logical,intent(inout) :: layer(:,:)
    logical,intent(inout) :: fragment(:,:)
    !> OUTPUT
    integer,intent(out) :: maxf
    integer,intent(out) :: indexf(nat) !> track highest fragment id for each atom
    !> LOCAL
    integer :: i,j,k,l
    integer :: maxl,ninlayer,minl,minf
    integer,allocatable :: tmp(:)
    integer :: ntmp

    !> initialize
    maxf = 0
    indexf(:) = 0
    allocate (tmp(nat))

    !> iterate through the layers and repair indices
    maxf = size(fragment,2)
    maxl = size(layer,2)
    if (debug) write (*,*) 'maxf',maxf,'maxl',maxl
    if (maxl > maxf) then
      write (stderr,'(a)') "**ERROR** Can't have more layers than subsystems!"
      write (stderr,'(10x,a)') 'This can occur if there is exact overlap between two layers.'
      error stop
    end if

    !> classical ONIOM case, one subsystem per layer (set up here)
    !> MC-ONIOM case, just counting stuff
    minf = maxf
    do k = 1,maxf
      j = count(fragment(:,k),1)
      if (j == 0) then
        minf = minf-1
        write (stderr,'(a,i0)') "**ERROR** Empty fragment detected: fragment ",k
      else
        do l = 1,nat
          if (fragment(l,k)) indexf(l) = k
        end do
      end if
    end do
    if (maxf > minf) then
      error stop
    end if

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

  subroutine construct_tree_ONIOM_multicenter(nat,layer,fragment,bond,maxf,parent)
!****************************************************************
!* Set layer hierarchy for MC-ONIOMn schemes, which needs to be
!* done recursively.
!****************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    logical,intent(in) :: layer(:,:)
    logical,intent(in) :: fragment(:,:)
    integer,intent(in) :: bond(nat,nat)
    integer,intent(in) :: maxf
    !> OUTPUT
    integer,intent(out) :: parent(maxf)
    !> LOCAL
    integer :: parentid
    integer :: maxl,thislayer
    integer :: i,j,k,l,nj,nk
    logical :: yes
    parent(:) = -1
    maxl = size(layer,2)
    !> check which of the fragments contains all atoms
    !> there must be exactly ONE such fragment which is our root
    parentid = 0
    do i = 1,maxf
      k = count(fragment(:,i),1)
      if (parentid == 0.and.k == nat) then
        if (debug) write (*,*) 'parent id ',i
        parentid = i
        parent(i) = 0
      end if
      !if(debug) write(*,*) 'atoms in fragment',i
      !if(debug)then
      !   do j=1,nat
      !     if(fragment(j,i)) write(*,'(1x,i0)',advance='no') j
      !   enddo
      !   write(*,*)
      !endif
    end do

    !> go through the fragments recursively (i.e. start with the highest ID)
    jloop: do j = maxf,2,-1
      if (all(.not.fragment(:,j))) then
        !> exclude empty fragments (there shouldn't be any at this point)
        if (debug) write (*,*) 'something is wrong, empty fragment',j
        parent(j) = -1
        cycle
      end if
      !> check all next lower-ID fragments.
      !> the next-lower-ID fragment that contains all the atoms of this fragment
      !> must be the parent
      nj = count(fragment(:,j),1)
      kloop: do k = j-1,1,-1
        yes = .true.
        if (debug) write (*,*) 'comparing fragment',j,'with',k
        nk = count(fragment(:,k),1)
        !> exlude smaller fragments
        if (nk < nj) cycle kloop
        if (all(fragment(:,j).eqv.fragment(:,k))) then
          !> these setups are inefficient, print a warning
          write (stdout,'(a,i0,1x,i0)') '**WARNING** exactly overlapping fragments: ',j,k
        end if
        do l = 1,nat
          if (fragment(l,j)) then
            if (debug) write (*,*) '   atom ',l
            yes = yes.and.fragment(l,k)
          end if
        end do
        if (yes) then
          parent(j) = k
          cycle jloop
        end if
      end do kloop
    end do jloop
    do i = 1,maxf
      if (parent(i) == -1) then
        write (stdout,'(a,i0,a)') '**WARNING** fragment ',i,' not in MC-ONIOM tree!'
      end if
    end do
  end subroutine construct_tree_ONIOM_multicenter

!========================================================================================!
  subroutine fragment_set_atoms_ONIOM(nat,at,xyz,layer,maxf,fragment,parent,k,frag)
!**************************************************************************
!* This subroutine creates a structure_data object frag for the
!* k-th fragment index in an ONIOM setup, based on the present information
!***************************************************************************
    implicit none
    character(len=*),parameter :: source = 'fragment_set_atoms'
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    logical,intent(in) :: layer(:,:)
    integer,intent(in) :: maxf
    logical,intent(in) :: fragment(:,:)
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
      taken(i) = fragment(i,k)
    end do
    layerk = match_layer(k,layer,fragment)

    !> determine how many atoms there are in the fragment total
    natf = count(taken(:),1)
    !> and finally, put the information into frag
    frag%id = k
    frag%layer = layerk
    frag%nat = natf
    allocate (frag%opos(natf),source=0)
    allocate (frag%at(natf),source=0)
    allocate (frag%xyz(3,natf),source=0.0_wp)
    if (frag%layer /= 1) then
      allocate (frag%grd_high(3,natf),source=0.0_wp)
    end if
    allocate (frag%grd_low(3,natf),source=0.0_wp)

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
  end subroutine fragment_set_atoms_ONIOM

!========================================================================================!
  subroutine fragment_set_atoms_QMMM(nat,at,xyz,layer,maxf,fragment,parent,k,frag)
!**************************************************************************
!* This subroutine creates a structure_data object frag for the
!* k-th fragment index in an QMMM setup, based on the present information
!***************************************************************************
    implicit none
    character(len=*),parameter :: source = 'fragment_set_atoms_QMMM'
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    logical,intent(in) :: layer(:,:)
    integer,intent(in) :: maxf
    logical,intent(in) :: fragment(:,:)
    integer,intent(in) :: parent(maxf)
    integer,intent(in) :: k
    type(structure_data),intent(out) :: frag
    logical,allocatable :: taken(:)
    integer :: i,natf,n,layerk,j,l

    call frag%deallocate()
    allocate (taken(nat),source=.false.)
    taken(:) = .false.

    if (k > maxf.or.k < 1) then
      write (stderr,'(a,a)') "**ERROR** invalid fragment index in ",source
      error stop
    end if

    !> add all atoms that belong to the fragment
    do i = 1,nat
      taken(i) = fragment(i,k)
    end do
    layerk = match_layer(k,layer,fragment)

    !> remove all atoms belonging to children (QMMM is additive compared to ONIOM)
    do i = 1,maxf
      if (parent(i) == k) then
        do j = 1,nat
          if (fragment(j,i)) taken(j) = .false.
        end do
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
    if (frag%layer /= 1) then
      allocate (frag%grd_high(3,natf),source=0.0_wp)
    end if
    allocate (frag%grd_low(3,natf),source=0.0_wp)

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

  function match_layer(f,layer,fragment) result(l)
!**************************************************
!* given the layer and fragment mapping,
!* determine which layer fragment f must belong to
!**************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: f !> the fragment-ID
    logical,intent(in) :: layer(:,:) !> the layer-atom mapping
    logical,intent(in) :: fragment(:,:) !> the fragment-atom mapping
    !> OUTPUT
    integer :: l
    !> LOCAL
    integer :: i,j,k
    integer :: maxl,maxf,nat
    logical :: yes
    l = 0 !> could always be the root
    nat = size(layer,1)
    maxf = size(fragment,2)
    maxl = size(layer,2)
    do i = maxl,1,-1 !> go through layers, lowest (innermost) to highest (outermost)
      yes = .true.
      do j = 1,nat
        if (fragment(j,f)) then !> check if all atoms of the fragment belong to layer i
          yes = yes.and.layer(j,i)
        end if
      end do
      if (yes) then !> the first matching is returned
        l = i
        exit
      end if
    end do
  end function match_layer

!========================================================================================!
  subroutine lwoniom_dump_fragments(self,xyz)
!****************************************************
!* Dump all fragments as .xyz
!*****************************************************
    implicit none
    class(lwoniom_data) :: self
    real(wp),intent(in),optional :: xyz(:,:)
    integer :: i

    do i = 1,self%nfrag
      call self%fragment(i)%dump_fragment(xyz)
    end do
  end subroutine lwoniom_dump_fragments


!========================================================================================!
  subroutine lwoniom_dump_layers(self,at,xyz)
!****************************************************
!* Dump all layers as .xyz
!*****************************************************
    implicit none
    class(lwoniom_data) :: self
    integer,intent(in) :: at(:)
    real(wp),intent(in) :: xyz(:,:) !> in Angstroem
    integer :: i,j,k,nat,natlay,l,ich
    logical,allocatable :: taken(:)
    character(len=40) :: atmp
    nat = size(xyz,2)
    allocate(taken(nat))
    do i=1,self%nlayer
       taken = .false.
       do j=1,self%nfrag
         if(self%fragment(j)%layer == i)then
           do k=1,self%fragment(j)%nat
             l = self%fragment(j)%opos(k)
             taken(l) = .true.
           enddo
         endif
       enddo
       natlay = count(taken(:),1)
       write(atmp,'(i0)') i 
       open(newunit=ich,file='layer.'//trim(atmp)//'.xyz')
       write(ich,*) natlay
       write(ich,*)
       do j=1,nat
         if(taken(j))then
           write(ich,'(a4,3f25.15)') PSE(at(j)),xyz(1:3,j)
         endif
       enddo
       close(ich)
    enddo
    deallocate(taken)
  end subroutine lwoniom_dump_layers


!========================================================================================!

  subroutine lwoniom_update_fragments(self,xyz,pc)
!*********************************
!* The geometry on all fragments.
!* Optionally, also point charges
!*********************************
    implicit none
    class(lwoniom_data) :: self
    real(wp),intent(in) :: xyz(:,:)
    real(wp),intent(in),optional :: pc(:)
    integer :: i,nat
    nat = size(xyz,2)
    do i = 1,self%nfrag
      call self%fragment(i)%update(nat,xyz,pc)
    end do
  end subroutine lwoniom_update_fragments

!========================================================================================!

  subroutine lwoniom_dump_bin(self)
!***********************************************************
!* Write routine for an lwoniom_data object to binary file
!* Makes easier I/O accesible
!***********************************************************
    implicit none
    class(lwoniom_data) :: self
    integer :: bin
    integer :: d1,d2,i,j
    logical :: bdum

    open (newunit=bin,file='lwoniom.data',status='replace',access='stream')!form='binary')

    !> number of layers
    !integer :: nlayer = 0
    write (bin) self%nlayer

    !integer,allocatable :: layer(:)
    bdum = allocated(self%layer)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%layer,1)
      d2 = size(self%layer,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%layer(i,j)
        end do
      end do
    end if

    !> number of fragments
    !integer :: nfrag = 0
    write (bin) self%nfrag

    !> list of fragments (subsystems)
    !type(structure_data),allocatable :: fragment(:)
    bdum = allocated(self%fragment)
    write (bin) bdum
    if (bdum) then
      d1 = self%nfrag
      write (bin) d1
      !write(*,*) 'dump fragments', d1
      do i = 1,d1
        call self%fragment(i)%dump_bin(bin)
      end do
    end if

    !> root fragment id (the original system)
    !integer :: root_id = 0
    write (bin) self%root_id

    !> replace atom-types in any of the layers?
    !logical :: replace_at = .false.
    write (bin) self%replace_at

    !> further system information
    !integer,allocatable :: bond(:,:)
    bdum = allocated(self%bond)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%bond,1)
      d2 = size(self%bond,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%bond(i,j)
        end do
      end do
    end if

    !integer,allocatable :: indexf(:)
    bdum = allocated(self%indexf)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%indexf,1)
      write (bin) d1
      do i = 1,d1
        write (bin) self%indexf(i)
      end do
    end if

    !> calculation mapping
    !integer :: ncalcs !> total number of required calculations
    write (bin) self%ncalcs

    !integer,allocatable :: calcids(:,:)  !> high (1,:) and low (2,:) level IDs for each fragment
    bdum = allocated(self%calcids)
    write (bin) bdum
    if (bdum) then
      d1 = size(self%calcids,1)
      d2 = size(self%calcids,2)
      write (bin) d1
      write (bin) d2
      do i = 1,d1
        do j = 1,d2
          write (bin) self%calcids(i,j)
        end do
      end do
    end if

    close (bin)
  end subroutine lwoniom_dump_bin

  subroutine lwoniom_read_bin(self,success)
!***********************************************************
!* Read routine for an lwoniom_data object to binary file
!* Makes easier I/O accesible
!***********************************************************
    implicit none
    class(lwoniom_data) :: self
    integer :: bin
    integer :: d1,d2,i,j
    logical :: bdum,ex
    logical,intent(out) :: success

    success = .false.
    inquire (file='lwoniom.data',exist=ex)
    if (.not.ex) return

    open (newunit=bin,file='lwoniom.data',status='old',access='stream')!,form='binary')

    !> number of layers
    !integer :: nlayer = 0
    read (bin) self%nlayer

    !integer,allocatable :: layer(:)
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      read (bin) d2
      allocate (self%layer(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%layer(i,j)
        end do
      end do
    end if

    !> number of fragments
    !integer :: nfrag = 0
    read (bin) self%nfrag

    !> list of fragments (subsystems)
    !type(structure_data),allocatable :: fragment(:)
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%fragment(d1))
      !write(*,*) 'read fragments',d1
      do i = 1,d1
        call self%fragment(i)%read_bin(bin)
        call self%fragment(i)%dump_fragment()
      end do
    end if

    !> root fragment id (the original system)
    !integer :: root_id = 0
    read (bin) self%root_id

    !> replace atom-types in any of the layers?
    !logical :: replace_at = .false.
    read (bin) self%replace_at

    !> further system information
    !integer,allocatable :: bond(:,:)
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      read (bin) d2
      allocate (self%bond(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%bond(i,j)
        end do
      end do
    end if

    !integer,allocatable :: indexf(:)
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      allocate (self%indexf(d1))
      do i = 1,d1
        read (bin) self%indexf(i)
      end do
    end if

    !> calculation mapping
    !integer :: ncalcs !> total number of required calculations
    read (bin) self%ncalcs

    !integer,allocatable :: calcids(:,:)  !> high (1,:) and low (2,:) level IDs for each fragment
    read (bin) bdum
    if (bdum) then
      read (bin) d1
      read (bin) d2
      allocate (self%calcids(d1,d2))
      do i = 1,d1
        do j = 1,d2
          read (bin) self%calcids(i,j)
        end do
      end do
    end if

    close (bin)
    success = .true.

  end subroutine lwoniom_read_bin

!========================================================================================!
!========================================================================================!
end module lwoniom_setup

