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

module lwoniom_parse
!********************************************
!* This module implements the toml-f parser
!* to provide an input parser for lwONIOM
!********************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
#ifdef WITH_TOMLF
  use tomlf
#endif
  use lwoniom_structures,only:zsym_to_at,lowercase,PSE
  implicit none
  private

!> exports
  public :: lwoniom_parse_inputfile

!> storage type for input
  public :: lwoniom_input
  type :: lwoniom_input
    logical :: read_input = .false.
    character(len=:),allocatable :: structurefile
    integer :: nat = 0
    integer :: maxfragments = 0
    integer :: maxlayers = 0
    logical,allocatable :: layer(:,:)
    logical,allocatable :: frag(:,:)
    character(len=:),allocatable :: wbofile
    real(wp),allocatable :: wbo(:,:)
    character(len=2),allocatable :: zsym(:)
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    character(len=:),allocatable :: cmd(:)
    integer,allocatable :: layerlvl(:)
    integer,allocatable :: fragmentchrg(:)
    integer,allocatable :: layerreplace(:,:)
    logical :: dump_frag = .false.
    logical :: dump_layer = .false.
    logical :: try_bin = .false.
  contains
    procedure :: deallocate => deallocate_lwoniom_input
    procedure :: parse_xyz => parse_structure_xyz
  end type lwoniom_input

!> printout param
  character(len=*),parameter,private :: ns = 'lwONIOM> '
  character(len=*),parameter,private :: clears = '         '

  logical,parameter :: debug = .false. 

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine deallocate_lwoniom_input(self)
    implicit none
    class(lwoniom_input) :: self
    if (allocated(self%structurefile)) deallocate (self%structurefile)
    self%nat = 0
    self%maxfragments = 0
    self%maxlayers = 0
    if (allocated(self%layer)) deallocate (self%layer)
    if (allocated(self%frag)) deallocate (self%frag)
    if (allocated(self%wbofile)) deallocate (self%wbofile)
    if (allocated(self%wbo)) deallocate (self%wbo)
    if (allocated(self%zsym)) deallocate (self%zsym)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%cmd)) deallocate (self%cmd)
    if (allocated(self%layerlvl)) deallocate (self%layerlvl)
    if (allocated(self%layerreplace)) deallocate (self%layerreplace)
  end subroutine deallocate_lwoniom_input

!=======================================================================================!

  subroutine lwoniom_parse_inputfile(tomlfile,input,required,natoms)
    implicit none
    character(len=*),intent(in)     :: tomlfile
    type(lwoniom_input),intent(out) :: input
    logical,intent(in),optional :: required
    integer,intent(in),optional :: natoms
#ifndef WITH_TOMLF
    write (stderr,*) "**ERROR** lwONIOM was compiled without TOML-f support. "// &
    & "To enable the parser, set up the build with -DWITH_TOMLF=true"
    error stop
#else /* WITH_TOMLF */
    integer :: io
    type(toml_table),allocatable       :: table
    type(toml_error),allocatable       :: error
    type(toml_table),pointer :: child
    logical :: ex,req

    req = .true.
    if (present(required)) then
      req = required
    end if

    inquire (file=tomlfile,exist=ex)
    if (.not.ex) then
      write (stderr,'("**ERROR** ",a,a,a)') "input file ",trim(tomlfile)," not found!"
      return
    end if
    if (req) write (stdout,'(a,a,a)') ns,"reading input file ",trim(tomlfile)

    !> the actual "reading" part
    open (newunit=io,file=tomlfile)
    call toml_parse(table,io,error)
    close (unit=io)
    if (allocated(error)) then
      print*,"Error parsing table:"//error%message
    end if

    !> look for [lwoniom] block
    call get_value(table,"lwoniom",child,requested=.false.)
    if (.not.associated(child)) then
      if (req) write (stderr,'("**ERROR** ",a,a)') "No [lwoniom] section found in input file ",trim(tomlfile)
      if (allocated(error)) deallocate (error)
      if (allocated(table)) deallocate (table)
      return
    end if
    input%read_input = .true.
    write (stdout,'(/,a,a,a)') ns,'parsing [lwoniom]-block from file ',trim(tomlfile)
    call read_lwoniom_block(error,input,child,natoms)

    if (allocated(error)) deallocate (error)
    if (allocated(table)) deallocate (table)
#endif
  end subroutine lwoniom_parse_inputfile

!========================================================================================!

#ifdef WITH_TOMLF
  subroutine read_lwoniom_block(error,input,table,natoms)
!***************************************************
!* Reader for [lwoniom] block in toml file
!***************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(out) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    integer,intent(in),optional :: natoms
    type(toml_table),pointer :: child
    integer :: io,i,j,k,l,tmpnat
    integer :: ikey
    type(toml_key),allocatable :: list(:)
    character(len=:),allocatable :: key
    character(len=:),allocatable :: val
    logical :: ex

    !> first get total number of atoms and allocate
    call get_value(table,"natoms",input%nat,stat=io)
    if (io /= toml_stat%success.and..not.present(natoms)) then
      write (stderr,'("**ERROR** ",a)') "Please provide the total number of atoms (natoms) in [lwoniom] block"
      return
    else if (present(natoms)) then
      input%nat = natoms
    end if

    !> iterate over keys (to check for input xyz, and other info)
    call table%get_keys(list)
    do ikey = 1,size(list)
      key = list(ikey)%key
      select case (key)
      case ('structure','struct','input','xyz','struc')
        !> try to parse an xyz file
        call get_value(table,key,input%structurefile,stat=io)
        if (io == toml_stat%success) then
          write (stdout,'(a,a,a)') ns,'reading XYZ from ',trim(input%structurefile)
          call read_xyz(input%structurefile,tmpnat,input%zsym,input%xyz)
          if (tmpnat /= input%nat) then
            write (stderr,'("**ERROR** ",3a)') "mismatch of number of atoms in [lwoniom] block ", &
            & 'and ',input%structurefile
            return
          end if
        end if

      case ('wbo','bo','topo')
        !> try to read a file containing bond info
        call get_value(table,key,input%wbofile,stat=io)
        if (io == toml_stat%success) then
          write (stdout,'(a,a,a)') ns,'reading bonds from ',trim(input%wbofile)
          call read_bo(input%wbofile,input%nat,input%wbo)
        end if

      case ('restart')
        call get_value(table,key,input%try_bin,stat=io)
        if (io == toml_stat%success) then
          if (input%try_bin) write (stdout,'(a,a)') ns,'trying to restart from from lwoniom.data'
        else
          input%try_bin = .false.
        end if

      end select
    end do

    !> check if we are in restart mode and potentially can restart
    if (input%try_bin) then
      inquire (file='lwoniom.data',exist=ex)
      if (ex) return !> we can leave the routine in this case. good luck.
    end if

    !> atom-wise definitions of fragments first
    call get_value(table,'fragment',child,requested=.false.)
    if (associated(child)) then
      call read_lwoniom_fragment(error,input,child)
    end if

    !> lwONIOM MUST have fragment information
    if (input%maxfragments < 1) then
      write (stderr,'("**ERROR** ",a)') 'lwONIOM requires definition of fragments!'
      error stop
    end if

    !> then fragment-wise definition of layers (= which fragment belongs to which layer)
    !> note: input%layer is still stored atom-wise!
    call get_value(table,'layer',child,requested=.false.)
    if (associated(child)) then
      call read_lwoniom_layer(error,input,child)
    else

      !> If NOT present, each fragment is assumed to be a layer on its own (FALLBACK)
      write (stdout,'(a,a)') ns,'No layer information provided; assuming fragments to be layers'
      input%maxlayers = input%maxfragments
      if (.not.allocated(input%layer)) allocate (input%layer(input%nat,input%maxlayers))
      input%layer(:,:) = .false.
      do i = input%maxfragments,1,-1 !> going from high-index to low-index!
        write (stdout,'(a,a)') clears,'fragment ',i,' --> layer ',i
        do j = 1,input%nat
          if (input%frag(j,i)) input%layer(j,i) = .true.
          !> in the fallback we are assuming hierarchical layer numbering!
          !> I.e., all higher-index layer (i-1) members are also assumed to be members of layer i
          if (i < input%maxfragments) then
            input%layer(j,i) = input%layer(j,i).or.input%layer(j,i+1)
          end if
        end do
      end do

    end if

    !> then fragment-wise definition of subprocesses to calculate energies and gradients
    !> note: this is only for the app usage
    call get_value(table,'cmd',child,requested=.false.)
    if (associated(child)) then
      if (input%maxlayers < 1) then
        write (stderr,'("**ERROR** ",a)') 'CMDs must not be defined without defining layers first'
        stop
      end if
      allocate (input%cmd(input%maxlayers),source=repeat(" ",200))
      call read_lwoniom_processcmd(error,input,child)
    end if

    !> then layer-wise definition of calculator IDs
    call get_value(table,'layerlevel',child,requested=.false.)
    if (associated(child)) then
      if (input%maxlayers < 1) then
        write (stderr,'("**ERROR** ",a)') 'Layer level IDs must not be defined without defining layers first'
        stop
      end if
      call read_lwoniom_layerlvl(error,input,child)
    end if

    !> replace elements in some layers
    call get_value(table,'replace',child,requested=.false.)
    if (associated(child)) then
      if (input%maxlayers < 1) then
        write (stderr,'("**ERROR** ",a)') 'Layer replacements must not be defined without defining layers first'
        stop
      end if
      call read_lwoniom_layerreplace(error,input,child)
    end if


    !> read charges for fragments
    call get_value(table,'chrg',child,requested=.false.)
    if (associated(child)) then
      call read_lwoniom_fragcharge(error,input,child)
    end if

    if(debug) stop
  end subroutine read_lwoniom_block

!========================================================================================!

  subroutine read_lwoniom_fragment(error,input,table)
!*********************************
!* Read all fragments (atom-wise)
!*********************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(inout) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    integer :: io,i,j,k,l,f,io2,io3
    integer :: ikey
    type(toml_key),allocatable :: list(:)
    type(toml_array),pointer    :: arr
    character(len=:),allocatable :: key
    character(len=:),allocatable :: val
    logical :: valbool
    integer,allocatable :: frag(:)
    logical,allocatable :: atlist(:)

    !> iterate over keys (which should be integer numbers)
    call table%get_keys(list)
    do ikey = 1,size(list)
      key = list(ikey)%key
      read (key,*,iostat=io) f
      if (io == 0) then
        if (allocated(frag)) deallocate (frag)
        call get_value(table,key,arr,stat=io2)

        !> if given as a list of integers, these are the corresponding atoms
        if (associated(arr)) then
          k = len(arr)
          allocate (frag(k))
          do l = 1,k
            call get_value(arr,l,frag(l))
          end do
          call set_fragment(input,f,frag)

        else
          !> otherwise, some string shortcuts can be defined

          call get_value(table,key,val,stat=io3)
          if (io3 == 0) then
            select case (val)
            case ("all")
              allocate (frag(input%nat))
              do i = 1,input%nat
                frag(i) = i
              end do
              call set_fragment(input,f,frag)
            case default
              !> atom list for defining ranges

              call lwoniom_get_atlist(input%nat,atlist,val)
              k = count(atlist,1)
              allocate (frag(k))
              k = 0
              do i = 1,input%nat
                if (atlist(i)) then
                  k = k+1
                  frag(k) = i
                end if
              end do
              call set_fragment(input,f,frag)

            end select
          end if
        end if
      else !> if the key is not a number, it might be accepted string
        select case( key )
        case ( 'dump' )
          call get_value(table,key,valbool,stat=io3)
          if(io3 == 0) input%dump_frag = valbool
        end select
      end if
    end do

  end subroutine read_lwoniom_fragment

!========================================================================================!

  subroutine read_lwoniom_layer(error,input,table)
!******************************************************************
!* Read all layers
!* (fragment-wise, i.e., different fragments may belong to
!*  the same layer --> MC-ONIOM)
!* Note, layer information is still saved atom-wise to input%layer
!******************************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(inout) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    integer :: io,i,j,k,l,f,io2,io3
    integer :: ikey
    type(toml_key),allocatable :: list(:)
    type(toml_array),pointer    :: arr
    character(len=:),allocatable :: key
    character(len=:),allocatable :: val
    logical :: valbool
    integer,allocatable :: lay(:)

    !> iterate over keys (which should be integer numbers)
    call table%get_keys(list)
    do ikey = 1,size(list)
      key = list(ikey)%key
      read (key,*,iostat=io) l
      if (io == 0) then
        if (allocated(lay)) deallocate (lay)
        call get_value(table,key,arr,stat=io2)
        !> if given as a list of integers, these are the corresponding fragments
        if (associated(arr)) then
          k = len(arr)
          allocate (lay(k))
          do j = 1,k
            call get_value(arr,j,lay(j))
          end do
          call set_layer(input,l,lay)
        end if
      else !> if the key is not a number, it might be accepted string
        select case( key )
        case ( 'dump' )
          call get_value(table,key,valbool,stat=io3)
          if(io3 == 0) input%dump_layer = valbool
        end select
      end if
    end do

  end subroutine read_lwoniom_layer

  subroutine read_lwoniom_processcmd(error,input,table)
!*******************************************************
!* Read all process commands
!* Note, there must be one process command per layer as
!* each layer corresponds to one level of theory
!*******************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(inout) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    integer :: io,i,j,k,l,f,io2,io3
    integer :: ikey
    type(toml_key),allocatable :: list(:)
    type(toml_array),pointer    :: arr
    character(len=:),allocatable :: key
    character(len=:),allocatable :: val
    integer,allocatable :: lay(:)

    !> iterate over keys (which should be integer numbers)
    call table%get_keys(list)
    if (size(list) .ne. input%maxlayers) then
      write (stderr,'("**ERROR** ",a)') 'Please define a CMD for all layers!'
      error stop
    end if
    do ikey = 1,size(list)
      key = list(ikey)%key
      read (key,*,iostat=io) l
      if (io == 0) then
        if (allocated(lay)) deallocate (lay)
        call get_value(table,key,val,stat=io2)
        !> if a command was specified, save it
        if (io2 == 0) then
          input%cmd(l) = val
          write (stdout,'(a,a,i0,a)') ns,'layer ',l,' subprocess command set to:'
          write (stdout,'(a,a )') repeat(' ',len(ns)),trim(input%cmd(l))
        end if
      end if
    end do

  end subroutine read_lwoniom_processcmd

  subroutine read_lwoniom_layerlvl(error,input,table)
!*******************************************************
!* Read all an id for the level of theory that
!* Shall be used for each of the layers
!*******************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(inout) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    integer :: io,i,j,k,l,f,io2,io3
    integer :: ikey
    type(toml_key),allocatable :: list(:)
    type(toml_array),pointer    :: arr
    character(len=:),allocatable :: key
    integer :: val
    integer,allocatable :: lay(:)
    integer :: lvl

    !> iterate over keys (which should be integer numbers)
    call table%get_keys(list)
    if (size(list) .ne. input%maxlayers) then
      write (*,*) size(list),input%maxlayers
      write (stderr,'("**ERROR** ",a)') 'Please define a level ID for ALL layers!'
      error stop
    end if
    if (.not.allocated(input%layerlvl)) allocate (input%layerlvl(input%maxlayers))
    do ikey = 1,size(list)
      key = list(ikey)%key
      read (key,*,iostat=io) l
      if (io == 0) then
        if (allocated(lay)) deallocate (lay)
        call get_value(table,key,val,stat=io2)
        if (io2 == 0) then
          input%layerlvl(l) = val
          write (stdout,'(a,a,i0,a,i0)') ns,'layer ',l,' associated with calculation level ',val
        end if
      end if
    end do

  end subroutine read_lwoniom_layerlvl

  subroutine read_lwoniom_layerreplace(error,input,table)
!*******************************************************
!* Read elements to replace in a given layer
!* Shall be used for each of the layers
!* the toml syntax will look like this:
!* replace.1.he = 'xe'   #to replace He with Xe
!*******************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(inout) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    type(toml_table) :: newtable
    integer :: io,i,j,k,l,f,io2,io3,i1,i2
    integer :: ikey,childikey
    type(toml_key),allocatable :: list(:)
    type(toml_key),allocatable :: childlist(:)
    type(toml_array),pointer    :: arr
    character(len=:),allocatable :: key
    character(len=:),allocatable :: val
    integer,allocatable :: lay(:)
    integer :: lvl

    !> iterate over keys (which should be integer numbers)
    call table%get_keys(list)
    if (.not.allocated(input%layerreplace)) &
    &  allocate (input%layerreplace(118,input%maxlayers),source=0)
    do ikey = 1,size(list)
      key = list(ikey)%key
      read (key,*,iostat=io) l !> the key must be a number referring to a layer
      if (io == 0.and.l <= input%maxlayers) then

        call get_value(table,key,child,requested=.false.)
        if (associated(child)) then
          call child%get_keys(childlist) !> the other keys are element symbols
          do childikey = 1,size(childlist)
            key = childlist(childikey)%key
            call get_value(child,key,val,stat=io2)
            if (io2 == 0) then
              i1 = zsym_to_at(key)
              i2 = zsym_to_at(val)
              input%layerreplace(i1,l) = i2
            end if
          end do
        end if
      end if
    end do

    k = 0
    do l = 1,input%maxlayers
      do i = 1,118
        j = input%layerreplace(i,l)
        if (j > 0) then
          k = k+1
          write (stdout,'(4a,i0,2a)') ns,'replacing "',trim(PSE(i)),'" atoms in layer ',l, &
          & ' with "',trim(PSE(j))//'"'
        end if
      end do
    end do
    if (k < 1) deallocate (input%layerreplace)

  end subroutine read_lwoniom_layerreplace


  subroutine read_lwoniom_fragcharge(error,input,table)
!*******************************************************
!* Read a charge for a given fragment
!* the toml syntax will look like this:
!* chrg.1 = 2   #to set charge of fragment 1 to 2
!*******************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(inout) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    type(toml_table) :: newtable
    integer :: io,i,j,k,l,f,io2,io3,i1,i2
    integer :: ikey,childikey
    type(toml_key),allocatable :: list(:)
    type(toml_key),allocatable :: childlist(:)
    type(toml_array),pointer    :: arr
    character(len=:),allocatable :: key
    integer :: val
    integer,allocatable :: lay(:)

    !> iterate over keys (which should be integer numbers)
    call table%get_keys(list)
    do ikey = 1,size(list)
      key = list(ikey)%key
      read (key,*,iostat=io) l !> the key must be a number referring to a layer
      if (io == 0.and.l <= input%maxfragments) then
        call get_value(table,key,val,stat=io2) !> read the charge
        if (io2 == 0) then
          if(.not.allocated(input%fragmentchrg))then
             allocate(input%fragmentchrg(input%maxfragments), source=-9990)
          endif
          input%fragmentchrg(l) = val
          write (stdout,'(a,i0,a,i0)') ns//'setting charge of fragment ',l,' to ',val
        endif
      end if
    end do
  end subroutine read_lwoniom_fragcharge

#endif
!========================================================================================!
!========================================================================================!

  subroutine set_fragment(input,f,frag)
!**********************************
!* Add all the atoms in "frag" to
!* fragment f in the input storage
!***********************************
    implicit none
    type(lwoniom_input),intent(inout) :: input
    integer,intent(in) :: f
    integer,intent(in) :: frag(:)
    integer :: i,j,k,c
    logical,allocatable :: frag_tmp(:,:)
    if (f > input%maxfragments) then
      allocate (frag_tmp(input%nat,f),source=.false.)
      do i = 1,input%maxfragments
        frag_tmp(:,i) = input%frag(:,i)
      end do
      input%maxfragments = f
      call move_alloc(frag_tmp,input%frag)
    end if
    k = size(frag(:),1)
    write (stdout,'(a,a,i0)') ns,'setup fragment ',f
    c = 0
    do i = 1,k
      j = frag(i)
      c = c+1
      input%frag(j,f) = .true.
    end do
    write (stdout,'(a,i0,a)') clears,c,' atoms were selected'
  end subroutine set_fragment

!========================================================================================!

  subroutine set_layer(input,l,lay)
!******************************************
!* Add all the atoms belonging to fragments
!* defined in "lay" to layer l
!*******************************************
    implicit none
    type(lwoniom_input),intent(inout) :: input
    integer,intent(in) :: l
    integer,intent(in) :: lay(:)
    integer :: i,j,k,m
    logical,allocatable :: layer_tmp(:,:)
    if (l > input%maxlayers) then
      allocate (layer_tmp(input%nat,l),source=.false.)
      do i = 1,input%maxlayers
        layer_tmp(:,i) = input%layer(:,i)
      end do
      input%maxlayers = l
      call move_alloc(layer_tmp,input%layer)
    end if
    k = size(lay(:),1)
    do i = 1,k
      j = lay(i)
      if (.not.any(input%frag(:,j))) then
        write (stdout,'(a,i0,a,i0)') '**WARNING** no atoms associated with fragment ',j, &
        & ' while trying to add layer ',l
      else
        do m = 1,input%nat
          if (input%frag(m,j)) then
            input%layer(m,l) = .true.
          end if
        end do
        write (stdout,'(a,2(a,i0))') ns,'added fragment ',j,'  -->  layer ',l
      end if
    end do
  end subroutine set_layer

!========================================================================================!

  subroutine read_xyz(fname,nat,zsym,xyz)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(out) :: nat
    character(len=2),intent(out),allocatable :: zsym(:)
    real(wp),intent(out),allocatable :: xyz(:,:)
    integer :: ich,i,j,k,l
    logical :: ex
    character(len=256) :: atmp
    inquire (file=fname,exist=ex)
    if (ex) then
      open (newunit=ich,file=fname)
      read (ich,*) nat
      allocate (zsym(nat))
      allocate (xyz(3,nat),source=0.0_wp)
      read (ich,'(a)') atmp
      do i = 1,nat
        read (ich,*) zsym(i),xyz(1:3,i)
      end do
      close (ich)
    else
      write (stderr,'(a)') '**ERROR** Input file '//trim(fname)//' does not exist'
    end if
  end subroutine read_xyz

!========================================================================================!

  subroutine parse_structure_xyz(input,fname)
    implicit none
    class(lwoniom_input) :: input
    character(len=*),intent(in) :: fname
    integer :: tmpnat
    input%structurefile = fname
    call read_xyz(input%structurefile,tmpnat,input%zsym,input%xyz)
    if (input%nat == 0) input%nat = tmpnat
  end subroutine parse_structure_xyz

!========================================================================================!

  subroutine read_bo(fname,nat,bo)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    real(wp),intent(out),allocatable :: bo(:,:)
    integer :: ich,i,j,k,l,io
    logical :: ex
    character(len=256) :: atmp
    real(wp) :: dum
    inquire (file=fname,exist=ex)
    if (ex) then
      open (newunit=ich,file=fname)
      allocate (bo(nat,nat))
      do
        read (ich,'(a)',iostat=io) atmp
        if (io /= 0) exit
        read (atmp,*) i,j,dum
        bo(i,j) = dum
        bo(j,i) = dum
      end do
      close (ich)
    end if
  end subroutine read_bo

!========================================================================================!

  subroutine lwoniom_get_atlist(nat,atlist,line,at)
!******************************************************
!* Analyze a string containing atom specifications.
!* "atlist" is a array of booleans for each atom,
!* which is set to .true. should the atom be contained
!* in atlist.
!******************************************************
    implicit none
    integer,intent(in) :: nat
    logical,intent(out),allocatable :: atlist(:)
    character(len=*),intent(in) :: line
    integer,intent(in),optional :: at(nat)
    character(len=:),allocatable :: substr(:)
    integer :: i,j,k,l,io,ns,ll,i1,i2,io1,io2,i3,i4
    character(len=:),allocatable :: atmp,btmp

    allocate (atlist(nat),source=.false.)
!>-- count stuff
    ll = len_trim(line)
    ns = 1
    do i = 1,ll
      if (line(i:i) .eq. ',') ns = ns+1
    end do
    allocate (substr(ns),source=repeat(' ',ll))
!>-- cut stuff
    if (ns > 1) then
      j = 1
      k = 1
      do i = 1,ll
        if (k == ns) then
          substr(k) = lowercase(adjustl(line(j:)))
          exit
        end if
        if (line(i:i) .eq. ',') then
          substr(k) = lowercase(adjustl(line(j:i-1)))
          k = k+1
          j = i+1
        end if
      end do
    else
      substr(1) = trim(line)
    end if
!>--- analyze stuff
    do i = 1,ns
      atmp = trim(substr(i))
      if (atmp .eq. 'all') then
        atlist(:) = .true.
        exit
      end if
      if (index(atmp,'.') .ne. 0) cycle !> exclude floats
      l = index(atmp,'-')
      if (l .eq. 0) then
        !> single atom
        read (atmp,*,iostat=io) i1
        if (io /= 0) then
          i2 = zsym_to_at(atmp)
          if (i2 /= 0.and.present(at)) then

          end if
        else
          atlist(i1) = .true.
        end if
      else
        !> range of atoms
        btmp = atmp(:l-1)
        read (btmp,*,iostat=io1) i1
        btmp = atmp(l+1:)
        read (btmp,*,iostat=io2) i2
        if (io1 .eq. 0.and.io2 .eq. 0) then
          i4 = max(i1,i2)
          i3 = min(i1,i2)
          do j = 1,nat
            if (i3 <= j.and.j <= i4) atlist(j) = .true.
          end do
        end if
      end if
    end do
    deallocate (substr)
  end subroutine lwoniom_get_atlist

!========================================================================================!
!========================================================================================!
end module lwoniom_parse

