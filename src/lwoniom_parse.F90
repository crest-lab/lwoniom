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
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
#ifdef WITH_TOMLF
  use tomlf
#endif
  implicit none
  private

!> exports
  public :: lwoniom_parse_inputfile

!> storage type for input
  public :: lwoniom_input
  type :: lwoniom_input
    character(len=:),allocatable :: structurefile
    integer :: nat = 0
    integer :: maxfragments = 0
    integer :: maxlayers = 0
    integer,allocatable :: layer(:)
    integer,allocatable :: frag(:)
    character(len=:),allocatable :: wbofile
  end type lwoniom_input

!> printout param
   character(len=*),parameter,private :: ns = 'lwONIOM> '
   character(len=*),parameter,private :: clears = '         '

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_parse_inputfile(tomlfile,input)
    implicit none
    character(len=*),intent(in)     :: tomlfile
    type(lwoniom_input),intent(out) :: input
#ifndef WITH_TOMLF
    write (stderr,*) "**ERROR** lwONIOM was compiled without TOML-f support. "// &
    & "To enable the parser, set up the build with -DWITH_TOMLF=true"
    error stop
#else /* WITH_TOMLF */
    integer :: io
    type(toml_table),allocatable       :: table
    type(toml_error),allocatable       :: error
    type(toml_table),pointer :: child
    logical :: ex

    inquire (file=tomlfile,exist=ex)
    if (.not.ex) then
      write (stderr,'("**ERROR** ",a,a,a)') "input file ",trim(tomlfile)," not found!"
      return
    end if
    write (stdout,'(a,a,a)') ns,"reading input file ",trim(tomlfile)

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
      write (stderr,'("**ERROR** ",a,a)') "No [lwoniom] section found in input file ",trim(tomlfile)
      return
    end if
    call read_lwoniom_block(error,input,child)

    if (allocated(error)) deallocate (error)
    if (allocated(table)) deallocate (table)
#endif
  end subroutine lwoniom_parse_inputfile

!========================================================================================!

#ifdef WITH_TOMLF
  subroutine read_lwoniom_block(error,input,table)
!***************************************************
!* Reader for [lwoniom] block in toml file
!***************************************************
    !> Error handler
    type(toml_error),allocatable :: error
    !> Hamiltonian input to be read
    type(lwoniom_input),intent(out) :: input
    !> Data structure
    type(toml_table),intent(inout) :: table
    type(toml_table),pointer :: child
    integer :: io,i,j,k,l
    integer :: ikey
    type(toml_key),allocatable :: list(:)

    !> first get total number of atoms and allocate
    call get_value(table,"natoms",input%nat,stat=io)
    if (io /= toml_stat%success) then
      write (stderr,'("**ERROR** ",a)') "Please provide the total number of atoms (natoms) in [lwoniom] block"
      return
    end if
    allocate (input%layer(input%nat),source=0)
    allocate (input%frag(input%nat),source=0)

    !> atom-wise definitions of fragments first
    call get_value(table,'fragment',child,requested=.false.)
    if (associated(child)) then
      call read_lwoniom_fragment(error,input,child)
    end if

    !> then fragment-wise definition of layers (= which fragment belongs to which layer)
    !> note: input%layer is still stored atom-wise!
    call get_value(table,'layer',child,requested=.false.)
    if (associated(child)) then
      if (input%maxfragments < 1) then
        write (stderr,'("**ERROR** ",a)') 'Layers must not be defined without defining fragments first'
        stop
      end if
      call read_lwoniom_layer(error,input,child)
    end if

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
    integer,allocatable :: frag(:)

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

          !> otherwise, some string shortcuts can be defined
        else
          call get_value(table,key,val,stat=io3)
          if (io3 == 0) then
            select case (val)
            case ("all")
              allocate (frag(input%nat))
              do i = 1,input%nat
                frag(i) = i
              enddo
              call set_fragment(input,f,frag)
            end select
          end if
        end if
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
      end if
    end do

  end subroutine read_lwoniom_layer

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
    if (f > input%maxfragments) input%maxfragments = f
    k = size(frag(:),1)
    write (stdout,'(a,a,i0)') ns,'setup fragment ',f
    c = 0
    do i = 1,k
      j = frag(i)
      if (input%frag(j) > f) then
        write (stdout,'(a,i0)') '**WARNING** atom ',j,' already defined in fragment ',input%frag(j)
      else
        c = c + 1
        input%frag(j) = f
      end if
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
    if (l > input%maxlayers) input%maxlayers = l
    k = size(lay(:),1)
    do i = 1,k
      j = lay(i)
      if (.not.any(input%frag(:) .eq. j)) then
        write (stdout,'(a,i0,a,i0)') '**WARNING** no atoms associated with fragment ',j, &
        & ' while trying to add layer ',l
      else
        do m = 1,input%nat
          if (input%frag(m) == j) then
            input%layer(m) = l
          end if
        end do
        write (stdout,'(a,2(a,i0))') ns,'added fragment ',j,'  -->  layer ',l
      end if
    end do
  end subroutine set_layer

!========================================================================================!
!========================================================================================!
end module lwoniom_parse

