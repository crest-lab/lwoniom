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
    integer :: maxlayer = 0
    integer,allocatable :: layer(:)
    integer,allocatable :: frag(:)
    character(len=:),allocatable :: wbofile
  end type lwoniom_input

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

    !> the actual "reading" part
    open (newunit=io,file=tomlfile)
    call toml_parse(table,io,error)
    close (unit=io)
    if (allocated(error)) then
      print*,"Error parsing table:"//error%message
    end if

    !> look for [lwoniom] block
    call get_value(table,"lwoniom",child, requested=.false.)
    if (.not.associated(child)) then
      write(stderr,'("**ERROR** ",a,a)') "No [lwoniom] section found in input file ",trim(tomlfile)
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
       write(stderr,'("**ERROR** ",a)') "Please provide the total number of atoms (natoms) in [lwoniom] block"
      return
    end if
    allocate(input%layer(input%nat), source = 0)
    allocate(input%frag(input%nat), source = 0)

    !> iterate over keys and check the max.layer and/or max.fragment number
    call table%get_keys(list)
    do ikey = 1, size(list)
      write(*,*) list(ikey)%key
    enddo

  end subroutine read_lwoniom_block
#endif

!========================================================================================!
!========================================================================================!
end module lwoniom_parse

