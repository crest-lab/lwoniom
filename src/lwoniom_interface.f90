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
module lwoniom_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  use lwoniom_parse
  use lwoniom_setup
  use lwoniom_structures,only:zsym_to_at
  use lwoniom_engrad
  implicit none
  private

  real(wp),parameter :: bohr = 0.52917726_wp

!> RE-exports
  public :: lwoniom_input,lwoniom_parse_inputfile
  public :: lwoniom_data,lwoniom_initialize
  public :: lwoniom_singlepoint

!> interface routines
  public :: lwoniom_new_calculator
  interface lwoniom_new_calculator
    module procedure :: lwoniom_new_calculator_file
    module procedure :: lwoniom_new_calculator_inp
  end interface lwoniom_new_calculator

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_new_calculator_inp(inp,dat)
!*****************************************************************
!* This routine will attempt create the lwoniom_data object 'dat'
!* from a previously parsed lwoniom_input object containing the
!* setup information
!****************************************************************
    implicit none
    !> input file name, TOML format
    type(lwoniom_input),intent(inout) :: inp
    !> the calculator/data type to be created
    type(lwoniom_data),intent(out) :: dat
    !> LOCAL
    integer :: i,j,k,l
    integer,allocatable :: at(:),bond(:,:)
    logical :: ex

    if (allocated(inp%xyz).and.allocated(inp%zsym).and. &
        allocated(inp%layer).and.allocated(inp%frag)) then

      allocate (at(inp%nat),source=0)
      do i = 1,inp%nat
        at(i) = zsym_to_at(inp%zsym(i))
      end do

      if (.not.allocated(inp%wbo)) then
        !> without bonding topology (will be set up from covalent radii)
        call lwoniom_initialize(inp%nat,at,inp%xyz/bohr,dat, &
        &                 inp%layer)

      else
        !> with user-defined bonding topology
        allocate (bond(inp%nat,inp%nat),source=0)
        do i = 1,inp%nat
          do j = 1,inp%nat
            bond(i,j) = nint(inp%wbo(i,j))
          end do
        end do
        call lwoniom_initialize(inp%nat,at,inp%xyz/bohr,dat, &
        &                 inp%layer,inp%frag,bond)

      end if
    else
      write (stderr,'("**ERROR** ",a)') "Required lwoniom_input setup info missing!"
      return
    end if

  end subroutine lwoniom_new_calculator_inp

!========================================================================================!

  subroutine lwoniom_new_calculator_file(inputfile,dat)
!**************************************************************
!* This routine will attempt to read a TMOL input file
!* specified as 'inputfile' and create the
!* lwoniom_data object 'dat'.
!* This is a wrapper for lwoniom_new_calculator_inp.
!* Obviously, ALL information must be present in the inputfile
!**************************************************************
    implicit none
    !> input file name, TOML format
    character(len=*),intent(in) :: inputfile
    !> the calculator/data type to be created
    type(lwoniom_data),intent(out) :: dat
    !> LOCAL
    type(lwoniom_input) :: inp  !> temporary storage
    integer :: i,j,k,l
    integer,allocatable :: at(:),bond(:,:)
    logical :: ex

    inquire (file=inputfile,exist=ex)
    if (.not.ex) then
      write (stderr,'("**ERROR** ",3a)') "Input file ",trim(inputfile)," was not found!"
      return
    end if
    call lwoniom_parse_inputfile(inputfile,inp)
    call lwoniom_new_calculator_inp(inp,dat)

    call inp%deallocate()
  end subroutine lwoniom_new_calculator_file

!========================================================================================!
!========================================================================================!
end module lwoniom_interface

