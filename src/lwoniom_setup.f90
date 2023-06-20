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
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use lwoniom_structures
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: lwoniom_data
  public :: lwoniom_initialize
  public :: print_lwoniom_results

!> this type bundles together most of the
!> data required for a lwONIOM calculation
  type :: lwoniom_data

    real(wp) :: lwoniom_energy
    
    !> number of layers
    integer :: nlayer = 0   
 
    !> number of fragments
    integer :: nfrag = 0
    type(structure_data),allocatable :: fragment(:)

  contains
    procedure :: deallocate => lwoniom_data_deallocate
  end type lwoniom_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!========================================================================================!

  subroutine print_lwoniom_results(iunit,dat)
    integer,intent(in) :: iunit ! file handle (usually output_unit=6)
    type(lwoniom_data),intent(in) :: dat
    character(len=*),parameter :: outfmt = &
                                  '(2x,a,f23.12,1x,a)'
    write (iunit,outfmt) "lwONIOM energy   ",dat%lwoniom_energy,"Eh   "
  end subroutine print_lwoniom_results
!========================================================================================!

  subroutine lwoniom_initialize(nat,at,xyz,dat,fname, &
  &                 print,verbose,iunit,iostat)
    character(len=*),parameter :: source = 'lwoniom_initialize'
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    character(len=*),intent(in) :: fname
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    integer,intent(out),optional :: iostat
    !> OUTPUT
    type(lwoniom_data),intent(inout) :: dat
    !> LOCAL
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun

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


!> lwONIOM calculator setup goes here
    if ((io /= 0).and.pr) then
      write (myunit,'("Could not create force field calculator ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if
  end subroutine lwoniom_initialize

!========================================================================================!
  subroutine lwoniom_data_deallocate(self)
    implicit none
    class(lwoniom_data) :: self
    self%lwoniom_energy = 0.0_Wp
    self%nlayer = 0
    self%nfrag = 0
    if (allocated(self%fragment)) deallocate (self%fragment)
  end subroutine lwoniom_data_deallocate

!========================================================================================!
end module lwoniom_setup

