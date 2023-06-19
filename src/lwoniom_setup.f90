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
  use iso_fortran_env,only:wp => real64,stdout => output_unit
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

!> variables go here

  contains
    procedure :: deallocate => lwoniom_data_deallocate
    procedure :: type_reset => lwoniom_data_reset_types
    procedure :: type_init => lwoniom_data_make_types
  end type lwoniom_data


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_singlepoint(nat,at,xyz,dat,energy,gradient,verbose,iostat)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: at(nat)    !> atom types
    real(wp),intent(in) :: xyz(3,nat) !> Cartesian coordinates in Bohr
    logical,intent(in),optional    :: verbose  !> printout activation 
    type(lwoniom_data),intent(inout) :: dat  !> collection of lwoniom datatypes and settings
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io
    logical :: pr

    !> printout activation via verbosity
    if(present(verbose))then
      pr = verbose
    else
      pr =.false. !> (there is close to no printout anyways)
    endif

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    io = 0


    if (present(iostat)) then
      iostat = io
    end if

  end subroutine lwoniom_singlepoint
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
    if(present(print))then
      pr = print
    else
      pr = .false.
    endif
    if(present(verbose))then
      pr2 = verbose
    else
      pr2 = .false.
    endif
    if(pr2) pr = pr2    
    if(present(iunit))then
      myunit = iunit
    else
      myunit = stdout
    endif

!> Reset datatypes
    call dat%type_init()


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

    !if (allocated(self%param)) deallocate (self%param)
  end subroutine lwoniom_data_deallocate
  subroutine lwoniom_data_reset_types(self)
    implicit none
    class(lwoniom_data) :: self

    !if (allocated(self%param)) deallocate (self%param)
  end subroutine lwoniom_data_reset_types
  subroutine lwoniom_data_make_types(self)
    implicit none
    class(lwoniom_data) :: self
    call self%type_reset()

    !allocate (self%param)
  end subroutine lwoniom_data_make_types

!========================================================================================!
end module lwoniom_setup

