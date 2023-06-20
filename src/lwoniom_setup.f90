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
    procedure :: add_fragment => lwoniom_add_fragment
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

  subroutine lwoniom_initialize(nat,at,xyz,dat, &
  &                 layer,center,               &
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
!*  B) If the "layer" and the "center" arrays are provided,
!*     one layer can have mutliple fragments (MC-ONIOM)
!*
!******************************************************************
    character(len=*),parameter :: source = 'lwoniom_initialize'
    !> INPUT
    integer,intent(in) :: nat          !> number of atoms    
    integer,intent(in) :: at(nat)      !> atom number for each atom
    real(wp),intent(in) :: xyz(3,nat)  !> Cartesian coordinates for each atom
    integer,intent(in) :: layer(nat)   !> layer instruction for each atom
    integer,intent(in),optional :: center(nat) !> center instruction (for MC-ONIOM)
    logical,intent(in),optional  :: print
    logical,intent(in),optional  :: verbose
    integer,intent(in),optional  :: iunit
    !> OUTPUT
    integer,intent(out),optional :: iostat
    type(lwoniom_data),intent(inout) :: dat
    !> LOCAL
    type(structure_data) :: tmp
    integer,allocatable :: center_tmp(:)
    integer :: ich,io,myunit
    logical :: ex,okbas,pr,pr2
    logical :: exitRun

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
    allocate(center_tmp(nat), source=0)
    if(present(center))then
        center_tmp(:) = center(:)
    endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  TO HERE                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> lwONIOM calculator setup goes here
    if ((io /= 0).and.pr) then
      write (myunit,'("Could not create lwONIOM object ",a)') source
    end if
    if (present(iostat)) then
      iostat = io
    end if
  end subroutine lwoniom_initialize

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
    if(frag%layer > self%nlayer)then
       self%nlayer = frag%layer
    endif
    if(self%nfrag == 0)then
      self%nfrag = 1
      allocate(self%fragment(1))
      self%fragment(1) = frag
    else
      n = self%nfrag
      n1 = n + 1
      allocate(tmp(n1))
      tmp(1:n) = self%fragment(1:n)
      tmp(n1) = frag
      call move_alloc(tmp,self%fragment) 
      self%nfrag = n1
    endif
  end subroutine lwoniom_add_fragment



!========================================================================================!
end module lwoniom_setup

