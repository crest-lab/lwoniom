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
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use lwoniom_setup
  use lwoniom_engrad
  implicit none
  private

!> re-exports
  public :: lwoniom_data,lwoniom_initialize

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!


! TODO we could write some input reader functions and put them here,
!      but this is not too urgent.


!========================================================================================!
!========================================================================================!
end module lwoniom_interface

