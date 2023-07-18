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
program lwoniom_app
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  implicit none

  integer :: nat
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  integer :: chrg
  integer :: uhf
  integer :: i,j,k,l

!========================================================================================!
  real(wp) :: energy
  real(wp),allocatable :: gradient(:,:)
  real(wp) :: gnorm
  integer,allocatable :: bo(:,:)

  logical :: fail,pr
  integer :: io
  type(lwoniom_data) :: dat
  integer,allocatable :: layer(:)
  integer,allocatable :: subsystem(:)

!========================================================================================!
!> Initialize
  fail = .false.
  pr = .true.
  energy = 0.0_wp
  gnorm = 0.0_wp

!=========================================================================================!
!> Parse arguments
!TODO

!=========================================================================================!
!> Read structure information

!TODO
!  nat = 
!  allocate (at(nat), source=0)
!  allocate (xyz(3,nat), source=0.0_wp)
!  allocate (gradient(3,nat),source=0.0_wp)
!  at = 
!  xyz = 
!  chrg = 
!  uhf = 

!> optionally read bond info
!TODO
! if(
!   allocate(bo(nat,nat))
!
! endif

!> Read layerand/or subsystem  information
  allocate(layer(nat) , source = 0) 
!TODO
! if(
!   allocate(subsystem(nat), source = 0)
!
! endif

!> Print
  write (*,*) nat
  write (*,*)
  do i = 1,nat
    write (*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
  end do

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

  write (*,*)
  write (*,*) '============================================================'
  write (*,*) '==================== lwONIOM SETUP ========================='
  write (*,*) '============================================================'
  write (*,*) ' Layer information'

  if(allocated(subsystem) .and. allocated(bo))then
    call lwoniom_initialize(nat,at,xyz,dat,layer,subsystem,bo,.true.,iostat=io)
  elseif(allocated(bo))then
    call lwoniom_initialize(nat,at,xyz,dat,layer,bond=bo,print=.true.,iostat=io)
  else
    call lwoniom_initialize(nat,at,xyz,dat,layer,print=.true.,iostat=io)
  endif
  call dat%info()

!=======================================================================================!

  write (*,*)
  write (*,*) '============================================================'
  write (*,*) '================== lwONIOM CALCULATIONS ===================='
  write (*,*) '============================================================'

!TODO


!=======================================================================================!

  write (*,*)
  write (*,*) '============================================================'
  write (*,*) '=============== lwONIOM GRADIENT CONSTRUCTION =============='
  write (*,*) '============================================================'

!TODO

!=======================================================================================!
  if(allocated(subsystem)) deallocate(subsystem)
  if(allocated(bo)) deallocate(bo)
  deallocate (gradient)
  deallocate (xyz,at)
!=======================================================================================!
end program lwoniom_app
