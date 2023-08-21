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
program lwoniom_main_tester
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use testmol
  use lwoniom_interface
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
  real(wp),allocatable :: refgrad(:,:)
  real(wp) :: gnorm

  logical :: fail,pr
  integer :: io
  type(lwoniom_data) :: dat
  integer,allocatable :: layer(:)

!========================================================================================!
  fail = .false.
  pr = .true.

  nat = testnat
  allocate (at(nat),xyz(3,nat))
  at = testat
  xyz = testxyz
  chrg = 0

  energy = 0.0_wp
  gnorm = 0.0_wp
  allocate (gradient(3,nat),source=0.0_wp)

  write (*,*) nat
  write (*,*)
  do i = 1,nat
    write (*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
  end do
  call writetestcoord()
  call writetestxyz()

!=======================================================================================!
!=======================================================================================!
!> STANDARD USAGE
!=======================================================================================!
!=======================================================================================!

  write (*,*)
  write (*,*) '============================================================'
  write (*,*) '==================== lwONIOM SETUP ========================='
  write (*,*) '============================================================'
  write (*,*)
    
    allocate(layer(nat))
    layer = at
    call lwoniom_initialize(nat,at,xyz,dat,layer)
    call dat%info()

  write (*,*)
  write (*,*) '========================== END ============================='
  write (*,*) '==================== lwONIOM SETUP ========================='
  write (*,*) '========================== END ============================='

!=======================================================================================!

  write(*,*) 'dumping fragments ... '
  call dat%dump_fragments()

!=======================================================================================!
  allocate(refgrad(3,nat),source=0.0_wp)
  do i=1,dat%nfrag
  dat%fragment(i)%grd = 1.0d0
  dat%fragment(i)%linkgrd = 1.0d0
   call dat%fragment(i)%jacobian(nat,gradient)
  write (*,*) nat,dat%fragment(i)%nat,size(gradient,2),sum(gradient)/3.0d0
    refgrad = 0.0d0
    do j=1,dat%fragment(i)%nat
      k = dat%fragment(i)%opos(j)
      refgrad(:,k) = 1.0d0
    enddo
    write(*,*) 'gradient difference',sum(gradient(:,:)-refgrad(:,:))
  enddo


!=======================================================================================!
  deallocate (gradient)
  deallocate (xyz,at)
!=======================================================================================!
end program lwoniom_main_tester
