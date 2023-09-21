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
  use lwoniom_parse
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
  logical,allocatable :: layer(:,:)
  character(len=*),parameter :: green = char(27)//'[92m'
  character(len=*),parameter :: red = char(27)//'[91m'
  character(len=*),parameter :: escape = char(27)//'[0m'
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
  write (*,*) '==========================================================='
  write (*,*) '==================== lwONIOM TEST ========================='
  write (*,*) '==========================================================='
  write (*,*)

  write (stdout,*) 'Setting up 5-layer ONIOM model system ... '
  write (stdout,*) 'Expecting a bunch of "Bad setup" warnings ...'
  write (stdout,*)

  !> some random ONIOM-5 setup
  allocate (layer(nat,5))
  do i = 1,4
    j = (i-1)*2+1
    k = 6-i
    layer(1:j,k) = .true.
  end do
  layer(:,1) = .true.
  call lwoniom_initialize(nat,at,xyz,dat,layer)
  write (stdout,*)
  call dat%info()

  write (stdout,*)
  write (stdout,*) repeat('-',50)
  write (stdout,*) 'Setting fragment energies gradients to test values ... '
  do i = 1,dat%nfrag
    dat%fragment(i)%energy_high = -10.0d0
    dat%fragment(i)%energy_low = -10.0d0
    dat%fragment(i)%grd_high = 1.0d0
    dat%fragment(i)%grd_low = 1.0d0
    if (allocated(dat%fragment(i)%linkgrd_high)) then
      dat%fragment(i)%linkgrd_high = 0.5d0
      dat%fragment(i)%linkgrd_low = 0.5d0
    end if
  end do
  write (stdout,*) 'Gradient projections via Jacobian ... '
  do i = 1,dat%nfrag
    write (stdout,'(2x,a,i0,a)',advance='no') 'fragment ',i,' ... '
    flush (stdout)
    call dat%fragment(i)%jacobian(nat)
    call dumpgrad(i,nat,dat%fragment(i)%gradient_low)
    write (stdout,'(a)') 'done.'
  end do

  write (stdout,'(1x,a)',advance='no') 'Reconstructing the full "ONIOM5" gradient ...'
  flush (stdout)
  call lwoniom_singlepoint(nat,dat,energy,gradient)
  !write(stdout,'(1x,a,f25.15)') 'energy: ',energy
  !write(stdout,'(1x,a)') 'gradient:'
  !do i=1,nat
  !  write(stdout,'(3f25.15)') gradient(1:3,i)
  !enddo
  write (stdout,*) 'done.'
  allocate (refgrad(3,nat),source=1.0d0)
  fail = any(abs(refgrad(:,:)-gradient(:,:)) .gt. 1.0d-10)
  call dat%deallocate()
  write (*,*)
  write (*,*) '========================= END ============================='
  write (*,*) '==================== lwONIOM TEST ========================='
  write (*,*) '========================= END ============================='

  write (stdout,*)
  write (stdout,'(a)',advance='no') '=> All tests passed:  '
  flush (stdout)
  if (fail) then
    write (stdout,*) red,'no!',escape
  else
    write (stdout,*) green,'yes!',escape
  end if
  write (stdout,*)
!=======================================================================================!
  deallocate (gradient)
  deallocate (xyz,at)
!=======================================================================================!
end program lwoniom_main_tester

subroutine dumpgrad(i,nat,grad)
  use iso_fortran_env,only:wp => real64
  implicit none
  integer :: i,nat
  real(wp) :: grad(3,nat)
  integer :: ich,j
  character(len=40) :: atmp
  write (atmp,'(a,i0,a)') 'fragment.',i,'.engrad'
  open (newunit=ich,file=trim(atmp))
  do j = 1,nat
    write (ich,'(3f25.15)') grad(1:3,j)
  end do
  close (ich)
end subroutine dumpgrad
