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
  use lwoniom_hessian
  implicit none
  private

  real(wp),parameter :: bohr = 0.52917726_wp

!> RE-exports
  public :: lwoniom_input,lwoniom_parse_inputfile
  public :: lwoniom_data,lwoniom_initialize
  public :: lwoniom_singlepoint
  public :: lwoniom_gethess,lwonion_placehess 

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
    integer :: i,j,k,l,p,ii,jj,kk,ll
    integer,allocatable :: at(:),bond(:,:)
    logical :: ex,success

    if (inp%try_bin) then
      call dat%read_bin(success)
      write (stdout,'("lwONIOM> ",a,l4)') 'successfully read lwoniom.data?  ',success
      if (success) return
    end if

    if (allocated(inp%xyz).and. &
        allocated(inp%layer).and.allocated(inp%frag)) then

      allocate (at(inp%nat),source=0)
      if (allocated(inp%zsym)) then
        do i = 1,inp%nat
          at(i) = zsym_to_at(inp%zsym(i))
        end do
      else if (allocated(inp%at)) then
        at = inp%at
      else
        write (stderr,'("**ERROR** ",a)') "Required lwoniom_input atom type info missing!"
        error stop
      end if

      if (.not.allocated(inp%wbo)) then
        !> without bonding topology (will be set up from covalent radii)
        call lwoniom_initialize(inp%nat,at,inp%xyz/bohr,dat, &
        &                 inp%layer,inp%frag)

      else
        !> with user-defined bonding topology, needs to be passed as int
        allocate (bond(inp%nat,inp%nat),source=0)
        do i = 1,inp%nat
          do j = 1,inp%nat
            bond(i,j) = nint(inp%wbo(i,j))
          end do
        end do
        call lwoniom_initialize(inp%nat,at,inp%xyz/bohr,dat, &
        &                 inp%layer,inp%frag,bond)

      end if
      !deallocate (at)

      !> set up mapping of theory levels
      if (allocated(inp%layerlvl)) then
        dat%ncalcs = 2*dat%nfrag-1
        allocate (dat%calcids(2,dat%nfrag),source=0)
        do i = 1,dat%nfrag
          j = dat%fragment(i)%layer
          dat%calcids(1,i) = inp%layerlvl(j) !> high level
          p = dat%fragment(i)%parent
          if (p .ne. 0) then
            j = dat%fragment(p)%layer
            dat%calcids(2,i) = inp%layerlvl(j) !> low level, same layer as parent
          else
            dat%calcids(2,i) = dat%calcids(1,i) !> highest layer has only one level
          end if
        end do
      end if

      !> atom replacement in layers?
      if (allocated(inp%layerreplace)) then
        dat%replace_at = .true.
        do i = 1,dat%nlayer
          if (.not.any(inp%layerreplace(:,i) > 0)) cycle
          do j = 1,dat%nfrag
            k = dat%fragment(j)%layer
            if (i == k) then
              do ii = 1,size(inp%layerreplace,1)
                kk = inp%layerreplace(ii,i)
                if (kk > 0) then
                  do jj = 1,dat%fragment(j)%nat
                    if (dat%fragment(j)%at(jj) == ii) then
                      dat%fragment(j)%at(jj) = kk
                    end if
                  end do
                end if
              end do
            end if
          end do
        end do
      end if

      !> fragment specific charges
      if(allocated(inp%fragmentchrg))then
        do j = 1,dat%nfrag
          if(inp%fragmentchrg(j) /= -9990 )then
            if(.not.allocated(dat%fragment(j)%chrg))then
             allocate(dat%fragment(j)%chrg)
            endif
            dat%fragment(j)%chrg = inp%fragmentchrg(j)
          endif
        enddo
      endif

      inquire (file='lwoniom.data',exist=ex)
      if (inp%try_bin) then
        if (.not.ex) then
          write (stdout,'("lwONIOM> ",a)') 'writing lwONIOM model to "lwoniom.data" binary file ... '
          flush (stdout)
          call dat%dump_bin
        else
          write (stdout,'("lwONIOM> ",a)') 'lwONIOM model file "lwoniom.data" already exists and will NOT be overwritten.'
        end if
      end if

      write (stdout,'("lwONIOM> ",a)') 'setup done.'
  
      !> some additional stuff that can be done
      if(inp%dump_frag)then
       write (stdout,'("lwONIOM> ",a)') 'dumping fragment xyz files'
       do i=1,dat%nfrag
         write(stdout,'(1x,a,i0,a)',advance='no') 'fragment.',i,'.xyz'
       enddo
       write(stdout,*)
       call dat%dump_fragments()
      endif

      if(inp%dump_layer)then
       write (stdout,'("lwONIOM> ",a)') 'dumping layer xyz files'
       do i=1,dat%nlayer
         write(stdout,'(1x,a,i0,a)',advance='no') 'layer.',i,'.xyz'
       enddo
       write(stdout,*)
       call dat%dump_layers(at,inp%xyz)
      endif

  
    else
      write (stderr,'("**ERROR** ",a)') "Required lwoniom_input setup info missing!"
      return
    end if
    if(allocated(at)) deallocate(at)
    call inp%deallocate()
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

