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
  use lwoniom_interface
  use lwoniom_structures,only:lowercase
  use app_example
  implicit none

!> CMD parsing
  intrinsic :: iargc,getarg
  character(len=:),allocatable :: arg(:)
  integer :: nargs
  character(len=1024) :: cmd
  character(len=:),allocatable :: argument

!> The input data object
  type(lwoniom_input),allocatable :: inp

!> The xyz input file name
  character(len=:),allocatable :: fname

!> The lwONIOM data object
  type(lwoniom_data),allocatable :: dat

!> Other
  integer :: i,j,k,l,io,ich
  logical :: ex,dumpfrags

  character(len=*),parameter :: namespace = 'lwONIOM> '

!=======================================================================================!
!=======================================================================================!
!> get cmd and arguments
  call get_command(cmd)
  nargs = iargc()
  l = len_trim(cmd)
  allocate (arg(nargs),source=repeat(' ',l))
  do i = 1,nargs
    call getarg(i,arg(i))
  end do
!> minimal version printout
  do i = 1,nargs
    argument = lowercase(trim(arg(i)))
    select case (argument)
    case ('-v','--version')
!>--- version number printout
      write (stdout,'(a)') 'lwoniom-app v0.0.1'
      call exit(0)
    case ('-h','--help')
!>--- help printout
      call help()
      call exit(0)
    end select
  end do

!=======================================================================================!
!> printout
  call header()
  write (stdout,'(/,1x,a)') 'Command line input:'
  write (stdout,'(1x,a,a,/)') '$ ',trim(cmd)

!> defaults
  dumpfrags = .false.

!=======================================================================================!
!> Parse arguments
  do i = 1,nargs
    argument = lowercase(trim(arg(i)))
    l = len_trim(argument)
    if (i == 1.and.argument(l-4:) .eq. '.toml') then
      allocate (inp)
      call lwoniom_parse_inputfile(argument,inp)
      cycle
    end if
    select case (argument)
    case ('-i','--input')
!>--- XYZ input file name
      inquire (file=trim(arg(i+1)),exist=ex)
      if (ex) fname = trim(arg(i+1))

    case ('-d','--dump')
!>--- dump fragments logical
      dumpfrags = .true.

    case ('-e','--example')
!>--- create the prophyrine example and exit
      call lwoniom_write_example()

    case default
      if (argument(1:1) .eq. '-') then
        write (stdout,'(2a)') 'Unknown argument: ',trim(arg(i))
        write (stdout,'(a)') 'Please refer to the --help page'
        call exit(1)
      end if

    end select
  end do

!=======================================================================================!
!> Set up the QM/MM data from the parsed input
  if (allocated(inp)) then

!> last chance to read the XYZ
    if (allocated(fname)) then
      call inp%parse_xyz(fname)
    end if

!> create the calculator (dat)
    allocate (dat)
    call lwoniom_new_calculator(inp,dat)

  else
    write (stdout,'(a)') '**WARNING** No lwONIOM setup defined'
    write (stdout,'(a)') 'Please refer to the --help page'
    stop
  end if

!=======================================================================================!
!> Do what the program does ...
  if (allocated(dat)) then
!>--- print some info
    write(stdout,*)
    call dat%info()
    write(stdout,*)

!>--- dump xyz fragments
    if (dumpfrags .and. .not.inp%dump_frag) then
      write(stdout,'(a,a)') namespace,'dumping all fragments ...'
      do i=1,dat%nfrag
       write(stdout,'(1x,a,i0,a)',advance='no') 'fragment.',i,'.xyz'
      enddo
      write(stdout,*)
      call dat%dump_fragments()
    end if



  end if

!=======================================================================================!
  write(stdout,'(/,a,a)') namespace,'normal termination.'
!=======================================================================================!
!=======================================================================================!
end program lwoniom_app


!=======================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=======================================================================================!


subroutine header()
  use iso_fortran_env,only:stdout => output_unit
  implicit none
  character(8)  :: date   ! Date in YYYYMMDD format
  character(10) :: time   ! Time in HHMMSS.sss format
  character(5)  :: zone
  call DATE_AND_TIME(date,time,zone)
  write (stdout,'(3x,a)')
  write (stdout,'(3x,a)') '*********************************************'
  write (stdout,'(3x,a)') '*                                           *'
  write (stdout,'(3x,a)') '*            lwONIOM-app v0.0.1             *'
  write (stdout,'(3x,a)') '*                                           *'
  write (stdout,'(3x,a)') '*   (C) 2023 P. Wesołowski and P. Pracht    *'
  write (stdout,'(3x,a)') '*                                           *'
  write (stdout,'(3x,a)') '*********************************************'
  write (stdout,'(3x,a)') '  run started on '//date(:4)//'/'//date(5:6)//'/'//date(7:)// &
  & ' at '//time(:2)//':'//time(3:4)//':'//time(5:)
end subroutine header

!=======================================================================================!

subroutine help()
  use iso_fortran_env,only:stdout => output_unit
  implicit none
  call header()
  write (stdout,*)
  write (stdout,'(a)') 'Usage:'
  write (stdout,'(a)') 'lwoniom-app <INPUT.toml> [options]'
  write (stdout,'(a)')
  write (stdout,'(a)') 'Description:'
  write (stdout,'(a)') 'ONIOM is a scheme to calculate the energy and properties of a molecule'
  write (stdout,'(a)') 'using a multi-layered (and multi-centered) approach combining, e.g., quantum'
  write (stdout,'(a)') 'mechanical (QM) and molecular mechanics (MM) calculations.'
  write (stdout,'(a)') 'lwONIOM-app is a standalone application that performs the corresponding subsystem partition.'
  write (stdout,'(a)') 'The app requires an input file in the .toml format and coordinates in .xyz format'
  write (stdout,'(a)')
  write (stdout,'(a)') 'Options:'
  write (stdout,'(a)') '  -h, --help                  Display this help message and exit.'
  write (stdout,'(a)') '  -v, --version               Display the program name and version and exit.'
  write (stdout,'(a)') '  -i, --input input_file      Specify the input file containing molecular coordinates (.xyz format).'
  write (stdout,'(a)') '                              May also be specified in the <input>.toml file'
  write (stdout,'(a)') '  -d, --dump                  Dump all subsystems as .xyz file.'
  write (stdout,'(a)') '  -e, --example               Create an example input and exit.'
  write (stdout,'(a)')

  write(stdout,*)
  write(stdout,'(1x,a)') 'lwONIOM is distributed in the hope that it will be useful,    ' 
  write(stdout,'(1x,a)') 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
  write(stdout,'(1x,a)') 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.          '
  write(stdout,'(1x,a)') 'See the GNU Lesser General Public License for more details.   ' 

  write(stdout,*)
  write(stdout,'(a)') '--help exit.'  
end subroutine help
!=======================================================================================!
!=======================================================================================!

