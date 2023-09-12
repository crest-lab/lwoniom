module app_subprocess
  use iso_fortran_env,only:wp => real64,stdout=>output_unit

  real(wp),parameter,private :: autoaa = 0.52917721067d0
  real(wp),parameter,private :: aatoau = 1.0d0/autoaa

  public

!========================================================================================!
!========================================================================================!
contains !> module procedures start here
!========================================================================================!
!========================================================================================!

  subroutine app_subprocess_launch(fname,cmd,fragment,layer)
!***********************************************
!* Launch a subprocess. 
!* If the template cmd contains 'XXXX', this
!* substring will be replaced by fname.
!* Otherwise fname will just be appended to cmd
!***********************************************
     implicit none
     character(len=*),intent(in) :: fname
     character(len=*),intent(in) :: cmd
     integer,intent(in),optional :: fragment
     integer,intent(in),optional :: layer
     character(len=:),allocatable :: subprocess
     character(len=*),parameter :: dev0 = ' 2>/dev/null'
     integer :: i,j,k,l,io
     character(len=500) :: thisdir
     character(len=50) :: atmp
     logical :: subdir

     subdir = .false.
     if(present(fragment).and.present(layer))then 
       subdir = .true. 
       call getcwd(thisdir)
       !> create a fresh subdirectory
       write(atmp,'(a,i0,a,i0)') 'frag.',fragment,'.layer.',layer
       call system('rm -rf '//trim(atmp)//dev0)
       call system('mkdir '//trim(atmp)//' && '//'cp '//fname//' '//trim(atmp)//'/')
       call chdir(trim(atmp))
     endif

     k = index(cmd,'XXXX')
     if(k.ne.0)then
       subprocess = cmd(:k-1)//' '//trim(fname)//' '//cmd(k+3:)//dev0
     else
       subprocess = trim(cmd)//' '//trim(fname)//dev0
     endif
     write(stdout,'(a,a)') 'lwONIOM> ','executing subprocess:'
     write(stdout,'(1x,a)') subprocess
     !==============================!
     !> The actual subprocess call
     call system(subprocess,io)
     !==============================!
     if(io.ne.0)then
      write(stdout,'(a,a,i0)') 'lwONIOM> ','**ERROR** subprocess exit status ',io
     endif

     if(subdir)then
       call chdir(trim(thisdir))
     endif

  end subroutine app_subprocess_launch


  subroutine app_read_engrad(fname,natoms,energy,grad)
!***********************************************
!* Read the file fname containing the energy
!* and gradient elements for each of the natoms
!* One entry per line.
!* Comments (#) and empty lines are skipped
!***********************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: natoms
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grad(3*natoms)
    integer :: ich,j1,io,k
    character(len=200) :: line

    energy = 0.0_wp
    grad(:) = 0.0_wp

    ! process the output file
    open (newunit=ich,file=fname,status='old',action='read')
    io = 0
    k = 0
    do
      read (ich,'(a)',iostat=io) line
      if (io /= 0) exit !> EOF
      ! check for comments and empty lines
      line = adjustl(line)
      if (line(1:1) .eq. '#') cycle
      if (len_trim(line) .eq. 0) cycle
      ! check for the gradient section
      k = k+1
      if (k == 1) read (line,*) energy
      if (k > 1) then
        read (line,*) grad(k-1)
      end if
    end do
    close (ich)

  end subroutine app_read_engrad

!========================================================================================!
!========================================================================================!
end module app_subprocess
