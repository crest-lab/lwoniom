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
module lwoniom_hessian
!************************************************************
!* This module implements the ONIOM Hessian reconstruction.
!************************************************************
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use lwoniom_setup
  use lwoniom_structures
  implicit none
  private

!> routines/datatypes that can be seen outside the module
  public :: lwoniom_gethess,lwonion_placehess

  logical,parameter :: debug = .false. 

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lwoniom_gethess(nat,dat,Hess,verbose,iostat)
!*************************************************
!* Obtain the MC-ONIOM gradient from dat.
!* The high- and low-level energies and gradients
!* for each fragment must be known at this point
!*************************************************
    implicit none
    !> INPUT
    integer,intent(in)  :: nat                 !> number of atoms in the original system
    logical,intent(in),optional    :: verbose  !> printout activation
    type(lwoniom_data),intent(inout) :: dat    !> collection of lwoniom datatypes and settings
    !> OUTPUT
    real(wp),intent(inout) :: Hess(3*nat,3*nat)
    integer,intent(out),optional  :: iostat
    !> LOCAL
    integer :: io,i,j,k,kk,root_id,n3
    integer :: nchilds
    logical :: pr

    !> printout activation via verbosity
    if (present(verbose)) then
      pr = verbose
    else
      pr = .false. !> (there is close to no printout anyways)
    end if

    io = 0
    n3 = 3*nat

    !> The Hessian reconstruction is done by calling the
    !> recursive hessian_recursion subroutine.
    !> We start at the root node (layer 1), and  the recursion terminates
    !> if a subsystem hasn't any child nodes
    root_id = dat%root_id
    call hessian_recursion(dat,root_id,nat)

    !> additive Hessian
    call unpack_symmetric_matrix(dat%fragment(root_id)%Hss_qq,Hess,n3)

    if (present(iostat)) then
      iostat = io
    end if

  end subroutine lwoniom_gethess

!========================================================================================!

  recursive subroutine hessian_recursion(dat,F,nat)
!*********************************************************************
!* Recursion construction of the MC-ONIOM Hessian
!* The Hessian of the current node H_i IN THE BASIS OF THE
!* ORIGINAL SYSTEM is calculated as:
!*
!*   H_i = H_i^h + ∑_j(H_j - H_j^l)
!*
!* where H_i^h is the high-level Hessian.
!* To this, (corrected) Hessian contibutions of all
!* children nodes H_j, minus their low-level contibution H_j^l.
!* Keep in mind that the high-level of the parent node (H_i^h) refers
!* to the same as the low-level of the child node (H_j^l).
!* H_j is itself calculated by the same equation as H_i, which is where
!* the recursion comes into play. The recusion terminates when
!* there are no child nodes. In this case H_i = H_i^h.
!*********************************************************************
    implicit none
    type(lwoniom_data),intent(inout) :: dat
    integer,intent(in) :: F
    integer,intent(in) :: nat
    integer :: nchilds,k,KK,n3,hdim

    n3 = 3*nat
    hdim = n3*(n3+1)/2

    if(debug) write(*,*) 'Constructing ONIOM Hessian for fragment',F
    if (.not.allocated(dat%fragment(F)%Hss_qq)) then
      allocate (dat%fragment(F)%Hss_qq(hdim))
    end if
    if (allocated(dat%fragment(F)%Hss_high)) then
      dat%fragment(F)%Hss_qq(:) = dat%fragment(F)%Hss_high(:)
    else
!>--- This fallback condition is implemented for the parent fragment/highest layer
!>--- which only has a low-level energy/gradient
      dat%fragment(F)%Hss_qq(:) = dat%fragment(F)%Hss_low(:)
    end if

    !> loop over all children nodes
    if (allocated(dat%fragment(F)%child)) then
      nchilds = size(dat%fragment(F)%child(:),1)
      do k = 1,nchilds
        KK = dat%fragment(F)%child(k) !> the child-node ID
        !==============================================!
        call hessian_recursion(dat,KK,nat) !> RECURSION
        !==============================================!

        dat%fragment(F)%Hss_qq(:) = dat%fragment(F)%Hss_qq(:) + &
        & dat%fragment(KK)%Hss_qq(:)-dat%fragment(kk)%Hss_low(:)

      end do
    else
      !> RECURSION TERMINATION: NO FURTHER CHILDREN NODES
      !> i.e., we take Hss_qq for this node 'as is'
      return
    end if

  end subroutine hessian_recursion

!========================================================================================!

  subroutine unpack_symmetric_matrix(packed_matrix,matrix,n)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: packed_matrix((n*(n+1))/2)
    real(wp),intent(inout) :: matrix(n,n)
    integer :: i,j,k
    k = 1
    do j = 1,n
      do i = 1,j
        matrix(i,j) = matrix(i,j)+packed_matrix(k)
        matrix(j,i) = matrix(i,j)!+packed_matrix(k)  !> lower triangle
        k = k+1
      end do
    end do
  end subroutine unpack_symmetric_matrix

!========================================================================================!

  subroutine pack_symmetric_matrix(packed_matrix,matrix,n,ulo)
!**************************************
!* Pack a symmetric matrix and only
!* store the upper or lower triangular 
!* form as a 1-d array
!**************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(out) :: packed_matrix((n*(n+1))/2)
    real(wp),intent(in) :: matrix(n,n)
    character(len=1),intent(in),optional :: ulo
    integer :: i,j,k
    logical :: uplo
    k = 1
    if (present(ulo)) then
      select case (ulo)
      case ('u','U'); uplo = .true.
      case ('l','L'); uplo = .false.
      end select
    else
      uplo = .true.
    end if
    do j = 1,n
      do i = 1,j
        if (uplo) then
          packed_matrix(k) = matrix(i,j) !> upper triangle
        else
          packed_matrix(k) = matrix(j,i) !> lower triangle
        end if
        k = k+1
      end do
    end do
  end subroutine pack_symmetric_matrix

!========================================================================================!

  subroutine lwonion_placehess(dat,truenat,nat,hess,F,highlow)
!****************************************************
!* Place a Hessian for a fragment with nat atoms
!* into fragment F of dat. highlow specfies if the
!* Hessian corresponds to the high or low level of F
!****************************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: dat
    integer,intent(in) :: truenat,nat
    real(wp),intent(in) :: hess(nat*3,nat*3)
    integer,intent(in) :: F,highlow
    !> LOCAL
    real(wp),allocatable :: JacoT(:,:)
    real(wp),allocatable :: hesstmp(:,:)
    integer :: i,j,k,l
    integer :: fragnat,m,n,hdim,n3

    m = nat
    n = truenat
    n3 = n*3
    hdim = n3*(n3+1)/2
    if (.not.allocated(dat%fragment(F)%Jaco)) then
      call dat%fragment(F)%construct_jacobian(m,n)
    end if
    allocate (JacoT(n3,m*3),source=0.0_wp)
    JacoT = transpose(dat%fragment(F)%Jaco)

    allocate (hesstmp(n3,n3),source=0.0_wp)
    hesstmp = matmul(matmul(JacoT,hess),dat%fragment(F)%Jaco)

    select case (highlow)
    case (1,3) !> high or root
      if (.not.allocated(dat%fragment(F)%Hss_high)) then
        allocate (dat%fragment(F)%Hss_high(hdim))
      end if
      call pack_symmetric_matrix( &
      & dat%fragment(F)%Hss_high, &
      & hesstmp,n3,'u')
      !if(highlow == 1) dat%fragment(F)%Hss_high = 0.0_wp
    case (2) !> low
      if (.not.allocated(dat%fragment(F)%Hss_low)) then
        allocate (dat%fragment(F)%Hss_low(hdim))
      end if
      call pack_symmetric_matrix( &
      & dat%fragment(F)%Hss_low,  &
      & hesstmp,n3,'u')
      !if(highlow == 2) dat%fragment(F)%Hss_high = 0.0_wp
    end select

    deallocate (hesstmp)
    deallocate (JacoT)
  end subroutine lwonion_placehess

!========================================================================================!
!========================================================================================!
end module lwoniom_hessian

