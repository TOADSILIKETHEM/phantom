!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testlinalg
!
! Unit tests of linear algebra
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: io, linalg, testutils, vectorutils
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_linalg

 private
 real :: A(3,3)

contains
!--------------------------------------------
!+
!  Unit tests of the polynomial solvers
!+
!--------------------------------------------
subroutine test_linalg(ntests,npass)
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a)") '--> TESTING LINEAR ALGEBRA'

 call test_matrix_inversion(ntests,npass)
 call test_jacobi_eigen_sym(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- LINEAR ALGEBRA TESTS COMPLETE'

end subroutine test_linalg

!--------------------------------------------
!+
!  Unit tests of the quartic solver
!+
!--------------------------------------------
subroutine test_matrix_inversion(ntests,npass)
 use linalg, only:solve_bicgstab,get_Ax_interface
 use vectorutils, only:matrixinvert3D
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfail(2),n
 real, parameter :: tol = 1.e-6
 real :: Ainv(3,3),b(3),x(3),xexact(3)
 procedure(get_Ax_interface), pointer :: get_Ax_routine

 if (id==master) write(*,"(/,a)") '--> checking BiCGSTAB solver'

 A = reshape((/10.,2.,3.,&
               4.,-10.,2.,&
               2.,1.,-10./), (/3,3/))
 b = (/5.,0.,1./)
 x = (/0.,0.,0./) ! initial guess
 n = 3

 ! exact solution
 call matrixinvert3D(A,Ainv,ierr)
 xexact = matmul(Ainv,b)

 get_Ax_routine=>get_Ax_local
 call solve_bicgstab(n,b,get_Ax_routine,x,ierr)
 call checkval(ierr,0,0,nfail(1),'ierr=0')
 call checkval(3,x,xexact,tol,nfail(2),'xexact')

 call update_test_scores(ntests,nfail,npass)

contains

function get_Ax_local(n,x) result(Ax)
 integer, intent(in) :: n
 real,    intent(in) :: x(n)
 real :: Ax(n)

 Ax = matmul(A,x)

end function get_Ax_local

end subroutine test_matrix_inversion

!--------------------------------------------
!+
!  Unit test of the symmetric Jacobi eigenvalue solver, using a
!  matrix with an analytically known eigendecomposition: a decoupled
!  x-axis (eigenvalue 5) plus a 2x2 y-z block [[2,0.5],[0.5,3]] whose
!  eigenvalues/eigenvectors follow the standard closed form.
!+
!--------------------------------------------
subroutine test_jacobi_eigen_sym(ntests,npass)
 use vectorutils, only:jacobi_eigen_sym
 integer, intent(inout) :: ntests,npass
 integer :: nfail(6),nrot,i,j,imatch
 real :: Asym(3,3),d(3),v(3,3)
 real :: eval_expected(3),evec_expected(3,3)
 real :: dtol,best,diff

 if (id==master) write(*,"(/,a)") '--> checking Jacobi eigenvalue solver'

 Asym = reshape((/5.0,0.0,0.0, &
                  0.0,2.0,0.5, &
                  0.0,0.5,3.0/),(/3,3/))

 eval_expected = (/1.792893218813452,3.207106781186548,5.0/)
 evec_expected = reshape((/0.0,-0.923879532511287,0.382683432365090, &
                           0.0, 0.382683432365090,0.923879532511287, &
                           1.0, 0.0,               0.0            /),(/3,3/))

 call jacobi_eigen_sym(Asym,3,3,d,v,nrot)

 dtol = 1.e-5
 nfail = 0
 do i=1,3
    imatch = 1
    best = abs(d(i)-eval_expected(1))
    do j=2,3
       diff = abs(d(i)-eval_expected(j))
       if (diff < best) then
          best = diff
          imatch = j
       endif
    enddo
    call checkval(d(i),eval_expected(imatch),dtol,nfail(i),'eigenvalue matches analytic value')
    ! eigenvector direction is only defined up to an overall sign
    call checkval(abs(dot_product(v(:,i),evec_expected(:,imatch))),1.0,dtol,nfail(3+i),&
                  'eigenvector direction matches analytic axis')
 enddo

 call update_test_scores(ntests,nfail,npass)

end subroutine test_jacobi_eigen_sym

end module testlinalg
