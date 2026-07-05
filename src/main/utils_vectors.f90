!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module vectorutils
!
! This module contains utilities for manipulating vectors
!  and 3x3 matrices
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none
 public :: minmaxave,cross_product,cross_product3D,curl3D_epsijk,det
 public :: matrixinvert3D,rotatevec,unitvec,mag,make_perp_frame,jacobi_eigen_sym,rotation_to_align

 private

contains
!-------------------------------------------------------------------
!+
!  find min, max and average of an array
!+
!-------------------------------------------------------------------
subroutine minmaxave(x,xmin,xmax,xav,npts)
 integer, intent(in)  :: npts
 real,    intent(in)  :: x(npts)
 real,    intent(out) :: xmin,xmax,xav
 integer :: i

 xav = 0.
 xmin = huge(xmin)
 xmax = -xmin
 do i=1,npts
    xav = xav + x(i)
    xmin = min(xmin,x(i))
    xmax = max(xmax,x(i))
 enddo
 xav = xav/real(npts)

end subroutine minmaxave

!-------------------------------------------------------------------
!+
!  vector cross product
!+
!-------------------------------------------------------------------
pure subroutine cross_product3D(veca,vecb,vecc)
 real, intent(in)  :: veca(3),vecb(3)
 real, intent(out) :: vecc(3)

 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end subroutine cross_product3D

!-------------------------------------------------------------------
!+
!  vector cross product as a function
!+
!-------------------------------------------------------------------
pure function cross_product(a, b) result(c)
 real, intent(in) :: a(3), b(3)
 real :: c(3)

 c(1) = a(2)*b(3) - a(3)*b(2)
 c(2) = a(3)*b(1) - a(1)*b(3)
 c(3) = a(1)*b(2) - a(2)*b(1)

end function cross_product

!-------------------------------------------------------------------
!+
!  curl from the 3 x 3 gradient matrix
!+
!-------------------------------------------------------------------
pure subroutine curl3D_epsijk(gradAvec,curlA)
 real, intent(in)  :: gradAvec(3,3)
 real, intent(out) :: curlA(3)

 curlA(1) = gradAvec(2,3) - gradAvec(3,2)
 curlA(2) = gradAvec(3,1) - gradAvec(1,3)
 curlA(3) = gradAvec(1,2) - gradAvec(2,1)

end subroutine curl3D_epsijk

!----------------------------------------------------------------
!+
!  Inverts a 3x3 matrix
!+
!----------------------------------------------------------------
subroutine matrixinvert3D(A,Ainv,ierr)
 real,    intent(in)  :: A(3,3)
 real,    intent(out) :: Ainv(3,3)
 integer, intent(out) :: ierr
 real :: x0(3),x1(3),x2(3),result(3)
 real    :: det, ddet

 ierr = 0

 x0 = A(1,:)
 x1 = A(2,:)
 x2 = A(3,:)

 call cross_product3D(x1,x2,result)
 det = dot_product(x0,result)

 if (abs(det) > tiny(det)) then
    ddet = 1./det
 else
    ddet = 0.
    Ainv = 0.
    ierr = 1
    return
 endif

 Ainv(:,1) = result(:)*ddet
 call cross_product3D(x2,x0,result)
 Ainv(:,2) = result(:)*ddet
 call cross_product3D(x0,x1,result)
 Ainv(:,3) = result(:)*ddet

end subroutine matrixinvert3D

!----------------------------------------------------------------
!+
!  Determinant of a 3x3 matrix
!+
!----------------------------------------------------------------
real function det(A)
 real, intent(in) :: A(3,3)
 real :: x0(3),x1(3),x2(3),result(3)

 x0 = A(1,:)
 x1 = A(2,:)
 x2 = A(3,:)

 call cross_product3D(x1,x2,result)
 det = dot_product(x0,result)

end function det

!------------------------------------------------------------------------
!+
!  rotate a vector (u) around an axis defined by another vector (v)
!  by an angle (theta) using the Rodrigues rotation formula
!+
!------------------------------------------------------------------------
pure subroutine rotatevec(u,v,theta)
 real, intent(inout) :: u(3)
 real, intent(in)    :: v(3)
 real, intent(in)    :: theta
 real :: k(3),w(3)

 !--normalise v
 k = v/sqrt(dot_product(v,v))
 !--Rodrigues rotation formula
 call cross_product3D(k,u,w)
 u = u*cos(theta) + w*sin(theta) + k*dot_product(k,u)*(1-cos(theta))
end subroutine rotatevec

!------------------------------------------------------------------------
!+
!  return unit vector given a vector
!+
!------------------------------------------------------------------------
pure function unitvec(u) result(uhat)
 real, intent(in) :: u(3)
 real :: uhat(3),u2

 u2 = dot_product(u,u)
 if (u2 > tiny(0.)) then
    uhat = u/sqrt(u2)
 else
    uhat = (/0.,0.,1./)  ! arbitrary if vector is zero
 endif

end function unitvec

!------------------------------------------------------------------------
!+
!  magnitude of a vector
!+
!------------------------------------------------------------------------
pure function mag(u) result(umag)
 real, intent(in) :: u(3)
 real :: umag

 umag = sqrt(dot_product(u,u))

end function mag

!--------------------------------------------------------------------------------
! +
!   Build two perpendicular unit vectors {b, c} that complete a right–handed frame
! +
!--------------------------------------------------------------------------------
pure subroutine make_perp_frame(a, b, c)
 real, intent(in)  :: a(3)     ! arbitrary non-zero vector
 real, intent(out) :: b(3), c(3)

 real :: aa(3), inv_norm

 ! normalise a
 inv_norm = 1.0 / sqrt(sum(a*a))
 aa        = a * inv_norm          ! temporarily store â in c (a is intent in so we can't modify it)

 ! pick the largest component in magnitude
 select case (maxloc(abs(aa), dim=1))
 case (1)                          ! |a_x| is largest -> use y-axis
    b = (/ 0.0, 1.0, 0.0 /)
 case (2)                          ! |a_y| is largest -> use z-axis
    b = (/ 0.0, 0.0, 1.0 /)
 case default                      ! |a_z| is largest -> use x-axis
    b = (/ 1.0, 0.0, 0.0 /)
 end select

 ! make b perpendicular to a via Gram–Schmidt process (https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process)
 b = b - dot_product(b, aa) * aa
 inv_norm = 1.0 / sqrt(sum(b*b))
 b = b * inv_norm

 ! c = a x b
 call cross_product3D(aa, b, c)
end subroutine make_perp_frame

!----------------------------------------------------------------
!+
!  Eigenvalues and eigenvectors of a real symmetric matrix via the
!  classical Jacobi rotation method with threshold pivoting. Same
!  algorithm already used privately in src/utils/analysis_NSmerger.f90
!  (Numerical Recipes, Press et al.), exposed here as a shared public
!  utility. On output, elements of a above the diagonal are destroyed.
!  d returns the eigenvalues (unordered); v's columns are the matching
!  unit eigenvectors. nrot is the number of rotations used (informational).
!  Source: Numerical Recipes in Fortran 77, section 11.1.
!+
!----------------------------------------------------------------
subroutine jacobi_eigen_sym(a,n,np,d,v,nrot)
 integer, intent(in)    :: n,np
 integer, intent(out)   :: nrot
 real,    intent(inout) :: a(np,np)
 real,    intent(out)   :: d(np),v(np,np)
 integer, parameter :: nmax = 500
 integer :: i,ip,iq,j
 real :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

 do ip=1,n
    do iq=1,n
       v(ip,iq) = 0.
    enddo
    v(ip,ip) = 1.
 enddo
 do ip=1,n
    b(ip) = a(ip,ip)
    d(ip) = b(ip)
    z(ip) = 0.
 enddo

 nrot = 0
 do i=1,50
    sm = 0.
    do ip=1,n-1
       do iq=ip+1,n
          sm = sm + abs(a(ip,iq))
       enddo
    enddo
    if (sm == 0.) return

    if (i < 4) then
       tresh = 0.2*sm/n**2
    else
       tresh = 0.
    endif

    do ip=1,n-1
       do iq=ip+1,n
          g = 100.*abs(a(ip,iq))
          if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) .and. (abs(d(iq))+g == abs(d(iq)))) then
             a(ip,iq) = 0.
          elseif (abs(a(ip,iq)) > tresh) then
             h = d(iq)-d(ip)
             if (abs(h)+g == abs(h)) then
                t = a(ip,iq)/h
             else
                theta = 0.5*h/a(ip,iq)
                t = 1./(abs(theta)+sqrt(1.+theta**2))
                if (theta < 0.) t = -t
             endif
             c = 1./sqrt(1+t**2)
             s = t*c
             tau = s/(1.+c)
             h = t*a(ip,iq)
             z(ip) = z(ip)-h
             z(iq) = z(iq)+h
             d(ip) = d(ip)-h
             d(iq) = d(iq)+h
             a(ip,iq) = 0.
             do j=1,ip-1
                g = a(j,ip); h = a(j,iq)
                a(j,ip) = g-s*(h+g*tau)
                a(j,iq) = h+s*(g-h*tau)
             enddo
             do j=ip+1,iq-1
                g = a(ip,j); h = a(j,iq)
                a(ip,j) = g-s*(h+g*tau)
                a(j,iq) = h+s*(g-h*tau)
             enddo
             do j=iq+1,n
                g = a(ip,j); h = a(iq,j)
                a(ip,j) = g-s*(h+g*tau)
                a(iq,j) = h+s*(g-h*tau)
             enddo
             do j=1,n
                g = v(j,ip); h = v(j,iq)
                v(j,ip) = g-s*(h+g*tau)
                v(j,iq) = h+s*(g-h*tau)
             enddo
             nrot = nrot+1
          endif
       enddo
    enddo

    do ip=1,n
       b(ip) = b(ip)+z(ip)
       d(ip) = b(ip)
       z(ip) = 0.
    enddo
 enddo

end subroutine jacobi_eigen_sym

!------------------------------------------------------------------------
!+
!  Determine the rotation (axis + angle, Rodrigues convention) that maps
!  unit vector a onto unit vector b. Handles the degenerate cases where a
!  and b are already parallel (aligned=.true., no rotation needed; caller
!  should skip applying one) or exactly antiparallel (the cross product
!  is undefined, so the rotation is a 180 degree flip about an arbitrary
!  axis perpendicular to b, via make_perp_frame — correct for any choice
!  of perpendicular axis since a 180 degree rotation about any axis
!  perpendicular to a antiparallel pair maps one onto the other).
!+
!------------------------------------------------------------------------
subroutine rotation_to_align(a,b,rot_axis,rot_angle,aligned,tol)
 use physcon, only:pi
 real,    intent(in)  :: a(3),b(3)
 real,    intent(out) :: rot_axis(3),rot_angle
 logical, intent(out) :: aligned
 real,    intent(in), optional :: tol
 real :: cross(3),cosang,atol,perp(3),third(3)

 atol = 1.e-6
 if (present(tol)) atol = tol

 call cross_product3D(a,b,cross)
 cosang = max(-1.,min(1.,dot_product(a,b)))
 rot_angle = acos(cosang)
 aligned = .false.

 if (mag(cross) < atol) then
    if (cosang < 0.) then
       call make_perp_frame(b,perp,third)
       rot_axis  = perp
       rot_angle = pi
    else
       rot_axis  = (/0.,0.,1./)
       rot_angle = 0.
       aligned   = .true.
    endif
 else
    rot_axis = cross
 endif

end subroutine rotation_to_align

end module vectorutils
