!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module read_obj
!
! Reads a Wavefront OBJ mesh and tests whether a point lies inside
! a closed triangulated surface using ray casting (Moller-Trumbore)
!
! :References: Moller & Trumbore (1997) J. Graphics Tools 2(1):21-28
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: read_obj_file, scale_and_centre_obj, point_in_polyhedron

 private

contains

!----------------------------------------------------------------
!+
!  Parse a Wavefront OBJ file, extracting vertices and triangular
!  faces. Polygons with more than 3 vertices are fan-triangulated.
!  Handles simple (f v1 v2 v3) and slash-delimited
!  (f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3) face formats.
!+
!----------------------------------------------------------------
subroutine read_obj_file(filename, verts, faces, nverts, nfaces)
 character(len=*),     intent(in)  :: filename
 real,    allocatable, intent(out) :: verts(:,:)   ! (3, nverts)
 integer, allocatable, intent(out) :: faces(:,:)   ! (3, nfaces), 1-based vertex indices
 integer,              intent(out) :: nverts, nfaces
 integer, parameter :: iunit = 42
 integer, parameter :: max_poly = 64
 character(len=512) :: line
 character(len=1)   :: c1, c2
 integer :: ios, nv, nf, nvids, vids(max_poly), i
 real    :: x, y, z

 nverts = 0
 nfaces = 0

 ! --- Pass 1: count vertices and triangles ---
 open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
 if (ios /= 0) then
    print "(a)", ' ERROR read_obj_file: cannot open '//trim(filename)
    allocate(verts(3,1), faces(3,1))
    return
 endif
 nv = 0
 nf = 0
 do
    read(iunit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    line = adjustl(line)
    if (len_trim(line) == 0) cycle
    c1 = line(1:1)
    if (len_trim(line) >= 2) then
       c2 = line(2:2)
    else
       c2 = ' '
    endif
    if (c1 == 'v' .and. c2 == ' ') then
       nv = nv + 1
    elseif (c1 == 'f') then
       call count_face_verts(line, nvids)
       if (nvids >= 3) nf = nf + (nvids - 2)
    endif
 enddo
 close(iunit)

 nverts = nv
 nfaces = nf
 allocate(verts(3, max(1, nverts)))
 allocate(faces(3, max(1, nfaces)))

 ! --- Pass 2: read data ---
 nv = 0
 nf = 0
 open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
 do
    read(iunit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    line = adjustl(line)
    if (len_trim(line) == 0) cycle
    c1 = line(1:1)
    if (len_trim(line) >= 2) then
       c2 = line(2:2)
    else
       c2 = ' '
    endif
    if (c1 == 'v' .and. c2 == ' ') then
       read(line(3:), *, iostat=ios) x, y, z
       if (ios == 0 .and. nv < nverts) then
          nv = nv + 1
          verts(1, nv) = x
          verts(2, nv) = y
          verts(3, nv) = z
       endif
    elseif (c1 == 'f') then
       call parse_face_line(line, vids, nvids, max_poly)
       do i = 2, nvids - 1
          if (nf < nfaces) then
             nf = nf + 1
             faces(1, nf) = vids(1)
             faces(2, nf) = vids(i)
             faces(3, nf) = vids(i+1)
          endif
       enddo
    endif
 enddo
 close(iunit)

 print "(a,i0,a,i0,a)", ' read_obj_file: ', nverts, ' vertices, ', nfaces, ' triangles'

end subroutine read_obj_file

!----------------------------------------------------------------
!+
!  Centre the mesh at the origin and scale so that the maximum
!  vertex radius equals r_target
!+
!----------------------------------------------------------------
subroutine scale_and_centre_obj(verts, nverts, r_target)
 real,    intent(inout) :: verts(:,:)
 integer, intent(in)    :: nverts
 real,    intent(in)    :: r_target
 integer :: i
 real    :: cx, cy, cz, rmax, r

 if (nverts <= 0) return

 cx = sum(verts(1, 1:nverts)) / real(nverts)
 cy = sum(verts(2, 1:nverts)) / real(nverts)
 cz = sum(verts(3, 1:nverts)) / real(nverts)

 do i = 1, nverts
    verts(1,i) = verts(1,i) - cx
    verts(2,i) = verts(2,i) - cy
    verts(3,i) = verts(3,i) - cz
 enddo

 rmax = 0.
 do i = 1, nverts
    r = sqrt(verts(1,i)**2 + verts(2,i)**2 + verts(3,i)**2)
    if (r > rmax) rmax = r
 enddo

 if (rmax > 0.) then
    verts(:, 1:nverts) = verts(:, 1:nverts) * (r_target / rmax)
 endif

end subroutine scale_and_centre_obj

!----------------------------------------------------------------
!+
!  Test whether a point xyz_local (expressed relative to the mesh
!  centroid) lies inside a closed triangulated surface, using the
!  Moller-Trumbore ray-casting algorithm along a perturbed +X ray.
!  Returns .true. if the point is inside.
!+
!----------------------------------------------------------------
logical function point_in_polyhedron(xyz_local, verts, faces, nverts, nfaces)
 real,    intent(in) :: xyz_local(3)
 real,    intent(in) :: verts(:,:)
 integer, intent(in) :: faces(:,:)
 integer, intent(in) :: nverts, nfaces
 ! Slightly off-axis ray avoids degenerate edge/vertex hits
 real, parameter :: ray_eps = 1.e-4
 real, parameter :: tol = 1.e-10
 real :: d(3), e1(3), e2(3), h(3), s(3), q(3)
 real :: a, f, u, v, t
 integer :: k, v0, v1, v2, nhits

 ! ray direction: nearly +X, avoids axis-aligned degeneracies
 d(1) = 1. + ray_eps
 d(2) = ray_eps
 d(3) = ray_eps

 nhits = 0
 do k = 1, nfaces
    v0 = faces(1,k)
    v1 = faces(2,k)
    v2 = faces(3,k)
    if (v0 < 1 .or. v0 > nverts) cycle
    if (v1 < 1 .or. v1 > nverts) cycle
    if (v2 < 1 .or. v2 > nverts) cycle

    e1(1) = verts(1,v1) - verts(1,v0)
    e1(2) = verts(2,v1) - verts(2,v0)
    e1(3) = verts(3,v1) - verts(3,v0)

    e2(1) = verts(1,v2) - verts(1,v0)
    e2(2) = verts(2,v2) - verts(2,v0)
    e2(3) = verts(3,v2) - verts(3,v0)

    ! h = d x e2
    h(1) = d(2)*e2(3) - d(3)*e2(2)
    h(2) = d(3)*e2(1) - d(1)*e2(3)
    h(3) = d(1)*e2(2) - d(2)*e2(1)

    a = e1(1)*h(1) + e1(2)*h(2) + e1(3)*h(3)
    if (abs(a) < tol) cycle   ! ray parallel to triangle

    f = 1.0 / a
    s(1) = xyz_local(1) - verts(1,v0)
    s(2) = xyz_local(2) - verts(2,v0)
    s(3) = xyz_local(3) - verts(3,v0)

    u = f * (s(1)*h(1) + s(2)*h(2) + s(3)*h(3))
    if (u < 0. .or. u > 1.) cycle

    ! q = s x e1
    q(1) = s(2)*e1(3) - s(3)*e1(2)
    q(2) = s(3)*e1(1) - s(1)*e1(3)
    q(3) = s(1)*e1(2) - s(2)*e1(1)

    v = f * (d(1)*q(1) + d(2)*q(2) + d(3)*q(3))
    if (v < 0. .or. u + v > 1.) cycle

    t = f * (e2(1)*q(1) + e2(2)*q(2) + e2(3)*q(3))
    if (t > tol) nhits = nhits + 1
 enddo

 point_in_polyhedron = (mod(nhits, 2) == 1)

end function point_in_polyhedron

!----------------------------------------------------------------
! internal helpers
!----------------------------------------------------------------

subroutine count_face_verts(line, nvids)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: nvids
 integer :: i, n
 logical :: in_token

 n = len_trim(line)
 nvids = 0
 in_token = .false.
 i = 2   ! skip 'f'
 do while (i <= n)
    if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
       in_token = .false.
    else
       if (.not.in_token) then
          nvids = nvids + 1
          in_token = .true.
       endif
    endif
    i = i + 1
 enddo
end subroutine count_face_verts

subroutine parse_face_line(line, vids, nvids, maxv)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: vids(:), nvids
 integer,          intent(in)  :: maxv
 integer :: i, n, istart, iend, ios, vid, islash
 character(len=32) :: token

 n = len_trim(line)
 nvids = 0
 i = 2   ! skip 'f'

 do while (i <= n .and. nvids < maxv)
    ! skip whitespace
    do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == char(9)))
       i = i + 1
    enddo
    if (i > n) exit

    ! gather non-whitespace token
    istart = i
    do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= char(9))
       i = i + 1
    enddo
    iend = i - 1

    token = line(istart:iend)
    islash = index(token, '/')
    if (islash > 1) then
       read(token(1:islash-1), *, iostat=ios) vid
    elseif (islash == 0) then
       read(token, *, iostat=ios) vid
    else
       ios = 1   ! malformed token (slash at position 1)
    endif

    if (ios == 0) then
       nvids = nvids + 1
       vids(nvids) = vid
    endif
 enddo

end subroutine parse_face_line

end module read_obj
