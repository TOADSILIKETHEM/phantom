!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup asteroid orbits using data from the IAU Minor Planet Center
!
! :References: https://minorplanetcenter.net/data
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - asteroids  : *add distant minor bodies as km-sized dust particles*
!   - dtmax_in   : *time between dumps (e.g. 1 hr)*
!   - epoch      : *epoch to query ephemeris, YYYY-MMM-DD HH:MM:SS.fff, blank = today*
!   - np_apophis : *number of particles used to represent apophis (0=none; 1=sink; n=gas)*
!   - tmax_in    : *end time of simulation (e.g. 3 days)*
!
! :Dependencies: centreofmass, eos_tillotson, infile_utils, io, kernel,
!   options, part, physcon, setbinary, setsolarsystem, setup_params,
!   spherical, timestep, units
!
 implicit none
 public :: setpart

 integer :: np_apophis
 logical :: asteroids
 character(len=20) :: epoch,tmax_in,dtmax_in
 logical :: use_dem,apophis_only
character(len=256) :: apophis_shape_file

 real :: scale_vel
 real :: scale_pos
 real :: scale_r_apophis
 real :: scale_rho

 ! Spin parameters: applied after DEM particles are placed (period=0 → no spin).
 ! Spin is DEM-only; single-sink runs are unaffected.
 real :: apophis_spin_period    ! rotation period in hours (0 = no spin)
 real :: apophis_spin_obliquity ! obliquity of spin axis from ecliptic north (degrees)
 real :: apophis_spin_azimuth   ! azimuth of spin axis in ecliptic plane (degrees)
 ! Flyby torque-alignment mode: rotate spin axis from +h to -h about Earth-Apophis separation.
 ! h = unit(r_apophis - r_earth) x unit(v_apophis - v_earth) at setup epoch.
 ! 0 deg = spin axis parallel to +h; 180 deg = parallel to -h (opposite to +h).
 ! Values < 0 disable this mode (use obliquity/azimuth instead).
 real :: apophis_spin_torque_align_deg

 private

contains
!----------------------------------------------------------------
!+
!  setup for solar system orbits
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,idust,set_particle_type,&
                        grainsize,graindens,ndustlarge,ndusttypes,ndustsmall,ihacc,igas
 use setbinary,     only:set_binary
 use units,         only:set_units,umass,udist,unit_density,unit_velocity,utime,in_code_units,in_units
 use physcon,       only:solarm,pi,au,km,solarr,ceresm,earthm,earthr,days
 use io,            only:master,fatal,warning
 use timestep,      only:tmax,dtmax
 use centreofmass,  only:reset_centreofmass
 use setsolarsystem,only:set_minor_planets,add_sun_and_planets,add_body
 use kernel,        only:hfact_default
 use eos_tillotson, only:rho_0,A
 use shape,         only:set_shape
 use options,       only:ieos
 use setup_params,  only:npart_total
 use infile_utils,  only:get_options
 use ptmass,        only:isink_potential
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout) :: massoftype(:)
 real,              intent(inout) :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: ierr,i,nerr,nptmass_dem_start
 character(len=32) :: apophis_shape_kind
 !integer :: values(8),year,month,day
 real    :: period,semia,mtot,dx
 real    :: r_apophis,m_apophis,rtidal,spsoundmin
!
! default runtime parameters
!
 tmax_in = '1000 yr'
 dtmax_in = '1 yr'
 asteroids = .true.
 np_apophis = 0
 use_dem = .false.
 apophis_only = .false.
 !call date_and_time(values=values)
 !year = values(1); month = values(2); day = values(3)
 !write(epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day
 epoch='2029-04-10'   ! encounter is on Friday 13th April
 scale_vel=1.
 scale_pos=1.
 scale_r_apophis=1.
 scale_rho=1.
 apophis_shape_file='apophis.shape'
 apophis_spin_period    = 0.
 apophis_spin_obliquity = 0.
 apophis_spin_azimuth   = 0.
 apophis_spin_torque_align_deg = -1.
!
! read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Superb Solar System Setup'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
!
! set units
!
 call set_units(mass=solarm,dist=km,G=1.d0)
!
! general parameters
!
 time  = 0.
 polyk = 0.
 gamma = 1.
 hfact = hfact_default
!
!--space available for injected gas particles
!
 npart = 0
 npart_total = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 semia  = 1.*au/udist  !  Earth
 mtot   = solarm/umass !  mass around which all bodies should orbit

 period = 2.*pi*sqrt(semia**3/mtot)
 tmax   = in_code_units(tmax_in,ierr,unit_type='time')
 if (ierr /= 0) call fatal('setup_solarsystem',' could not parse tmax')
 dtmax  = in_code_units(dtmax_in,ierr,unit_type='time')
 if (ierr /= 0) call fatal('setup_solarsystem',' could not parse dtmax')

 if (asteroids) then
    call set_minor_planets(npart,npartoftype,massoftype,xyzh,vxyzu,&
                           mtot,itype=idust,sample_orbits=.false.)
    print*,'npart = ',npart,' npartoftype = ',npartoftype(idust)
    !
    ! treat minor bodies as km-sized dust particles
    !
    ndustlarge = 1
    ndustsmall = 0
    ndusttypes = 1
    grainsize(ndustlarge) = km/udist         ! assume km-sized bodies
    graindens(ndustlarge) = 2./unit_density  ! 2 g/cm^3
 endif
 !
 ! add the planets
 !
 ierr = 0
 call add_sun_and_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
 if (nerr > 0) ierr = ierr + nerr

 if (apophis_only) then
    xyzmh_ptmass(:,1:nptmass) = 0.
    vxyz_ptmass(:,1:nptmass) = 0.
    nptmass = 0
 endif
 !
 ! add the bringer of death
 !
 if (np_apophis > 0) then
    call add_body('apophis',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
    if (nerr > 0) call warning('apophis','missing some information')

    r_apophis = xyzmh_ptmass(5,nptmass) * scale_r_apophis
    xyzmh_ptmass(5,nptmass) = r_apophis
    print "(a,1pg10.3)",' apophis radius scaled by ',scale_r_apophis

    m_apophis = 4./3.*pi*(rho_0*scale_rho/unit_density)*r_apophis**3
    xyzmh_ptmass(4,nptmass) = m_apophis
    print "(a,2(es10.3,a))",' mass of apophis is ',m_apophis*umass,&
                            ' g or ',m_apophis*umass/ceresm,' ceres masses'
    print "(a,1pg10.3,a)",' density is ',m_apophis/(4./3.*pi*r_apophis**3)*unit_density,' g/cm^3'
    print "(a,1pg10.3,a)",' Apophis-Earth relative velocity = ',&
                            in_units(sqrt(sum((vxyz_ptmass(1:3,nptmass)-vxyz_ptmass(1:3,4))**2)),'km/s'),' km/s'

    rtidal = r_apophis*(earthm/umass/m_apophis)**(1./3.)
    print "(3(a,1pg10.3),a)",' r_tidal is ',rtidal,' au,',rtidal*udist/km,' km, or ',rtidal*udist/earthr,' earth radii'

    vxyz_ptmass(1:3,nptmass) = vxyz_ptmass(1:3,nptmass)*scale_vel
    print "(a,1pg10.3)",' velocity of apophis scaled by ',scale_vel

    xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass)*scale_pos
    print "(a,1pg10.3)",' initial position of apophis scaled by ',scale_pos


    if (np_apophis > 1) then
       !
       ! replace the sink particle with a ball of stuff
       !
       call set_shape('closepacked',id,master,np_apophis,xyzmh_ptmass(1:3,nptmass),r_apophis,&
                      hfact,npart,xyzh,npart_total,objfile=apophis_shape_file,&
                      shape_kind_out=apophis_shape_kind)
       !call set_sphere('closepacked',id,master,0.,r_apophis,dx,hfact,npart,xyzh,npart_total,&
       !                xyz_origin=xyzmh_ptmass(1:3,nptmass),exactN=.true.,np_requested=np_apophis)

       do i=1,npart
          vxyzu(1:3,i) = vxyz_ptmass(1:3,nptmass)
       enddo
       massoftype(igas) = m_apophis / npart
       npartoftype(igas) = npart
       nptmass = nptmass - 1

       if (use_dem) then
          ! Record how many sinks exist before DEM replacement so we know which
          ! indices are the new Apophis rubble-pile particles afterward.
          nptmass_dem_start = nptmass
          call replace_gas_with_dem(id,npart,npartoftype(igas),massoftype(igas),&
                                    xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,hfact)
          isink_potential = 2
          ! Apply rigid-body spin after particles are placed (period=0 → skipped).
          if (apophis_spin_period > 0.) then
             ! Sphere DEM only: reorient the lattice so its I_max principal
             ! axis coincides with the spin axis before imposing spin
             ! (docs/SPHERE_LATTICE_FIX.md). OBJ/mesh shapes are untouched.
             if (trim(apophis_shape_kind) == 'sphere') then
                call align_dem_to_principal_axis(nptmass_dem_start+1,nptmass,&
                                                 xyzmh_ptmass,vxyz_ptmass,&
                                                 apophis_spin_obliquity,apophis_spin_azimuth,&
                                                 apophis_spin_torque_align_deg,i_earth=4)
             endif
             call apply_apophis_spin(nptmass_dem_start+1,nptmass,&
                                     xyzmh_ptmass,vxyz_ptmass,&
                                     apophis_spin_period,apophis_spin_obliquity,&
                                     apophis_spin_azimuth,utime,&
                                     apophis_spin_torque_align_deg,i_earth=4)
          endif
       endif
       !
       ! print quantities from the equation of state to give an idea of the timestep
       !
       if (ieos==23) then
          spsoundmin = sqrt(A/rho_0)/unit_velocity
          print "(a,1pg11.4,a)",'     sound speed min = ',spsoundmin*unit_velocity/km,' km/s'
          print "(a,1pg10.3,a)",' sound crossing time = ',(r_apophis/spsoundmin)*utime,' seconds'
       endif
    endif
 endif
 !
 ! set centre of mass as the origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (ierr /= 0) call fatal('setup','ERRORS during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  replace gas with discrete element method particles
!+
!----------------------------------------------------------------
subroutine replace_gas_with_dem(id,npart,ngas,pmass,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,hfact)
 use part, only:iReff,ihacc
 use units, only:udist
 use physcon, only:km
 use io, only:master
 integer, intent(in)    :: id
 integer, intent(inout) :: npart,ngas,nptmass
 real, intent(inout) :: pmass,xyzh(:,:),vxyzu(:,:)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(in)    :: hfact
 integer :: i,j
 real :: dxij,dyij,dzij,dmin,reff,sep(npart)
 real :: sep_med,sep_min,sep_max

 if (npart < 1) return

 !
 ! DEM sphere radius from cropped lattice geometry (not r_apophis/40).
 ! Use per-particle Reff = 0.5*nearest-neighbour distance so ellipsoid/mesh
 ! surfaces do not share one radius (median) that overlaps close pairs on step 1.
 !
 sep_med = 0.
 sep_min = huge(sep_min)
 sep_max = 0.
 if (npart == 1) then
    sep(1) = xyzh(4,1)/hfact
 else
    do i=1,npart
       dmin = huge(dmin)
       do j=1,npart
          if (j == i) cycle
          dxij = xyzh(1,i) - xyzh(1,j)
          dyij = xyzh(2,i) - xyzh(2,j)
          dzij = xyzh(3,i) - xyzh(3,j)
          dmin = min(dmin,sqrt(dxij*dxij + dyij*dyij + dzij*dzij))
       enddo
       sep(i) = dmin
       sep_min = min(sep_min,dmin)
       sep_max = max(sep_max,dmin)
    enddo
    do i=2,npart
       dmin = sep(i)
       j = i - 1
       do while (j >= 1 .and. sep(j) > dmin)
          sep(j+1) = sep(j)
          j = j - 1
       enddo
       sep(j+1) = dmin
    enddo
    if (mod(npart,2) == 1) then
       sep_med = sep((npart+1)/2)
    else
       sep_med = 0.5*(sep(npart/2) + sep(npart/2+1))
    endif
 endif

 if (id == master) then
    print "(a,1pg12.4,a)",' DEM NN spacing (min/median/max) = ',sep_min*udist/km,&
          ' / ',sep_med*udist/km,' / ',sep_max*udist/km,' km'
 endif

 xyzmh_ptmass(:,nptmass+1:) = 0.
 do i=1,npart
    nptmass = nptmass + 1
    vxyz_ptmass(1:3,nptmass) = vxyzu(1:3,i)
    xyzmh_ptmass(1:3,nptmass) = xyzh(1:3,i)
    xyzmh_ptmass(4,nptmass) = pmass
    reff = 0.5*sep(i)
    xyzmh_ptmass(iReff,nptmass) = reff
    xyzmh_ptmass(ihacc,nptmass) = reff
 enddo
 npart = 0
 ngas = 0
 pmass = 0.

end subroutine replace_gas_with_dem

!----------------------------------------------------------------
!+
!  Compute the prescribed Apophis spin axis unit vector, either from
!  the flyby torque-alignment parameterisation (rotate the
!  Earth-Apophis orbital angular momentum direction about the
!  separation vector by torque_align_deg) or from ecliptic
!  obliquity/azimuth. ierr=0 on success; ierr/=0 means the caller must
!  not apply any spin (mirrors the original apply_apophis_spin abort
!  conditions: torque-align requested but no Earth sink index given,
!  or the given index is out of range). verbose (default true)
!  controls whether the diagnostic h_hat/r_hat/torque-align lines are
!  printed, so callers that need the axis without re-printing
!  diagnostics (e.g. the lattice pre-alignment step) can pass .false.
!+
!----------------------------------------------------------------
subroutine compute_apophis_spin_axis(i_start,i_end,xyzmh_ptmass,vxyz_ptmass,&
                                      obliquity_deg,azimuth_deg,torque_align_deg,&
                                      nx,ny,nz,ierr,i_earth,verbose)
 use physcon, only:pi
 integer, intent(in)  :: i_start,i_end
 real,    intent(in)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)  :: obliquity_deg,azimuth_deg,torque_align_deg
 real,    intent(out) :: nx,ny,nz
 integer, intent(out) :: ierr
 integer, intent(in), optional :: i_earth
 logical, intent(in), optional :: verbose
 integer :: i,n,ie
 real    :: obl_rad,az_rad
 real    :: rcm(3),v_spin(3)
 real    :: rrel(3),vrel(3),hvec(3),rhat(3),hnorm,hn(3)
 real    :: theta,ct,st,kdotv
 logical :: use_torque_align,say

 say = .true.
 if (present(verbose)) say = verbose
 ierr = 0
 use_torque_align = .false.
 nx = 0.; ny = 0.; nz = 1.
 n = i_end - i_start + 1

 if (torque_align_deg >= 0.) then
    if (.not.present(i_earth)) then
       if (say) print "(a)",' ERROR: apophis_spin_torque_align_deg requires Earth sink index'
       ierr = 1
       return
    endif
    ie = i_earth
    if (ie < 1 .or. ie > size(xyzmh_ptmass,2)) then
       ierr = 1
       return
    endif
    ! Apophis DEM centre-of-mass position and velocity at setup epoch.
    rcm = 0.
    v_spin = 0.
    do i = i_start, i_end
       rcm(1:3) = rcm(1:3) + xyzmh_ptmass(1:3,i)
       v_spin(1:3) = v_spin(1:3) + vxyz_ptmass(1:3,i)
    enddo
    rcm = rcm / real(n)
    v_spin = v_spin / real(n)
    rrel = rcm - xyzmh_ptmass(1:3,ie)
    vrel = v_spin - vxyz_ptmass(1:3,ie)
    hvec = (/ rrel(2)*vrel(3) - rrel(3)*vrel(2), &
              rrel(3)*vrel(1) - rrel(1)*vrel(3), &
              rrel(1)*vrel(2) - rrel(2)*vrel(1) /)
    hnorm = sqrt(sum(rrel**2))
    if (hnorm <= 0.) then
       if (say) print "(a)",' WARN: zero Earth-Apophis separation; using ecliptic spin axis'
    else
       rhat = rrel / hnorm
       hnorm = sqrt(sum(hvec**2))
       if (hnorm <= 0.) then
          if (say) print "(a)",' WARN: collinear Earth-Apophis r,v; using ecliptic spin axis'
       else
          hn = hvec / hnorm
          ! Rodrigues rotation of hn about rhat by torque_align_deg (0=+h, 180=-h).
          theta = torque_align_deg * pi / 180.
          ct = cos(theta)
          st = sin(theta)
          kdotv = rhat(1)*hn(1) + rhat(2)*hn(2) + rhat(3)*hn(3)
          nx = hn(1)*ct + (rhat(2)*hn(3) - rhat(3)*hn(2))*st + rhat(1)*kdotv*(1.-ct)
          ny = hn(2)*ct + (rhat(3)*hn(1) - rhat(1)*hn(3))*st + rhat(2)*kdotv*(1.-ct)
          nz = hn(3)*ct + (rhat(1)*hn(2) - rhat(2)*hn(1))*st + rhat(3)*kdotv*(1.-ct)
          if (say) then
             print "(a,1pg10.3,a)",' Apophis spin torque-align = ',torque_align_deg,' deg (0=+h, 180=-h)'
             print "(a,3(1pg10.3,1x))",' Earth-Apophis h_hat (orbit normal) = ',hn(1),hn(2),hn(3)
             print "(a,3(1pg10.3,1x))",' Earth-Apophis r_hat (separation) = ',rhat(1),rhat(2),rhat(3)
          endif
          use_torque_align = .true.
       endif
    endif
 endif

 if (.not.use_torque_align) then
    obl_rad = obliquity_deg * pi / 180.
    az_rad  = azimuth_deg   * pi / 180.
    ! Spin axis unit vector from obliquity and azimuth (ecliptic frame)
    nx = sin(obl_rad) * cos(az_rad)
    ny = sin(obl_rad) * sin(az_rad)
    nz = cos(obl_rad)
 endif

end subroutine compute_apophis_spin_axis

!----------------------------------------------------------------
!+
!  Reorient a closepacked-lattice DEM sphere about its centre of mass
!  so its maximum-inertia principal axis aligns with the spin axis
!  apply_apophis_spin is about to impose. A finite HCP lattice cropped
!  to a sphere is not exactly isotropic (~2% spread in principal
!  moments, fixed in the simulation x/y/z frame regardless of the
!  requested spin axis); without this step, sweeping the spin axis
!  spins the same slightly-triaxial body about different axes
!  relative to its own principal axes, biasing the intrinsic spin
!  period estimate (docs/METRICS.md, docs/SPHERE_LATTICE_FIX.md).
!  No-op if fewer than 2 grains, if the axis computation aborts
!  (mirrors compute_apophis_spin_axis's ierr), or if already aligned.
!+
!----------------------------------------------------------------
subroutine align_dem_to_principal_axis(i_start,i_end,xyzmh_ptmass,vxyz_ptmass,&
                                        obliquity_deg,azimuth_deg,torque_align_deg,i_earth)
 use vectorutils, only:jacobi_eigen_sym,rotatevec,rotation_to_align
 use physcon,     only:pi
 integer, intent(in)    :: i_start,i_end
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 real,    intent(in)    :: vxyz_ptmass(:,:)
 real,    intent(in)    :: obliquity_deg,azimuth_deg,torque_align_deg
 integer, intent(in), optional :: i_earth
 integer :: i,n,nrot,ierr,imax
 real    :: rcm(3),dr(3)
 real    :: inertia(3,3),evec(3,3),eval(3)
 real    :: nx,ny,nz,ntarget(3),emax(3)
 real    :: rot_axis(3),rot_angle
 logical :: already_aligned

 n = i_end - i_start + 1
 if (n < 2) return

 ! Target spin axis: identical computation apply_apophis_spin will use
 ! (verbose=.false. so setup.log doesn't get the h_hat/r_hat lines twice).
 call compute_apophis_spin_axis(i_start,i_end,xyzmh_ptmass,vxyz_ptmass,&
                                 obliquity_deg,azimuth_deg,torque_align_deg,&
                                 nx,ny,nz,ierr,i_earth=i_earth,verbose=.false.)
 if (ierr /= 0) return
 ntarget = (/nx,ny,nz/)

 ! Centre of mass of DEM grains
 rcm = 0.
 do i = i_start, i_end
    rcm(1:3) = rcm(1:3) + xyzmh_ptmass(1:3,i)
 enddo
 rcm = rcm / real(n)

 ! Equal-mass grain inertia tensor about the CoM (mass factor cancels
 ! in the eigenvectors, so it is omitted; only directions are needed)
 inertia = 0.
 do i = i_start, i_end
    dr = xyzmh_ptmass(1:3,i) - rcm
    inertia(1,1) = inertia(1,1) + dr(2)**2 + dr(3)**2
    inertia(2,2) = inertia(2,2) + dr(1)**2 + dr(3)**2
    inertia(3,3) = inertia(3,3) + dr(1)**2 + dr(2)**2
    inertia(1,2) = inertia(1,2) - dr(1)*dr(2)
    inertia(1,3) = inertia(1,3) - dr(1)*dr(3)
    inertia(2,3) = inertia(2,3) - dr(2)*dr(3)
 enddo
 inertia(2,1) = inertia(1,2)
 inertia(3,1) = inertia(1,3)
 inertia(3,2) = inertia(2,3)

 call jacobi_eigen_sym(inertia,3,3,eval,evec,nrot)
 imax = maxloc(eval,dim=1)
 emax = evec(:,imax)

 ! Rotation that maps emax onto ntarget (Rodrigues, about their cross product;
 ! handles the already-aligned and exactly-antiparallel degenerate cases).
 call rotation_to_align(emax,ntarget,rot_axis,rot_angle,already_aligned)
 if (already_aligned) then
    print "(a)",' Apophis principal-axis alignment: already aligned (no rotation applied)'
    return
 endif

 do i = i_start, i_end
    dr = xyzmh_ptmass(1:3,i) - rcm
    call rotatevec(dr,rot_axis,rot_angle)
    xyzmh_ptmass(1:3,i) = rcm + dr
 enddo

 print "(a,3(1pg10.3,1x))",' Apophis principal moments of inertia = ',eval(1),eval(2),eval(3)
 print "(a,1pg10.3,a)",' Apophis principal-axis alignment: rotated grains by ',&
       rot_angle*180./pi,' deg so I_max axis || prescribed spin axis'

end subroutine align_dem_to_principal_axis

!----------------------------------------------------------------
!+
!  Apply rigid-body spin to DEM sink particles after placement.
!  Computes the angular velocity from the requested spin period,
!  then adds v_spin = omega_vec x (r_i - r_cm) to each particle so
!  the rubble pile enters the simulation already rotating.
!  No-op when period_hr <= 0 or the index range is empty.
!+
!----------------------------------------------------------------
subroutine apply_apophis_spin(i_start,i_end,xyzmh_ptmass,vxyz_ptmass,&
                               period_hr,obliquity_deg,azimuth_deg,utime,&
                               torque_align_deg,i_earth)
 use physcon, only:pi
 integer, intent(in)    :: i_start,i_end
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)    :: period_hr,obliquity_deg,azimuth_deg,utime
 real,    intent(in)    :: torque_align_deg
 integer, intent(in), optional :: i_earth
 integer :: i,n,ierr
 real    :: omega,nx,ny,nz
 real    :: rcm(3),dr(3),v_spin(3)

 n = i_end - i_start + 1
 if (n < 1 .or. period_hr <= 0.) return

 ! Angular velocity in code units: omega = 2*pi / T, T = period_hr * 3600 s / utime
 omega = 2.*pi / (period_hr * 3600. / utime)

 call compute_apophis_spin_axis(i_start,i_end,xyzmh_ptmass,vxyz_ptmass,&
                                 obliquity_deg,azimuth_deg,torque_align_deg,&
                                 nx,ny,nz,ierr,i_earth=i_earth)
 if (ierr /= 0) return

 ! Centre of mass of DEM particles
 rcm = 0.
 do i = i_start, i_end
    rcm(1:3) = rcm(1:3) + xyzmh_ptmass(1:3,i)
 enddo
 rcm = rcm / real(n)

 ! Add omega_vec x (r_i - r_cm) to each particle's translational velocity
 do i = i_start, i_end
    dr(1) = xyzmh_ptmass(1,i) - rcm(1)
    dr(2) = xyzmh_ptmass(2,i) - rcm(2)
    dr(3) = xyzmh_ptmass(3,i) - rcm(3)
    v_spin(1) = omega * (ny*dr(3) - nz*dr(2))
    v_spin(2) = omega * (nz*dr(1) - nx*dr(3))
    v_spin(3) = omega * (nx*dr(2) - ny*dr(1))
    vxyz_ptmass(1:3,i) = vxyz_ptmass(1:3,i) + v_spin(1:3)
 enddo

 print "(a,1pg10.3,a)",' Apophis spin period    = ',period_hr,' hr'
 print "(a,3(1pg10.3,1x))",' Apophis spin axis (nx,ny,nz) = ',nx,ny,nz
 print "(a,1pg10.3,a)",' Apophis spin omega     = ',omega/utime,' rad/s'

end subroutine apply_apophis_spin

!----------------------------------------------------------------
!+
!  write setup parameters to file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for solar system setup routines'
 call write_inopt(tmax_in,'tmax_in','end time of simulation (e.g. 3 days)',iunit)
 call write_inopt(dtmax_in,'dtmax_in','time between dumps (e.g. 1 hr)',iunit)
 call write_inopt(asteroids,'asteroids','add distant minor bodies as km-sized dust particles',iunit)
 call write_inopt(np_apophis,'np_apophis','number of particles used to represent apophis (0=none; 1=sink; n=gas)',iunit)
 call write_inopt(epoch,'epoch','epoch to query ephemeris, YYYY-MMM-DD HH:MM:SS.fff, blank = today',iunit)

 call write_inopt(use_dem,'use_dem','use the discrete element method for sink-sink interactions',iunit)
 call write_inopt(apophis_only,'apophis_only','only add apophis',iunit)

 call write_inopt(scale_vel,'scale_vel','scaling factor for apophis velocity',iunit)
 call write_inopt(scale_pos,'scale_pos','scaling factor for apophis initial position',iunit)
 call write_inopt(scale_r_apophis,'scale_r_apophis','scaling factor for apophis radius',iunit)
 call write_inopt(scale_rho,'scale_rho','scaling factor for apophis bulk density',iunit)
call write_inopt(apophis_shape_file,'apophis_shape_file','shape config file for lattice cropping',iunit)
 call write_inopt(apophis_spin_period,   'apophis_spin_period',   'Apophis spin period in hours (0 = no spin, DEM only)',iunit)
 call write_inopt(apophis_spin_obliquity,'apophis_spin_obliquity','obliquity of spin axis from ecliptic north (degrees)',iunit)
 call write_inopt(apophis_spin_azimuth,  'apophis_spin_azimuth',  'azimuth of spin axis in ecliptic plane (degrees)',iunit)
 call write_inopt(apophis_spin_torque_align_deg,'apophis_spin_torque_align_deg',&
      'flyby torque alignment: rotate spin from +h to -h about Earth-Apophis r_hat; <0=use obl/az',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read setup parameters from file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(tmax_in, 'tmax_in',db,errcount=nerr)
 call read_inopt(dtmax_in,'dtmax_in',db,errcount=nerr)
 call read_inopt(asteroids,'asteroids',db,errcount=nerr)
 call read_inopt(np_apophis,'np_apophis',db,min=0,errcount=nerr)
 call read_inopt(epoch,'epoch',db,errcount=nerr)

 call read_inopt(use_dem,'use_dem',db,errcount=nerr)
 call read_inopt(apophis_only,'apophis_only',db,errcount=nerr)

 call read_inopt(scale_vel,'scale_vel',db,default=1.0,errcount=nerr)
 call read_inopt(scale_pos,'scale_pos',db,default=1.0,errcount=nerr)
 call read_inopt(scale_r_apophis,'scale_r_apophis',db,default=1.0,errcount=nerr)
 call read_inopt(scale_rho,'scale_rho',db,default=1.0,errcount=nerr)
call read_inopt(apophis_shape_file,'apophis_shape_file',db,default='apophis.shape',errcount=nerr)
 call read_inopt(apophis_spin_period,   'apophis_spin_period',   db,default=0.0,errcount=nerr)
 call read_inopt(apophis_spin_obliquity,'apophis_spin_obliquity',db,default=0.0,errcount=nerr)
 call read_inopt(apophis_spin_azimuth,  'apophis_spin_azimuth',  db,default=0.0,errcount=nerr)
 call read_inopt(apophis_spin_torque_align_deg,'apophis_spin_torque_align_deg',db,default=-1.0,errcount=nerr)

 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
