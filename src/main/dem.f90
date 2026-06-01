!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dem
!
! This module implements the soft sphere discrete element method (DEM)
! for sink-sink interactions.
!
! epsilon is a parameter [0,1] that controls the restitution of the collision.
! epsilon = 1 is perfectly elastic, epsilon = 0 is perfectly inelastic.
! epsilon = 0.5 is a reasonable default.
!
! :References: Schwartz+2012, Granular Matter 14, 363-380
!
! :Owner: Daniel Price
!
 implicit none
 private

 public :: get_ssdem_force

 real, public :: ct_dem = 0.0         ! Tangential damping coefficient (0=frictionless; dashpot-only model is unphysical for rotating bodies)
 real, public :: epsilon_n_dem = 0.5  ! Normal coefficient of restitution (user-settable)
 real, public :: kn_cgs = 1e7         ! Spring constant (e.g. 10^4 kg/s^2 = 10^7 g/s^2)
 real, public :: kc_cgs = 0.0         ! Cohesive spring constant (dyne/cm); 0 = no cohesion
 real, public :: dn_cohes_factor = 0.1 ! Cohesion range as fraction of combined radii

contains

!----------------------------------------------------------------
!+
!  Soft-sphere DEM normal force (Hooke's law)
!  Implements Eq. (3) from Schwartz+2012 for overlapping spheres
!+
!----------------------------------------------------------------
subroutine get_ssdem_force(Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,fx,fy,fz,veli,velj,wi,wj,dtmin)
 use physcon,     only:pi
 use vectorutils, only:cross_product
 use units,       only:umass,utime
 real, intent(in)    :: Rsinki,Rsinkj,mi,mj,ddr,dx,dy,dz,veli(3),velj(3),wi(3),wj(3)
 real, intent(inout) :: fx,fy,fz,dtmin
 real :: r,overlap,kn,kn_dem
 real :: cn,ct,reduced_mass,log_epsilon_n_dem,li,lj
 real :: nvec(3),vrel(3),n_cross_wi(3),n_cross_wj(3),u_dot_n,u_n(3),u_t(3)
 real :: gap,kc_dem_val
 logical, save :: cohes_print_done = .false.

 !----------------------------------------------------------------
 ! Normal force
 !----------------------------------------------------------------
 r = 1.0 / ddr
  ! Normal unit vector
 nvec(1) = dx * ddr
 nvec(2) = dy * ddr
 nvec(3) = dz * ddr

 overlap = Rsinki + Rsinkj - r
 kn = 0.
 kn_dem = kn_cgs / (umass/utime**2)  ! convert to code units
 if (overlap > 0.0) then
    ! Spring force (Schwartz+2012 Eq. 3)
    kn = kn_dem
    fx = fx + kn * overlap * nvec(1) / mj
    fy = fy + kn * overlap * nvec(2) / mj
    fz = fz + kn * overlap * nvec(3) / mj

    !----------------------------------------------------------------
    ! Damping (only during contact: Schwartz+2012 Eqs. 8-15)
    ! Applying ct to non-contacting spinning pairs creates O(omega*dx)
    ! relative velocities on all pairs, producing enormous artificial
    ! forces that crush the timestep.
    !----------------------------------------------------------------
    n_cross_wi = cross_product(nvec,wi)
    n_cross_wj = cross_product(nvec,wj)
    li = (Rsinki**2 - Rsinkj**2 + r**2) / (2.0 * r)
    lj = (Rsinkj**2 - Rsinki**2 + r**2) / (2.0 * r)
    vrel = veli - velj + li * n_cross_wi - lj * n_cross_wj
    u_dot_n = dot_product(vrel, nvec)
    u_n = u_dot_n * nvec
    u_t = vrel - u_n
    reduced_mass = mj * mi / (mj + mi)
    log_epsilon_n_dem = log(epsilon_n_dem)
    cn = -2.0 * sqrt(reduced_mass * kn) * log_epsilon_n_dem / sqrt(pi**2 + log_epsilon_n_dem**2)
    ct = ct_dem
    fx = fx - cn * u_n(1) / mj - ct * u_t(1)
    fy = fy - cn * u_n(2) / mj - ct * u_t(2)
    fz = fz - cn * u_n(3) / mj - ct * u_t(3)

    ! Spring timescale: only relevant when spring is active
    dtmin = min(dtmin,sqrt(reduced_mass/kn_dem))
 endif

 ! Cohesive attraction: acts when particles are separated but within dn_cohes_factor*(Ri+Rj)
 if (overlap < 0. .and. kc_cgs > 0.) then
    gap = -overlap
    if (gap < dn_cohes_factor * (Rsinki + Rsinkj)) then
       kc_dem_val = kc_cgs / (umass/utime**2)
       fx = fx - kc_dem_val * gap * nvec(1) / mj
       fy = fy - kc_dem_val * gap * nvec(2) / mj
       fz = fz - kc_dem_val * gap * nvec(3) / mj
       if (.not. cohes_print_done) then
          print*, '[DEM cohesion] ACTIVE: kc_cgs=', kc_cgs, &
                  ' F_cohes=', kc_dem_val*gap, ' (code units/mj)'
          cohes_print_done = .true.
       endif
    endif
 endif

end subroutine get_ssdem_force

end module dem