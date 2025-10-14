!
! Copyright (c) 1989-2024 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine run_phsft_m(lmax,lloc,nproj,epa,epsh1,epsh2,depsh,vkb,evkb, &
&                     rr,vfull,vtau,vp,vtaups,zz,mmax,mxprj,irc,srel,ircmin)

! computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
! at rr(irphs) comparing all-electron with Vanderbilt-Kleinman-Bylander
! results for 1 and 2 projectors, or the semi-local pseudpotential
! when that is the local potential for some l
! the computed quantity is reminiscent of a scattering phase shift, but isn't

!lmax  maximum angular momentum
!lloc  l for local potential
!nproj  number ov V / KB projectors for  each l
!ep  bound-state or scattering state reference energies for vkb potentials
!epsh1  low energy limit for "phase shift" calculation
!epsh2  high energy limit for "phase shift" calculation
!depsh  energy increment
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!rr  log radial grid
!vfull  all-electron potential
!vtau  meta-gga dExc/d(KE density) all-electron ground state
!vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
!vtaups  meta-gga dExc/d(KE density) pseudo ground state
!zz  atomic number
!mmax  size of radial grid
!mxprj dimension of number of projectors
!irphs  index of rr beyond which all vp==vlocal
!srel .true. for scalar-relativistic, .false. for non-relativistic

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,mmax,mxprj,ircmin
 integer :: nproj(6),irc(6)
 real(dp) :: epsh1,epsh2,depsh,zz
 real(dp) :: rr(mmax),vp(mmax,5),epa(mxprj,6)
 real(dp) :: vfull(mmax),vkb(mmax,mxprj,4),evkb(mxprj,4)
 real(dp) :: vtau(mmax), vtaups(mmax)
 logical :: srel

!Output variables - printing only

!Local variables
 integer :: ii,irphs,ll,l1,npsh
 real(dp) :: epsh

 real(dp),allocatable :: pshf(:),pshp(:)

 npsh=int(((epsh2-epsh1)/depsh)-0.5d0)+1

 allocate(pshf(npsh),pshp(npsh))

! loop for phase shift calculation -- full, then Kleinman-
! Bylander / Vanderbilt
 
 do l1 = 1, 4

   ll = l1 - 1
   if(ll<=lmax) then
     irphs=irc(l1)+2
   else
     irphs=irc(lloc+1)
   end if

   call fphsft_m(ll,epsh2,depsh,pshf,rr,vfull,vtau,zz,mmax,irphs,npsh,srel,ircmin)

   call  vkbphsft_m(ll,nproj(l1),epsh2,depsh,epa(1,l1),pshf,pshp, &
&                 rr,vp(1,lloc+1),vtaups,vkb(1,1,l1),evkb(1,l1), &
&                 mmax,irphs,npsh)

   write(6,'(/a,i2)') 'log derivativve data for plotting, l=',ll
   write(6,'(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=',rr(irphs)
   write(6,'(a/)') 'l, energy, all-electron, pseudopotential'
   do ii = 1, npsh
     epsh = epsh2 - depsh * dfloat(ii - 1)
     write(6,'(a, i6, 3f12.6)') '! ',ll, epsh, pshf(ii), pshp(ii)
   end do
 end do
 deallocate(pshf,pshp)
 return
 end subroutine run_phsft_m
