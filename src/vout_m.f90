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
! calculates Hartree and exchange-correlation potentials and total-energy
! term to add to eigenvalue sum

 subroutine vout_m(mode,rho,tau,vh,vxc,vtau,zz,sf,eeel,eexc,etot, &
&                rr,mmax,iexc)

!mode  0 for full calculation, 1 for xc only (eg., with rhomod)
!rho  total charge density or valence/pseudovalence charge density
!tau  kinetic energy density
!vo  hartree potential
!vxc  output exchange-correlation potential
!vtau  dexc/dtau for meta-gga
!zz  atomic number (or pseudo)
!sf  electron number
!eeel  electron-electron interaction energy
!eexc  exchange-correlation total energy
!rr  log radial mesh
!output electrostatic and exchange-correlation potential

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0_dp*pi

!Input vaiables
 integer :: iexc,mmax,mode
 real(dp) :: rho(mmax),tau(mmax),rr(mmax)
 real(dp) :: zz,sf

!Output variables
 real(dp) :: vh(mmax),vxc(mmax),vtau(mmax)
 real(dp) :: eeel,eexc,etot

!Local function
 real(dp) :: tfapot,fcrossover

!Local variables
 integer ii,nsm1,nsm2
 real(dp) :: al,tv,xcr

!Function
 real(dp) :: aii

 real(dp), allocatable :: rvp(:),rv(:)
 real(dp), allocatable :: vxcd(:),exca(:)
 real(dp), allocatable :: rhosm(:),tausm(:),vxcsm(:),vtausm(:)
 real(dp), allocatable :: foutv(:),foutx(:),foutfc(:)

 allocate(rvp(mmax),rv(mmax))
 allocate(vxcd(mmax),exca(mmax))
 allocate(rhosm(mmax),tausm(mmax),vxcsm(mmax),vtausm(mmax))
 allocate(foutv(mmax),foutx(mmax),foutfc(mmax))


 al = 0.01d0 * dlog(rr(101) / rr(1))

 foutfc(:)=0.0d0

 vh(:)=0.0d0 ; vxc(:)=0.0d0 ; vtau(:)=0.0d0 ; eeel=0.0d0

 if(mode==0) then

! integration for electrostatic potential for the
! electron-electron interaction for the total energy

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**3
 end do

 rv(mmax)=sf
 rv(mmax-1)=sf
 rv(mmax-2)=sf

 do ii=mmax-2,2,-1
   rv(ii-1)=rv(ii)+aii(rvp,ii)
 end do

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**2
 end do

 tv=0.0d0
 do ii=mmax-2,2,-1
   tv=tv+aii(rvp,ii)
   rv(ii-1)=rv(ii-1)-rr(ii-1)*tv
 end do

 do ii=1,mmax
   vh(ii)=rv(ii)/rr(ii)
 end do
 do ii=1,mmax
  foutv(ii)=rho(ii)*vh(ii)*rr(ii)**3
 end do
 eeel=(9.0d0*foutv(1) + 28.0d0*foutv(2) &
&   + 23.0d0*foutv(3))/24.0d0
 do ii=4,mmax
   eeel=eeel + foutv(ii)
 end do
 eeel=al*eeel + foutv(1)/3.0d0

! integration for the hartree potential
 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**3
 end do

 rv(mmax)=zz-sf
 rv(mmax-1)=zz-sf
 rv(mmax-2)=zz-sf

 do ii=mmax-2,2,-1
   rv(ii-1)=rv(ii)+aii(rvp,ii)
 end do

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**2
 end do

 tv=0.0d0
 do ii=mmax-2,2,-1
   tv=tv+aii(rvp,ii)
   rv(ii-1)=rv(ii-1)-rr(ii-1)*tv
 end do

 do ii=1,mmax
   vh(ii)=rv(ii)/rr(ii)
 end do

 end if !mode==0

! exchange-correlation potential added

 vtau(:)=0.0d0
 vxc(:)=0.0d0
 if(iexc .eq. 1) then
   call excwig(rho,vxc,exca,mmax)
 else if(iexc .eq. 2) then
   call exchdl(rho,vxc,exca,mmax)
 else if(iexc .eq. 3) then
   call excpzca(rho,vxc,exca,mmax)
 else if(iexc .eq. 4) then
   call excggc(rho,vxc,exca,rr,mmax)

 else if(iexc .eq. 5) then

   call excr2scan01_m(al,rho,tau,vxc,vtau,exca,rr,mmax)

!if(.false.) then

! smoothing cutoffs
 do ii=1,mmax
!  if(rr(ii) > 1.0d2*rr(1) .or. rr(ii) >5.0d-3) then
   if(rr(ii) > 1.0d1*rr(1) .or. rr(ii) >5.0d-4) then
     nsm1=ii
     exit
   end if
 end do

 do ii=1,mmax
   if(rr(ii) > 10.0d0*rr(nsm1)) then
     nsm2=ii
     exit
   end if
 end do

   call gauss_smooth(rr,vxc,vxcsm,nsm2+20,25)
   call gauss_smooth(rr,vtau,vtausm,nsm2+20,25)


  do ii=1,mmax

   xcr=2.0d0*(rr(ii)-rr(nsm1))/(rr(nsm2)-rr(nsm1)) - 1.0d0

   vtau(ii)=fcrossover(xcr)*vtau(ii) + (1.0d0-fcrossover(xcr))*vtausm(ii)

   vxc(ii)=fcrossover(xcr)*vxc(ii) + (1.0d0-fcrossover(xcr))*vxcsm(ii)
  end do
!end if !.false.

!else if(iexc <0) then
!  call exc_libxc_m(iexc,al,rho,tau,vxc,vtau,exca,rr,mmax)
 else
   write(6,'(/a,i7)') 'ERROR vout_m: bad input iexc =',iexc
   stop
 end if

! exchange-correlation total energy

 do ii=1,mmax
  foutx(ii)=rho(ii)*exca(ii)*rr(ii)**3
 end do
 eexc=(9.0d0*foutx(1) + 28.0d0*foutx(2) &
&   + 23.0d0*foutx(3))/24.0d0
 do ii=4,mmax
   eexc=eexc + foutx(ii)
 end do
!eexc=al*eexc + foutx(1)/3.0d0
 eexc=al*eexc + foutx(1)/2.0d0

 do ii=1,mmax
  foutx(ii)=(tau(ii) + rho(ii)*(vh(ii) + exca(ii)))*rr(ii)**3
! foutx(ii)=(tau(ii) + rho(ii)*(vh(ii) + vxc(ii)))*rr(ii)**3
 end do
 etot=(9.0d0*foutx(1) + 28.0d0*foutx(2) &
&   + 23.0d0*foutx(3))/24.0d0
 do ii=4,mmax
   etot=etot + foutx(ii)
 end do
!etot=al*etot + foutx(1)/3.0d0
 etot=al*etot + foutx(1)/2.0d0

! total energy output based on tau integral, not eigenvalues
! this is required for the metagga; the two are not the same
 etot=etot - 0.5d0*eeel


 deallocate(rvp,rv)
 deallocate(vxcd,exca)
 deallocate(rhosm,tausm,vxcsm,vtausm)
 deallocate(foutv,foutx,foutfc)
 return

 end subroutine vout_m
