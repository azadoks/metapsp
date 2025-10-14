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
 subroutine sratom_m(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&           rr,vi,vh,vxc,tauc,tau,vitau,zz,mmax,iexc, &
&           etotk,dltetot,ierr,srel,ircmin)

! self-consistent scalar-relativistic all-electron atom
! calculation using log mesh (non-relativistic when srel=.false.)

!na  principal quantum number array, dimension ncv
!la  angular-momenta
!ea  eigenvalues (output)
!fa  occupancies
!rpk  radius of outermost peak of wave function
!nc  number of core states
!ncv  number of core+valence states
!it  number of iterations (output)
!rr  log radial mesh
!vi  all-electron potential (in/out)
!zz  atomic number
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!etot  all-electron total energy (output)
!dltetot  final etot - first iteration etot
!ierr  error flag
!srel  .true. for scalar-relativistic, .false. for non-relativistic
!ircmin  minimum core radii index

 implicit none
 integer, parameter :: dp=kind(1.0d0)

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4

!Input variables
 
 integer :: mmax,iexc,ircmin,nc,ncv
 integer :: na(30),la(30)
 real(dp) :: zz
 real(dp) :: fa(30),rr(mmax)
 logical :: srel

!Output variables
 integer :: it,ierr
 real(dp) :: etotk,dltetot
 real(dp) :: ea(30),rpk(30),eatst(30,2)
 real(dp) :: rho(mmax),vi(mmax),tau(mmax),vitau(mmax)
 real(dp) :: rhoc(mmax),tauc(mmax)
 real(dp) :: vh(mmax),vxc(mmax)
 real(dp) :: perf(3)

!Local function
 real(dp) :: tfapot,fcrossover

!Local variables
 integer :: mch
 real(dp) :: amesh,al
 real(dp) :: dr,eeel,eexc,et,rl,rl1,sd,sf,sn,st,eeig,etot0,etot1,etot,detot
 real(dp) :: thl,vn,vntau,zion
 real(dp) :: sls,rcr,xcr
 integer :: ii,jj,nsm1,nsm2
 logical :: convg

 real(dp), allocatable :: u(:),up(:),u0(:,:),u0p(:,:)
 real(dp), allocatable :: vo(:),vi1(:),vo1(:)
 real(dp), allocatable :: work(:)
 real(dp), allocatable :: rh(:), rhgrd(:),rhpp(:),rhlap(:)
 real(dp), allocatable :: votau(:),vi1tau(:),vo1tau(:)
 real(dp), allocatable :: dpr(:), dppr(:), dlap(:), dvtau(:), hu(:)
 real(dp), allocatable :: vxcsm(:),votausm(:)

!real(dp), parameter ::  bl=0.5d0 ! blend parameter
!real(dp), parameter ::  bl=0.4d0 ! blend parameter
 real(dp), parameter ::  bl=0.3d0 ! blend parameter
 real(dp), parameter ::  eps=1.0d-12 ! etot convergence scale

 allocate(u(mmax),up(mmax),u0(mmax,ncv),u0p(mmax,ncv))
 allocate(vo(mmax),vi1(mmax),vo1(mmax))
 allocate (votau(mmax),vi1tau(mmax),vo1tau(mmax))
 allocate (vxcsm(mmax),votausm(mmax))
 allocate(work(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

 vitau(:)=0.0d0 ; vi1tau(:)=0.0d0 ; votau(:)=0.0d0 ; vo1tau=0.0d0
 etot1=0.0d0

! starting approximation for energies
 sf=0.0d0
 do jj=1,ncv
   sf=sf+fa(jj)
   zion=zz+1.0d0-sf
!  ea(jj)=-0.5d0*(zion/na(jj))**2
!  if(ea(jj)>vi(mmax)) ea(jj)=2.0d0*vi(mmax)
 end do !jj

! smoothing cutoffs
 do ii=1,mmax
   if(rr(ii) > 1.0d2*rr(1) .or. rr(ii) >1.0d-2) then
     nsm1=ii
     exit
   end if
 end do

 do ii=1,mmax
   if(rr(ii) > 1.0d-1) then
!  if(rr(ii) > 5.0d-2) then
     nsm2=ii
     exit
   end if
 end do

! big self  self-consietency loop

 do it=1,100
!do it=1,600
   
   convg=.true.

   rho(:)=0.0d0
   rhoc(:)=0.0d0
   tau(:)=0.0d0
   tauc(:)=0.0d0

! solve for bound states in turn
   eeig=0.0d0
   do jj=1,ncv
     et=ea(jj)
     ierr = 0
     call lschfb_m(na(jj),la(jj),ierr,et, &
&                  rr,vi,vitau,u,up,zz,mmax,mch,srel,ircmin)
     if(ierr/=0) then
       write(6,'(/a,4i4,1p,e18.6)') 'sratom_m: lschfb_m convergence error n,l,iter,ierr', &
&       na(jj),la(jj),it,ierr,et
       stop
!      exit
     end if
     ea(jj)=et

! accumulate charge and eigenvalues
     eeig = eeig + fa(jj) * et
     do ii=1,mmax
       if(abs(u(ii))<1.0d-145) exit
       rho(ii)=rho(ii) + fa(jj)*(u(ii)/rr(ii))**2
       sls=la(jj)*(la(jj)+1)
       tau(ii)=tau(ii) + 0.5d0*fa(jj)*(((up(ii)/al - u(ii))/rr(ii)**2)**2 &
&          + (sls/rr(ii)**2)*(u(ii)/rr(ii))**2)

       if(jj<=nc) then
         rhoc(ii)=rhoc(ii) + fa(jj)*(u(ii)/rr(ii))**2
         tauc(ii)=tauc(ii) + 0.5d0*fa(jj)*(((up(ii)/al - u(ii))/rr(ii)**2)**2 &
&            + (sls/rr(ii)**2)*(u(ii)/rr(ii))**2)
       end if
     end do

! find outermost peak of wavefunction
     do ii=mch-1,1,-1
       if(up(ii)*up(ii+1)<0.0d0) then
         rpk(jj)=rr(ii)
         exit
       end if
     end do

! save full set of wave functions
     u0(:,jj)=u(:)
     u0p(:,jj)=up(:)

   end do !jj

   if(ierr/=0) then
    exit
   end if

! output potential
   call  vout_m(0,rho,tau,vh,vxc,votau,zz,sf,eeel,eexc,etotk, &
                rr,mmax,iexc)


vo(:)=vh(:) + vxc(:)
!  charge times input potential for kinetic energy
   work(:)=0.0d0 
   do ii=1,mmax
    work(ii)=rho(ii)*(vh(ii)-vi(ii))*rr(ii)**3
   end do
   sn=(9.0d0*work(1) + 28.0d0*work(2) &
      + 23.0d0*work(3))/24.0d0
   do ii=4,mmax
     sn=sn + work(ii)
   end do
   sn=al*sn + work(1)/3.0d0

!etot from eeig calculation better for convergence testing

!  etot =  eeig + sn + eexc - 0.5d0*eeel
   etot =  eeig + sn + eexc
   detot=etot-etot1


! test for total energy converrgence
  if(abs(detot)>eps*abs(etot)) convg=.false.

   etot1=etot
   if(it==1) etot0=etot

   if(convg .and. it>5) exit

!simple mixing for iteration
   vi(:)=bl*vo(:)+(1.0d0-bl)*vi(:)
   vitau(:)=bl*votau(:)+(1.0d0-bl)*vitau(:)


   if(it==600 .and. .not. convg) then
     write(6,'(/a)') 'sratom_m: WARNING failed to converge'
   end if
 end do !it

 dltetot=etot-etot0

 if(.not. convg .and. ierr==0) then
   ierr=100
 end if

 deallocate(u,up,u0,u0p)
 deallocate(vo,vi1,vo1)
 deallocate (vxcsm,votausm)
 deallocate(work)
 return

 end subroutine sratom_m
