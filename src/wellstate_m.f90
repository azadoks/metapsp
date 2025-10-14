!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
 subroutine wellstate_m(nnin,ll,irca,ep,rr,vfull,vxc,vtau,uu,up,zz, &
&                       mmax,mch,srel,ircmin)

!creates quantum well to confine positive-energy state, and calculates
!the resulting all-electron wave function

!nn  principal quantum number of well state
!ll  angular momentum
!irc  index of core radius
!ep  target energy for well state (>0)
!rr  log radial mesh
!vfull  all-electron potential
!vxc all-electron xc potential
!vtau  meta-gga dExc/d(KE density)
!uu  all-electron well-bound wave function
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log radial mesh 
!mch matching mesh point for inward-outward integrations
!srel .true. for scalar-relativistic, .false. for non-relativistic
 
 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr(mmax),vfull(mmax),vxc(mmax),vtau(mmax)
 real(dp) :: ep,zz
integer :: irca(6)
 integer :: nnin,ll,irc,mmax,ircmin
 logical :: srel

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 integer :: mch
 
!Local variables
  real(dp), allocatable :: vfullw(:),vtauw(:),vwell(:)
 real(dp) :: al,cwell,et,xx,rwell,rwmax,rwmin,rwscale,umax,sls
 real(dp) :: rwmin0,rwmax0

!external function
 real(dp) :: fcrossover

 real(dp),parameter :: eps=1.0d-8
 integer :: ii,ircmax,iter,ierr,iumax,l1,itrwell,ivmx,nn,nnloop
 logical :: convg

 allocate(vfullw(mmax),vtauw(mmax),vwell(mmax))
 
 l1=ll+1
 irc=irca(l1)
 ircmax=max(irca(1),irca(2),irca(3),irca(4))
 al = 0.01d0 * dlog(rr(101) / rr(1))

!get rid of large-radius weird behavior of atomic vxc and vtau for vwell
!construction
 do ii=1,mmax
   xx=2.0d0*(1.5d0*rr(ircmax+5)-rr(ii))/(rr(ircmax+5))
   vtauw(ii)=fcrossover(xx)*vtau(ii)
   vfullw(ii)=vfull(ii)-fcrossover(-xx)*vxc(ii)
 end do

 uu(:)=0.0d0 ; up(:)=0.0d0

 et=-0.1d0
 call lschfb_m(nnin,ll,ierr,et,rr,vfull,vtau,uu,up,zz,mmax,mch,srel,ircmin)

!if bound state is found, check its localization
 if(ierr==0) then
   umax=0.0d0
   do ii=mmax,1,-1
     if(dabs(uu(ii))>=umax) then
       umax=dabs(uu(ii))
     else
       iumax=ii+1
       exit
     end if
   end do
!if bound state is localized compared to rc, use it and its energy
   if(rr(iumax)<0.75d0*rr(irc)) then
     ep=et
     write(6,'(/a,i2,a,i2/a,f12.8/a)') &
&          'WARNING wellstate: localized bound state found for n=', &
&           nnin,' l=', ll,'WARNING this state with energy',ep, &
&          'WARNING will be used for the first projector'
     return
   end if
 end if

 uu(:)=0.0d0 ; up(:)=0.0d0

!enforce non-valence projector energy minimum
!if(ep<0.25d0) then
 if(ep<0.150) then
!  write(6,'(/a/)') 'Minimum well-state energy reset to 0.25)'
   write(6,'(/a/)') 'Minimum well-state energy reset to 0.15'
! ep=0.25
  ep=0.15
 end if

 cwell=ep+0.5d0

!loop which will change number of nodes before well wall gets too steep or
!shallow

 nn=nnin

 do nnloop=0,10

   rwmin0=0.5d0*rr(irca(l1))
   rwmax0=15.0d0
   
   rwmin=rwmin0
   rwmax=rwmax0

   convg=.false.

   rwell=0.5d0*(rwmin0+rwmax0)

   do itrwell=1,100

!if well is getting too small, increment nn and drop out of itrwell loop
   if(rwell<1.05d0*rwmin0) then
      nn=nn+1
      exit
   end if

!if well is getting too large,decrement cwell and drop out of itrwell loop
   if(rwell>0.95d0*rwmax0) then
      cwell=ep+0.9d0*(cwell-ep)
      exit
   end if

!create well potential
     do ii=1,mmax
       vwell(ii)=vfullw(ii)
     end do

     rewind 30
!   start outside rc to keep numerical derivatives at rc accurate
     do ii=ircmax+6,mmax
       xx=(rr(ii)-rr(irc+5))/rwell
      vwell(ii)=vwell(ii)+cwell*xx**4/(1.0d0+xx**4)
     end do

       
!find bound state in well
     flush 6
     et=ep
     uu(:)=0.0d0 ; up(:)=0.0d0
     call lschfb_m(nn,ll,ierr,et,rr,vwell,vtauw,uu,up,zz,mmax,mch,srel,ircmin)

     if(abs(et-ep)<eps) then
!     ep=et
      convg=.true.
      exit
     end if

!Interval-halving search after proper rwell has been bracketed
!Increment or decrement nn if search gets too close to initial
!limits
     if(rwmin>0.0d0 .and. rwmax>0.0d0) then
       if(et>ep) then
         rwmin=rwell
       else
         rwmax=rwell
      end if
      rwell=0.5d0*(rwmax+rwmin) 
      cycle
     end if

     if(convg) exit

   end do !itrwell

   if(convg) then
!    nnin=nn
     exit
   end if
 end do !nnloop

 if(.not. convg) then
   write(6,'(a,a,i3,a,i3,a,f8.4)') &
&   'ERROR wellstate: well potential iteration failed', &
&   ' to converge, n=',nn,' l=',ll,' ep=',ep
  stop
 end if

 if(cwell>0.0d0) then
  write(6,'(/a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
&          '   Wellstate for l =',ll,'  n =',nn, &
&          '     eigenvalue = ',et, &
&          '     asymptotic potential = ',cwell, &
&          '     half-point radius =',rr(irc)+rwell
 end if

  deallocate(vfullw,vtauw,vwell)
 return
 end subroutine wellstate_m
