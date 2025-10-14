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
 subroutine lschfs_m(nn,ll,ierr,ee,rr,vv,vtau,uu,up,zz,mmax,mch,srel,ircmin)

! integrates radial Pauli-type scalar-relativistic equation
! on a logarithmic mesh
! modified routine to be used in finding norm-conserving
! pseudopotential

!nn  effective principal quantum number based on nodes inside mch (output)
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!rr  log radial mesh
!vv  local atomic potential
!vtau  meta-gga dExc/d(KE density)
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations
!srel .true. for scalar-relativistic, .false. for non-relativistic
!ircmin  minimum core radii index

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: mmax
 real(dp) :: rr(mmax),vv(mmax),vtau(mmax)
 real(dp) :: zz
 integer :: ll,mch,ircmin
 logical :: srel

!external function
 real(dp) :: fcrossover
 
!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee
 integer :: ierr,nn


!Local variables
 real(dp) :: amesh,al
 real(dp) :: aei,aeo,aii,aio,als,cn !functions in aeo.f90
!real(dp) :: aeo, aio, als, cn
 real(dp) :: fss, tfss, gamma, ro, xx, ff, rcmin
 real(dp) :: sls, sn, uout, upout, tmp
 integer :: ii, it,node

 real(dp), allocatable :: upp(:),cf(:),dv(:),fr(:),frp(:),dvtau(:)
 real(dp), allocatable :: cfm(:),cfmp(:)
 
 allocate(upp(mmax),cf(mmax),dv(mmax),fr(mmax),frp(mmax),dvtau(mmax))
 allocate(cfm(mmax),cfmp(mmax))


 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 ierr = 0

! relativistic - non-relativistic switch
 if(srel) then
   fss=(1.0d0/137.036d0)**2
 else
   fss=1.0d-24
 end if


 if(ll==0) gamma=sqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*sqrt(ll**2-fss*zz**2) &
& +(ll+1)*sqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

 sls=ll*(ll+1)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0
 cfm(:)=0.0d0
 cfmp(:)=0.0d0
 cf(:)=0.0d0

 node=0

 als=al**2

! coefficient array for u in differential eq.
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
 end do

! calculate dv/dr for darwin correction
 dv(:)=0.0d0
 
 dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
&       -6.d0*vv(5))/(24.d0*al*rr(1))
 dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
&       +2.d0*vv(5))/(24.d0*al*rr(2))

 do ii=3,mmax-2
   dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
&         -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
 end do
! calculate dvtau/dr for meta-gga operator
   dvtau(:)=0.0d0

   dvtau(1)=(-50.d0*vtau(1)+96.d0*vtau(2)-72.d0*vtau(3)+32.d0*vtau(4) &
&         -6.d0*vtau(5))/(24.d0*al*rr(1))
   dvtau(2)=(-6.d0*vtau(1)-20.d0*vtau(2)+36.d0*vtau(3)-12.d0*vtau(4) &
&         +2.d0*vtau(5))/(24.d0*al*rr(2))

   do ii=3,mmax-2
     dvtau(ii)=(2.d0*vtau(ii-2)-16.d0*vtau(ii-1)+16.d0*vtau(ii+1) &
&          -2.d0*vtau(ii+2))/(24.d0*al*rr(ii))
   end do

!  relativistic coefficient arrays for u (fr) and up (frp).
 do ii=1,mmax
   tfss=fss
   fr(ii)=als*(rr(ii)**2)*(-tfss*(vv(ii)-ee)**2 + 0.5d0*tfss*dv(ii)/ &
&   (rr(ii)*(1.0d0+0.5d0*tfss*(ee-vv(ii)))))
   frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(ee-vv(ii)))
 end do

!  metagga-corrected primary coefficient arrays for u (fr) and up (frp)

!THESE FORMS AGREE WITH MY 2016 lschfb_m (now _mt) - DO NOT DISTURG
!      cfm(ii)  = als*(sls + ((2.0d0*(vv(ii)-ee)*rr(ii)**2 &
!&             + rr(ii)*dvtau(ii))/(1.0d0+vtau(ii))))
!     cfmp(ii) = al*(1.0d0 - rr(ii)*dvtau(ii)/(1.0d0+vtau(ii)))

    do ii=1,mmax
      cfm(ii)  = als*(sls + ((2.0d0*(vv(ii)-ee)*rr(ii)**2 &
&             + rr(ii)*dvtau(ii))/(1.0d0+vtau(ii))))

      cfmp(ii) = al*(1.0d0 - rr(ii)*dvtau(ii)/(1.0d0+vtau(ii)))
    end do

! kill fr,frp smoothly approaching minimim psp rc

   rcmin=rr(ircmin)
!  ff=0.5d0
   ff=0.1d0
   do ii=1,mmax
     xx=(rcmin-rr(ii))/(ff*rcmin) - 1.0d0
     fr(ii) = fr(ii)*fcrossover(xx)
     frp(ii)=frp(ii)*fcrossover(xx)
   end do

! start wavefunction with series
   do ii=1,4
    if(ll==0) then
     uu(ii)=rr(ii)**gamma - zz*rr(ii)**2
     up(ii)=al*gamma*rr(ii)**gamma - 2.0d0*al*zz*rr(ii)**2
    else
     uu(ii)=rr(ii)**gamma
     up(ii)=al*gamma*rr(ii)**gamma
    end if
    upp(ii)=(cfmp(ii) + frp(ii))*up(ii) + (cfm(ii) + fr(ii))*uu(ii)
   end do


!     tmp=rr(ii)*dvtau(ii)*(al*up(ii) - als*uu(ii)) &
!&         +vtau(ii)*(        - al*up(ii) - als*sls*uu(ii))
!     upp(ii)=frp(ii)*up(ii)+fr(ii)*uu(ii) &
!&           +(al*up(ii) + cf(ii)*uu(ii) - tmp)/(1.0d0+vtau(ii))
!   end do

! New approach to obtain better continuity at the ii=4 and ii=5 junction
! Integrate out to ii=10, then back into ii=1 and use the new ii=1-4
! values to re-start the outward integration

   do ii=4,9
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2
       upp(ii+1)=(cfmp(ii+1) + frp(ii+1))*up(ii+1) + (cfm(ii+1) + fr(ii+1))*uu(ii+1)
       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
!    if(uu(ii+1)*uu(ii) .le. 0.0d0) node=node+1
   end do

! integrate inward
  
     do ii=5,2,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=(cfmp(ii-1) + frp(ii-1))*up(ii-1) + (cfm(ii-1) + fr(ii-1))*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do

  
! outward integration using predictor once, corrector
! twice
   node=0
  
   do ii=4,mch-1
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2

!       tmp=rr(ii+1)*dvtau(ii+1)*(al*up(ii+1) - als*uu(ii+1)) &
!&           +vtau(ii+1)*(        - al*up(ii+1) - als*sls*uu(ii+1))
!       upp(ii+1)=frp(ii+1)*up(ii+1)+fr(ii+1)*uu(ii+1) &
!&             +(al*up(ii+1) + cf(ii+1)*uu(ii+1) - tmp)/(1.0d0+vtau(ii+1))

       upp(ii+1)=(cfmp(ii+1) + frp(ii+1))*up(ii+1) + (cfm(ii+1) + fr(ii+1))*uu(ii+1)


       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
!    if(uu(ii+1)*uu(ii) .le. 0.0d0) node=node+1
!    if(uu(ii+1)*uu(ii) .le. 0.0d0 .and. rr(ii)>0.5d0/zz) node=node+1
     if(uu(ii+1)*uu(ii) .le. 0.0d0 .and. rr(ii)>0.05d0/zz) node=node+1
   end do
  
 uout=uu(mch)
 upout=up(mch)

!perform normalization sum

 ro=rr(1)/dsqrt(amesh)
 sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)

 do ii=1,mch-3
   sn=sn+al*rr(ii)*uu(ii)**2
 end do

 sn =sn + al*(23.0d0*rr(mch-2)*uu(mch-2)**2 &
&           + 28.0d0*rr(mch-1)*uu(mch-1)**2 &
&          +  9.0d0*rr(mch  )*uu(mch  )**2)/24.0d0

!normalize u

 cn=1.0d0/dsqrt(sn)
 uout=cn*uout
 upout=cn*upout

 do ii=1,mch
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
 end do
 do ii=mch+1,mmax
   uu(ii)=0.0d0
 end do

!calculate effective principal quantum number as if this were a bound
!state with a barrier at mch
 nn=node+ll+1

 deallocate(upp,cf,dv,fr,frp,dvtau)
 deallocate(cfm,cfmp)

 return
 end  subroutine lschfs_m
