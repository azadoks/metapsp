
! Copyright (c) 1989-2022 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
 subroutine lschfb_m(nn,ll,ierr,ee,rr,vv,vtau,uu,up,zz,mmax,mch,srel,ircmin)

!Finds bound states of an all-electron atomic potential using
!Pauli-type  scalar-relativistic Schroedinger equation

!nn  principal quantum number
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
 real(dp) :: rr(mmax),vv(mmax),vtau(mmax)
 real(dp) :: zz
 integer :: ircmin,nn,ll
 integer :: mmax
 logical :: srel

!external function
 real(dp) :: fcrossover

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee
 integer :: ierr,mch

!Local variables

 real(dp) :: aei,aeo,aii,aio,als !functions in aeo.f90
 real(dp) :: de,emax,emin
 real(dp) :: eps,fss,tfss,gamma,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap,tmp
 real(dp) :: ff, rcmin
 real(dp) :: amesh,al,xx
 integer :: ii,it,jj,imch,nint,node,nin

 real(dp), allocatable :: dv(:),fr(:),frp(:)
 real(dp), allocatable :: cfm(:),cfmp(:)
 real(dp), allocatable :: dvtau(:),upp(:)
 allocate(dv(mmax),fr(mmax),frp(mmax))
 allocate(cfm(mmax),cfmp(mmax))
 allocate(dvtau(mmax),upp(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
 eps=1.0d-14
 ierr = 100

! relativistic - non-relativistic switch
 if(srel) then
   fss=(1.0d0/137.036d0)**2
 else
   fss=1.0d-24
 end if

 if(ll==0) gamma=dsqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*dsqrt(ll**2-fss*zz**2) + &
& (ll+1)*dsqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

 sls=ll*(ll+1)

 emax=vv(mmax)+0.5d0*sls/rr(mmax)**2
 emin=0.0d0
 do ii=1,mmax
   emin=dmin1(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
 end do
 if(ee>emax) ee=1.25d0*emax
 if(ee<emin) ee=0.75d0*emin
 if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0
 cfm(:)=0.0d0
 cfmp(:)=0.0d0
 

 als=al**2

! return point for bound state convergence
!do nint=1,60
 do nint=1,100

! calculate dv/dr for darwin correction
   dv(:)=0.0d0

   dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
&         -6.d0*vv(5))/(24.d0*al*rr(1))
   dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
&         +2.d0*vv(5))/(24.d0*al*rr(2))
  
   do ii=3,mmax-2
     dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
&          -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
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
&     (rr(ii)*(1.0d0+0.5d0*tfss*(ee-vv(ii)))))
     frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(ee-vv(ii)))
   end do

    
! kill fr,frp smoothly approaching minimim psp rc
  rcmin=rr(ircmin)
! ff=0.5d0
  ff=0.1d0
  do ii=1,mmax
    xx=(rcmin-rr(ii))/(ff*rcmin) - 1.0d0
    fr(ii) = fr(ii)*fcrossover(xx)
    frp(ii)=frp(ii)*fcrossover(xx)
  end do

! metagga-corrected primary coefficient arrays for u (fr) and up (frp)
!THESE FORMS AGREE WITH MY 2016 lschfb_m (now _mt) - DO NOT DISTURG
!      cfm(ii)  = als*(sls + ((2.0d0*(vv(ii)-ee)*rr(ii)**2 &
!&             + rr(ii)*dvtau(ii))/(1.0d0+vtau(ii))))
!     cfmp(ii) = al*(1.0d0 - rr(ii)*dvtau(ii)/(1.0d0+vtau(ii)))

    do ii=1,mmax
      cfm(ii)  = als*(sls + ((2.0d0*(vv(ii)-ee)*rr(ii)**2 &
&             + rr(ii)*dvtau(ii))/(1.0d0+vtau(ii))))

      cfmp(ii) = al*(1.0d0 - rr(ii)*dvtau(ii)/(1.0d0+vtau(ii)))
    end do

! find classical turning point for matching
  mch=0

!use simpler version for wellstate_m, where large-r garbage has been removed
 if(ee>0.0d0) then
  do ii=mmax,2,-1
    if(cfm(ii-1)<=0.d0 .and. cfm(ii)>0.d0) then
      mch=ii
      exit
    end if
  end do

 else

! revised search strategy for classical turning point 
! to deal with mgga iissues at large r
   imch=0
!  do ii=1,mmax
   do ii=1,mmax-1
!    if(rr(ii)<50.0d0*rr(1)) cycle
     if(rr(ii)<0.5d0/zz) cycle
     if(cfm(ii)<0.0d0) then
        imch=ii
        exit
     end if
   end do

   if(imch==0) then
      ierr=-1
      return
   end if

   mch=0
   do ii=imch+1,mmax
     if(cfm(ii)>0.0d0) then
        mch=ii
        exit
     end if
   end do

 end if !ee < or > 0 switch

   if(mch==0 .or. mch>mmax-20) then
    ierr=-1
!   exit
    return
   end if

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
   mch=min(mch,mmax-4)
   do ii=4,mch-1
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2
       upp(ii+1)=(cfmp(ii+1) + frp(ii+1))*up(ii+1) + (cfm(ii+1) + fr(ii+1))*uu(ii+1)
       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
!    if(uu(ii+1)*uu(ii) .le. 0.0d0) node=node+1
     if(uu(ii+1)*uu(ii) .le. 0.0d0 .and. rr(ii)>0.5d0/zz) node=node+1
   end do
  
   uout=uu(mch)
   upout=up(mch)
  
  
   if(node-nn+ll+1==0) then
! start inward integration at 10*classical turning
! point with simple exponential
     nin=mmax-4
     do ii=mmax-4,1,-1
       if(rr(ii)<15.0d0*rr(mch)) then
         nin=ii
         exit
       end if
     end do
  
     if(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee) <0.0d0) then
       return
     end if
     xkap=dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))
  
     do ii=nin,nin+4
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)
       upp(ii)=(cfmp(ii) + frp(ii))*up(ii) + (cfm(ii) + fr(ii))*uu(ii)

     end do
  
! integrate inward
  
     do ii=nin,mch+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=(cfmp(ii-1) + frp(ii-1))*up(ii-1) + (cfm(ii-1) + fr(ii-1))*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do

     do ii=nin-5,mmax
       uu(ii)=0.0d0
       up(ii)=0.0d0
       upp(ii)=1.0d0
     end do
  
! scale outside wf for continuity
  
     sc=uout/uu(mch)
  
     do ii=mch,nin
       upp(ii)=sc*upp(ii)
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do
  
     upin=up(mch)
  
! perform normalization sum
  
     ro=rr(1)/dsqrt(amesh)
     sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
  
     do ii=1,nin-3
       sn=sn+al*rr(ii)*uu(ii)**2
     end do
  
     sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
&              + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
&              +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0
  
! normalize u
  
     cn=1.0d0/dsqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin
  
     do ii=1,nin
       upp(ii)=cn*upp(ii)
       up(ii)=cn*up(ii)
       uu(ii)=cn*uu(ii)
     end do
     do ii=nin+1,mmax
       uu(ii)=0.0d0
       up(ii)=0.0d0
       upp(ii)=0.0d0
     end do
  
! perturbation theory for energy shift
  
     de=0.5d0*uout*(upout-upin)/(al*rr(mch))
  
! convergence test and possible exit
  
     if(dabs(de)<dmax1(dabs(ee),0.2d0)*eps) then
       ierr = 0
       exit
     end if
  
     if(de>0.0d0) then 
       emin=ee
     else
       emax=ee
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)
  
   else if(node-nn+ll+1<0) then
! too few nodes
     emin=ee
     ee=0.5d0*(emin+emax)
  
   else
! too many nodes
     emax=ee
     ee=0.5d0*(emin+emax)
   end if
  
 
 end do

!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
   upp(:)=-upp(:)
 end if

 deallocate(dv,fr,frp)
 deallocate(cfm,cfmp)
 deallocate(dvtau,upp)
 return

 end subroutine lschfb_m
