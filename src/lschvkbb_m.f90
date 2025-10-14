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
 subroutine lschvkbb_m(nn,ll,nvkb,ierr,ee,emin,emax,ebar, &
&   rr,vloc,vtau,vkb,evkb,uu,up,mmax,mch)

! Finds bound states of a  pseudopotential with
! Vanderbilt-Kleinman-Bylander non-local projectors


!nn  principal quantum number
!ll  angular-momentum quantum number
!nvkb  = number of VKB projectors to be used
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!emin  externally generaated estimate of lower bound for ee
!emax  externally generaated estimate of upper bound for ee
!ebar  <kinetic + local + non-local> /= ee if vtau/=0
!rr  log radial mesh
!vloc  local part of psp
!vtau  meta-gga dExc/d(KE density) pseudo ground state
!vkb  VKB projectors
!evkb coefficients of BKB projectors
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input Variables
 real(dp) :: emin,emax
 real(dp) :: rr(mmax),vloc(mmax),vtau(mmax),vkb(mmax,nvkb),evkb(nvkb)
 integer :: nn,ll,nvkb,mmax

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee  !in/out, really - needs starting guess
 real(dp) :: ebar
 integer :: ierr,mch


!Local variables
 real(dp) :: aei,aii,cn
 real(dp) :: de
 real(dp) :: eps,ro,sc
 real(dp) :: sls,sn,uout,upin,upout
 real(dp) :: amesh,al,als,tmp
 real(dp) :: rc,xkap
 integer :: ii,it,jj,nin,nint,node

 real(dp), allocatable :: upp(:),cf(:),dvtau(:)
!New for vkb total energy contribution
 real(dp) :: gg0(4)
 real(dp) :: etest,eloc,enloc,ekin
 real(dp), allocatable :: vkbt(:),work(:),tau(:)

 allocate(upp(mmax),cf(mmax),dvtau(mmax))
 allocate(vkbt(mmax),work(mmax),tau(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)
 node=0
! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
!eps=1.0d-10
 eps=1.0d-8
 ierr = 100

 sls=ll*(ll+1)

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

!if(ee>emax) ee=1.25d0*emax
!if(ee<emin) ee=0.75d0*emin
!if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0

 als=al**2


! return point for bound state convergence
 do nint=1,100

! coefficient array for u in differential eq.
   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vloc(ii)-ee)*rr(ii)**2
   end do
  
! find classical turning point for matching
   mch=1
   do ii=mmax,2,-1
     if(cf(ii-1)<=0.d0 .and. cf(ii)>0.d0) then
       mch=ii
       exit
     end if
   end do

   if(mch==1 .and. nvkb==0) then
    write(6,'(/a)') 'lschvkbb: ERROR no classical turning point'
    write(6,'(a,i2,a,i6)') 'lschvkbb: ERROR nvkb=',nvkb,'mch=',mch
    write(6,'(a,i2,a,f8.4,a)') 'lschvkbb: ERROR l=',ll,'  e=',ee
!   stop
   end if

   if(nvkb>0) then
! find cutoff radius for projectors
     rc=0.0d0
     do ii=mmax,1,-1
       if(abs(vkb(ii,1))>0.0d0) then
         rc=rr(ii)
         exit
       end if
     end do

     if(mch==1) rc=1.25d0*rc

! adjust matching radius if necessary
     if(rr(mch)<rc) then
      do ii=mch,mmax
       if(rr(ii)>rc) then
        mch=ii
        exit
       end if
      end do
     end if
   end if
! outward integration
   call vkboutwf_m(ll,nvkb,ee,vkb,evkb,rr,vloc,vtau,dvtau, &
&                  uu,up,node,mmax,mch)
  
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
  
     xkap=dsqrt(sls/rr(nin)**2 + 2.0d0*(vloc(nin)-ee))
  
     do ii=nin,nin+4
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)

       tmp=rr(ii)*dvtau(ii)*(al*up(ii) - als*uu(ii)) &
&           +vtau(ii)*(        - al*up(ii) - als*sls*uu(ii))
       upp(ii)=(al*up(ii) + cf(ii)*uu(ii) - tmp)/(1.0d0+vtau(ii))
     end do
  
! integrate inward
  
     do ii=nin,mch+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2

         tmp=rr(ii-1)*dvtau(ii-1)*(al*up(ii-1) - als*uu(ii-1)) &
&             +vtau(ii-1)*(        - al*up(ii-1) - als*sls*uu(ii-1))
         upp(ii-1)=(al*up(ii-1) + cf(ii-1)*uu(ii-1) - tmp)/(1.0d0+vtau(ii-1))

         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do
  
! scale outside wf for continuity
  
     sc=uout/uu(mch)
  
     do ii=mch,nin+4
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do
  
     upin=up(mch)
  
! perform normalization sum
  
     ro=rr(1)/dsqrt(amesh)
     sn=ro**(2*ll+3)/dfloat(2*ll+3)
  
!It's very silly to have to put this exit statement, but I got floating-
!point error aborts when this loop tried, for example, to square 4.572036-154
     do ii=1,nin-3
      if(abs(uu(ii))<1.0d-145) then
         nin=ii-1
         exit
      end if
       sn=sn+al*rr(ii)*uu(ii)**2
     end do
  
     if(abs(uu(nin))>1.0d-145) then
       sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
&                + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
&                +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0
      end if
  
! normalize u
  
     cn=1.0d0/dsqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin

     
  
     do ii=1,nin
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
!      emin=ee
       emin=0.5d0*(emin+ee)
     else
!      emax=ee
       emax=0.5d0*(emax+ee)
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)
  
   else if(node-nn+ll+1<0) then
! too few nodes
!    emin=ee
     emin=0.5d0*(emin+ee)
     ee=0.5d0*(emin+emax)
  
   else
! too many nodes
!    emax=ee
     emax=0.5d0*(emax+ee)
     ee=0.5d0*(emin+emax)
   end if
  
 end do

 if(ierr>0) return

!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

!average kinetic, potential, and non-local energies, not exactly the left
!side of the Schroedinger equation unless vtau=0.  It can be set to zero
!when this is called in which case ebar should be identical to the
!calculated eigenvalue ee.

 do ii=1,nin
   tau(ii)= 0.5d0*(((up(ii)/al - uu(ii))/rr(ii)**2)**2 &
&          + (sls/rr(ii)**2)*(uu(ii)/rr(ii))**2)
 end do

 call vpinteg(rr**2,tau,nin,2*ll+2,ekin,rr,mmax)

 work(:)=vloc(:)*uu(:)
 call vpinteg(uu,work,nin,2*ll+2,eloc,rr,mmax)

 vkbt(:)=0.0d0
 do jj=1,nvkb
   call vpinteg(uu,vkb(1,jj),mch,2*ll+2,gg0(jj),rr,mmax)
   vkbt(:)=vkbt(:)+gg0(jj)*evkb(jj)*vkb(:,jj)
 end do

 call vpinteg(uu,vkbt,nin,2*ll+2,enloc,rr,mmax)

 ebar=ekin+eloc+enloc


 deallocate(upp,cf,dvtau)
 deallocate(vkbt,work,tau)

 return
 end subroutine lschvkbb_m
