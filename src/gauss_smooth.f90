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
 subroutine gauss_smooth(rr, yy, zz, nn, ncon)

! Cheap imitation of Gaussian smoothing using convolution of input yy
! on uniformly spaced mesh xx with function (1 - x**2)**2 from -1 to 1
! intended for log mesh with approximate correction in weights for
! uneven real-space grid spacing

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr(nn),yy(nn)
 integer nn,ncon

!Output variables
 real(dp) :: zz(nn)

!Local variables

 real(dp), allocatable :: yyp(:),yypsm(:),work(:),work2(:)
 real(dp) :: wts(-100:100)
 real(dp) :: zs1(100),zs2(100),zs3(100),zs4(100),zs5(100)
 real(dp) :: al,xx,dx,sm,snorm,tt,csm,yypmin
 integer :: ii,jj,kk,nlim,next


 allocate(yyp(nn),yypsm(nn),work(nn),work2(nn))

 al = 0.01d0 * dlog(rr(101) / rr(1))

!check that ncon is odd
 if(mod(ncon,2) == 0) then
   write(6,'(/a/)') 'ERROR polysmooth must have odd ncon'
   stop
 end if

 csm=-log(1.0d-2)

!set up weights
 dx=2.0d0/dfloat(ncon+1)
 nlim=(ncon+1)/2
 snorm=0.0d0
 wts(:)=0.0d0
 do jj=-nlim,nlim
   xx=dx*dfloat(jj)
   wts(jj)=exp(-csm*xx**2 + al*jj)
!  wts(jj)=(1.0d0-xx**2)**2
   snorm=snorm+wts(jj)
 end do

 do jj=-nlim, nlim
   wts(jj)=wts(jj)/snorm
 end do
 
 zz(:)=0.0d0 
!do ii=1,nn
 do ii=nn,1,-1
   sm=0.0d0
   do jj=-nlim,nlim
     if(ii+jj <1) then
        sm=sm+wts(jj)*yy(1-jj)
     else if(ii+jj>nn) then
        sm=sm+wts(jj)*yy(nn-jj)
     else
        sm=sm+wts(jj)*yy(ii+jj)
     end if
   end do
   zz(ii)=sm
 end do

!set up first derivative for smoothing
 call derivs(nn,yy,al,rr,yyp,work,work2)

 yypsm(:)=0.0d0 
!do ii=1,nn
 do ii=nn,1,-1
   sm=0.0d0
   do jj=-nlim,nlim
     if(ii+jj <1) then
        sm=sm+wts(jj)*yyp(1-jj)
     else if(ii+jj>nn) then
        sm=sm+wts(jj)*yyp(nn-jj)
     else
        sm=sm+wts(jj)*yyp(ii+jj)
     end if
   end do
   yypsm(ii)=sm
 end do

!replace zz at smallest radii with linear extrapolation from last nominally
!reliable smoothed points using smoothed radial derivative
!for good measure, average over 5 inward-extrapolation starting points

!look for flattest point in the near region above nlim
 yypmin=1.0d10
 next=nlim+5
 do ii=nlim+5, nlim+20
   if(abs(yypsm(ii))<yypmin) then
     yypmin=abs(yypsm(ii))
     next=ii
   end if
 end do

 do ii=1,next+10
   zs1(ii)=zz(ii)
   zs2(ii)=zz(ii)
   zs3(ii)=zz(ii)
   zs4(ii)=zz(ii)
   zs5(ii)=zz(ii)
 end do

 kk=next-2
 do ii=kk-1,1,-1
   zs1(ii)=zz(kk)-yypsm(kk)*(rr(kk)-rr(ii))
 end do

 kk=next-1
 do ii=kk-1,1,-1
   zs2(ii)=zz(kk)-yypsm(kk)*(rr(kk)-rr(ii))
 end do

 kk=next
 do ii=kk-1,1,-1
   zs3(ii)=zz(kk)-yypsm(kk)*(rr(kk)-rr(ii))
 end do

 kk=next+1
 do ii=kk-1,1,-1
   zs4(ii)=zz(kk)-yypsm(kk)*(rr(kk)-rr(ii))
 end do

 kk=next+2
 do ii=kk-1,1,-1
   zs5(ii)=zz(kk)-yypsm(kk)*(rr(kk)-rr(ii))
   zz(ii)=0.2d0*(zs1(ii)+zs2(ii)+zs3(ii)+zs4(ii)+zs5(ii))
 end do

 deallocate(yyp,yypsm,work,work2)
 return
 
 end subroutine gauss_smooth
