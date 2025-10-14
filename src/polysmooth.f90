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
 subroutine polysmooth(yy, zz, nn, ncon)

! Cheap imitation of Gaussian smoothing using convolution of input yy
! on uniformly spaced mesh xx with function (1 - x**2)**2 from -1 to 1
! intended for log mesh ignoring uneven real-space grid spacing

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: yy(nn)
 integer nn,ncon

!Output variables
 real(dp) :: zz(nn)

!Local variables

 real(dp) :: wts(-50:50)
 real(dp) :: xx,dx,sm,snorm,tt
 integer :: ii,jj,nlim

!check that ncon is odd
 if(mod(ncon,2) == 0) then
   write(6,'(/a/)') 'ERROR polysmooth must have odd ncon'
   stop
 end if

!set up weights
 dx=2.0d0/dfloat(ncon+1)

 nlim=(ncon-1)/2

 snorm=0.0d0
 wts(:)=0.0d0
 do jj=-nlim,nlim
  xx=dx*dfloat(jj)
  wts(jj)=(1.0d0-xx**2)**2
  snorm=snorm+wts(jj)
 end do

 do jj=-nlim, nlim
   wts(jj)=wts(jj)/snorm
 end do
 
 zz(:)=0.0d0
 do ii=1,nn
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
 
 return
 
 end subroutine polysmooth
