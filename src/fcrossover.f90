!
! Copyright (c) 1989-2015 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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

! fcrossover varies from 0 to 1 on [-1,1] with 3 zero derivatives at the ends
! value at xx=0 is 0.5, slope is 1.23048675 = 945/(2*384)
!
 function fcrossover(xx)

 implicit none
 integer, parameter :: dp=kind(1.0d0)

 real(dp) :: fcrossover,xx,yy

 if(xx<-1.0d0) then
   fcrossover = 0.0d0
 else if(xx>1.0d0) then
   fcrossover = 1.0d0
 else

   yy=xx**9/9.0d0 - 4.0d0*xx**7/7.0d0 + 6.0d0*xx**5/5.0d0 &
&     - 4.0d0*xx**3/3.0d0 + xx

   yy=(945.0d0/384.0d0)*yy

   fcrossover=0.5d0*(1.0d0 + yy)

 end if

 return
 end function fcrossover
