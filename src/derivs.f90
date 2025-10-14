!
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
 
subroutine derivs(mmax, rho, al, rr, dpr, dppr, dlap)
! NB: presumes incoming rho is properly scaled without 4pi or r**2 factors.
 
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 
 integer, intent(in) :: mmax
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: rho(mmax)
 real(dp), intent(out) :: dpr(mmax), dppr(mmax), dlap(mmax)
 
! local vars
 integer :: i
 real(dp) :: dpn(mmax), dppn(mmax)
 real(dp) :: c11,c12,c13,c14,c15
 real(dp) :: c21,c22,c23,c24,c25
 
 c11 =   2.0d0 / 24.0d0
 c12 = -16.0d0 / 24.0d0
 c13 =   0.0d0 / 24.0d0
 c14 =  16.0d0 / 24.0d0
 c15 =  -2.0d0 / 24.0d0
!
 c21 =   -1.0d0 / 12.0d0
 c22 =   16.0d0 / 12.0d0
 c23 =  -30.0d0 / 12.0d0
 c24 =   16.0d0 / 12.0d0
 c25 =   -1.0d0 / 12.0d0
!
! n derivatives of d
!     
 i=1
 dpn(i) = -25.d0/12.d0*rho(i) +4.d0*rho(i+1) -3.d0*rho(i+2) &
&         +4.d0/3.d0*rho(i+3) -1.d0/4.d0*rho(i+4)
 dppn(i) = 15.d0/4.d0*rho(i) -77.d0/6.d0*rho(i+1) +107.d0/6.d0*rho(i+2) &
&         -13.d0*rho(i+3) +61.d0/12.d0*rho(i+4) -5.d0/6.d0*rho(i+5)
 i=2
 dpn(i) = -25.d0/12.d0*rho(i) +4.d0*rho(i+1) -3.d0*rho(i+2)  &
&         +4.d0/3.d0*rho(i+3) -1.d0/4.d0*rho(i+4)
 dppn(i) = 15.d0/4.d0*rho(i) -77.d0/6.d0*rho(i+1) +107.d0/6.d0*rho(i+2) &
&         -13.d0*rho(i+3) +61.d0/12.d0*rho(i+4) -5.d0/6.d0*rho(i+5)
 
 do i = 3, mmax - 2
   dpn(i) =  c11*rho(i-2) + c12*rho(i-1) + c14*rho(i+1) + c15*rho(i+2)
   dppn(i) = c21*rho(i-2) + c22*rho(i-1) + c23*rho(i)   + c24*rho(i+1) &
&           +c25*rho(i+2)
 end do
 
 i=mmax-1
 dpn(i) = +25.d0/12.d0*rho(i) -4.d0*rho(i-1) +3.d0*rho(i-2) &
&         -4.d0/3.d0*rho(i-3) +1.d0/4.d0*rho(i-4)
 dppn(i) = -15.d0/4.d0*rho(i) +77.d0/6.d0*rho(i-1) -107.d0/6.d0*rho(i-2) &
&          +13.d0*rho(i-3) -61.d0/12.d0*rho(i-4) +5.d0/6.d0*rho(i-5)
 i=mmax
 dpn(i) = +25.d0/12.d0*rho(i) -4.d0*rho(i-1) +3.d0*rho(i-2) &
&         -4.d0/3.d0*rho(i-3) +1.d0/4.d0*rho(i-4)
 dppn(i) = -15.d0/4.d0*rho(i) +77.d0/6.d0*rho(i-1) -107.d0/6.d0*rho(i-2) &
&          +13.d0*rho(i-3) -61.d0/12.d0*rho(i-4) +5.d0/6.d0*rho(i-5)
 
!
! r derivatives of d
!
 do i = 1, mmax
   dpr(i) = dpn(i) / (al * rr(i))
   dppr(i) = (dppn(i) - al * dpn(i)) / (al * rr(i))**2
   dlap(i) = (dppn(i) + al * dpn(i)) / (al * rr(i))**2
 end do
 
end subroutine derivs
