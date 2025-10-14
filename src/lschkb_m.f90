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
 subroutine lschkb_m(nn,ll,ierr,ee,vkb,rr,vv,vtau,dvtau,uu,up,mmax,mch)

! outward integration of the inhomogeneous radial Schroedinger equation
! on a logarithmic mesh with local potential and one proector term

!nn  principal quantum number (not used)
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!vkb Vanderbilt-Kleinman-bylander projector
!rr  log radial mesh
!vv  local pseudopotential
!vtau  meta-gga dExc/d(KE density) pseudo ground state
!uu  wave function
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: mmax,mch
 real(dp) :: rr(mmax),vv(mmax),vkb(mmax),vtau(mmax),dvtau(mmax)
 real(dp) :: zz
 integer :: nn,ll
 
!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee
 integer :: ierr


!Local variables
 real(dp) :: amesh,al
 real(dp) :: aeo, aio, als, cn, tmp
 real(dp) :: sls, uout, upout
 real(dp) :: akb,ckb
 integer :: ii, it

 real(dp), allocatable :: upp(:),cf(:)
 
 allocate(upp(mmax),cf(mmax))


 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

 ierr = 0

 sls=ll*(ll+1)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0

 als=al**2

! coefficient array for uu in differential eq.
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
 end do

! start wavefunction with series based on projector

 ckb = 2.0d0*vkb(1)/rr(1)**(ll+1)
 akb = ckb/(6.0d0+4.0d0*ll)
 do ii=1,4
   uu(ii)=akb*rr(ii)**(ll+3)
   up(ii)= al*(ll+3)*uu(ii)
!  upp(ii)= als*(ll+3)**2*uu(ii)

   tmp=rr(ii)*dvtau(ii)*(al*up(ii) - als*uu(ii)) &
&       +vtau(ii)*(        - al*up(ii) - als*sls*uu(ii))
   upp(ii)=(al*up(ii) + cf(ii)*uu(ii) - tmp &
&        + 2.0d0*als*vkb(ii)*rr(ii)**2)/(1.0d0+vtau(ii))
 end do

! outward integration using predictor once, corrector
! twice

 do ii=4,mch-1
   uu(ii+1)=uu(ii)+aeo(up,ii)
   up(ii+1)=up(ii)+aeo(upp,ii)
   do it=1,2
       tmp=rr(ii+1)*dvtau(ii+1)*(al*up(ii+1) - als*uu(ii+1)) &
&           +vtau(ii+1)*(        - al*up(ii+1) - als*sls*uu(ii+1))
       upp(ii+1)=(al*up(ii+1) + cf(ii+1)*uu(ii+1) - tmp &
&            + 2.0d0*als*vkb(ii+1)*rr(ii+1)**2)/(1.0d0+vtau(ii+1))

!     upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1) &
!&            + 2.0d0*als*vkb(ii+1)*rr(ii+1)**2

     up(ii+1)=up(ii)+aio(upp,ii)
     uu(ii+1)=uu(ii)+aio(up,ii)
   end do
 end do

 uout=uu(mch)
 upout=up(mch)

 deallocate(upp,cf)

 return
 end  subroutine lschkb_m
