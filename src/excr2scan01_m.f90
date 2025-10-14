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
! calculates Hartree and exchange-correlation potentials and total-energy
! term to add to eigenvalue sum

! adapted from sections of Natalie's excor.F90

! My rho(:) is r(:)**2 times hers  (Rydbergs -> Hartrees)

! My tau(:) is 0.5*r(:)**2 times hers



  subroutine  excr2scan01_m(al,rho,tau,vxc,vtau,exc,rr,mmax)

 use r2scanmod
 implicit none
 
 integer, parameter :: dp=kind(1.0d0)

!real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 
 real(dp), intent(in) :: al
  real(dp), intent(in) :: rho(mmax),tau(mmax),rr(mmax)
  real(dp), intent(out) :: vxc(mmax),vtau(mmax),exc(mmax)
 integer, intent(in) :: mmax
 
! local
!real(dp) :: eta
 real(dp), allocatable :: v(:), vt(:), e(:)
 real(dp), allocatable :: dpr(:), dppr(:), dlap(:)
 real(dp), allocatable :: rho_loc(:),tau_loc(:)
 real(dp), allocatable :: dexcdn(:),dexcds(:)
 real(dp), allocatable :: tmp1(:)
 integer :: ii

 allocate(v(mmax), vt(mmax), e(mmax))
 allocate(dpr(mmax), dppr(mmax), dlap(mmax))
 allocate(rho_loc(mmax),tau_loc(mmax))
 allocate(dexcdn(mmax),dexcds(mmax))
 allocate(tmp1(mmax))

  v(:)=0.0d0 ;  vt(:)=0.0d0 ;  e(:)=0.0d0
  dpr(:)=0.0d0 ;  dppr(:)=0.0d0 ;  dlap(:)=0.0d0
  rho_loc(:)=0.0d0 ; tau_loc(:)=0.0d0
  dexcdn(:)=0.0d0 ; dexcds(:)=0.0d0
  tmp1(:)=0.0d0
!
!conf = (3.d0*pi**2)**thrd
!conrs = (3.d0/(4.d0*pi))**thrd

!both Natalie's and my rho and tau seem to have the 4pi factors in them
 
 rho_loc(:) = rho(:) /(4.0d0*pi)

 tau_loc(:) = tau(:) /(4.0d0*pi)

 call derivs(mmax, rho_loc, al, rr, dpr, dppr, dlap)

!r2scan smoothing parameter
!moved to main routine oncvpsp_m.f90
!eta = 0.01d0
!eta = 0.002d0
!eta = 0.001d0

!call r2scaninit(eta)

!call r2scanfun(rho,grad,tau,exc,vtau,vxcn,vxcs)

 do ii=1,mmax
   call r2scanfun(rho_loc(ii),dpr(ii),tau_loc(ii), &
&                 exc(ii),vtau(ii),dexcdn(ii),dexcds(ii))

 end do


!# quotes from Natalie's excor.F90
!# dum(1:n)=dexcds(1:n)
!# call derivative(Grid,dum,dum1,1,n)
!# rvxc(1:n)=(dexcdn(1:n)-dum1(1:n))*Grid%r(1:n)-2*dexcds(1:n)
 
!my interpretations of the above

 call derivs(mmax, dexcds, al, rr, tmp1, dppr, dlap)

!DRH change (Natalie's vxc is my  r*vxc)
!vxc(:) = dexcdn(:)-rr(:)*tmp1(:) - 2.0d0*dexcds(:)
!vxc(:) = dexcdn(:)-tmp1(:) - 2.0d0*dexcds(:)/rr(:)
 vxc(:) = dexcdn(:)-2.0d0*tmp1(:) - 4.0d0*dexcds(:)/rr(:)



! go to convention for Don Hamman code - I do not understand:
! in principle libxc outputs energy per particle,
! and here we divide by rho again to get an LDA exc
! identical to the other routine...
 do ii=1,mmax
   exc(ii) = exc(ii) / max(rho_loc(ii),1.d-20)
 end do


 deallocate(v, vt, e)
 deallocate(dpr, dppr, dlap)
 deallocate(rho_loc,tau_loc)
 deallocate(dexcdn,dexcds)
 deallocate(tmp1)

 end subroutine  excr2scan01_m
