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
 program oncvpsp_m
!
! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
! and references therein.
!
!   D. R. Hamann
!   Mat-Sim Research LLC
!   P.O. Box 742
!   Murray Hill, NJ 07974
!   USA
!
!   Developed from original "gncpp" code of March 8,1987
!
!   Output format for ABINIT pspcod=8 and upf format for quantumespresso
!
 Use m_psmlout, only: psmlout, copy_input_file_for_psml
 use m_psmlout, only: PSML_PATCH
 use m_psmlout, only: get_ec_hints

 Use globalmath_loc
 Use r2scanmod


 implicit none
 integer, parameter :: dp=kind(1.0d0)

!
 integer :: ii,ierr,iexc,iexct,ios,iprint,irps,it,icmod,lpopt
 integer :: jj,kk,ll,l1,lloc,lmax,lt,inline
 integer :: mch,mchf,mmax,n1,n2,nc,nlim,nlloc,nlmax,irpsh,nrl
 integer :: nv,irct,ncnf,nvt
 integer :: iprj,mxprj
 integer :: ircmax,ircmin,ipsmax,irx,irxt,irxmin,irxmax,irxn,iter,iter2
 integer,allocatable :: npa(:,:)
!
 integer :: dtime(8),na(30),la(30),np(6)
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 integer :: nat(30),lat(30),indxr(30),indxe(30)
 integer :: irc(6),nodes(4)
 integer :: nproj(6),npx(6),lpx(6)
 integer :: ncon(6),nbas(6)

 real(dp) :: fcrossover

 real(dp) :: al,amesh,csc,csc1,deblt,depsh,depsht,drl,eeel
 real(dp) :: ekvnl,eexc
 real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t
 real(dp) :: et,etest,ebar,emin,sls
 real(dp) :: fcfact,rcfact,dvloc0
 real(dp) :: rr1,rcion,rciont,rcmax,rlmax,rpkt
 real(dp) :: sf,zz,zion,zval,etot
 real(dp) :: xx,fcr,drtau,rdmin,tdmin,mrho0,mtau0
 real(dp) :: epstot
 real(dp) :: xdummy,fpt
 real(dp) :: rcr,xcr
 real(dp) :: dltetot,taupsmax,taudiff0,rhodiff0,rhopsmax
!
 real(dp) :: cl(6),debl(6),ea(30),ep(6),fa(30),facnf(30,5)
 real(dp) :: eapbe(30),eacopy(30),srerr(30)
 real(dp) :: fat(30,2)
 real(dp) :: fnp(6),fp(6)
 real(dp) :: qcut(6),qmsbf(6),rc(6),rc0(6)
 real(dp) :: rpk(30)
 real(dp) :: epx(6),fpx(6)

 real(dp), allocatable :: evkb(:,:),cvgplt(:,:,:,:),qq(:,:)
 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: rho(:),rhoc(:),rhot(:)
 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: vp(:,:),vfull(:),vkb(:,:,:),pswf(:,:,:)
 real(dp), allocatable :: pspwf(:,:,:),psppwf(:,:,:)
 real(dp), allocatable :: vwell(:)
 real(dp), allocatable :: vpuns(:,:)
 real(dp), allocatable :: vo(:),vxc(:)
 real(dp), allocatable :: vh(:),vtau(:),voep(:),tau(:),tauc(:)
 real(dp), allocatable :: taumod(:),taumodps(:),vtaumodps(:)
 real(dp), allocatable :: vhps(:),vxcps(:),vtaups(:),rhops(:),taups(:)
 real(dp), allocatable :: rhomod(:,:),rhoae(:,:),tautae(:)
 real(dp), allocatable :: rhomodps(:),rhodiff(:),taudiff(:)
 real(dp), allocatable :: rhov(:),tauv(:)
 real(dp), allocatable :: rholoc(:),tauloc(:)
 real(dp), allocatable :: uupsa(:,:) !pseudo-atomic orbitals array
 real(dp), allocatable :: epa(:,:),fpa(:,:)
 real(dp), allocatable :: uua(:,:),upa(:,:)
 real(dp), allocatable :: u0(:),u0p(:)
 real(dp), allocatable :: work(:)

 character*2 :: atsym,tstsr
 character*4 :: psfile

 logical :: srel,cset
 logical :: nderr

 write(6,'(a/a/a/)') &
&      'METAPSP  (Metagga psuedopotential code based on ONCVPSp ', &
&      'with generalized norm-consiervatio and convergence optimization )', &
&      'alpha version 1.0.1 05/16/2024'

 write(6,'(a/a/a/)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&      'in any publication utilizing these pseudopotentials.'

!srel=.true.
 srel=.false.

!test if scalar relativistic calculation is wanted
 inquire(file='SR', exist=srel)

 nproj(:)=0
 fcfact=0.0d0
 rcfact=0.0d0
 rc(:)=0.0d0
 ep(:)=0.0d0


 if(srel) then
   write(6,'(/a)') 'Scalar-relativistic Calculation'
 else
     write(6,'(a/)') 'Non-relativistic Calculation'
 end if


! read input data
 inline=0

! atom and reference configuration
 call cmtskp(inline)
 read(5,*,iostat=ios) atsym,zz,nc,nv,iexc,psfile
 call read_error(ios,inline)
!DEFAULT SUBSTITUTONS FOR META-INCOMPATABLE DATA
 iexc=5

!TEMPORARY
!psfile='none'
!psfile='psp8'
!psfile='upf'
!psfile='psml'

 call cmtskp(inline)
 do ii=1,nc+nv
   read(5,*,iostat=ios) na(ii),la(ii),fa(ii)
   call read_error(ios,inline)
 end do

! pseudopotential and optimization
 call cmtskp(inline)
 read(5,*,iostat=ios) lmax
 call read_error(ios,inline)

 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,rc(l1),ep(l1),ncon(l1),nbas(l1),qcut(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
 end do

! local potential
 call cmtskp(inline)
 read(5,*, iostat=ios) lloc,lpopt,rc(5),dvloc0
 call read_error(ios,inline)

! Vanderbilt-Kleinman-Bylander projectors
 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,nproj(l1),debl(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
 end do

! model core charge
 icmod=0 ; fcfact=0.0d0 ; rcfact=0.0d0
 call cmtskp(inline)
 read(5,*, iostat=ios) icmod, fcfact
 if(ios==0 .and. icmod==5) then
  backspace(5)
  read(5,*, iostat=ios) icmod, fcfact, rcfact
 end if
 call read_error(ios,inline)

! log derivative analysis
 call cmtskp(inline)
 read(5,*,iostat=ios) epsh1, epsh2, depsh
 call read_error(ios,inline)

! output grid
 call cmtskp(inline)
 read(5,*,iostat=ios) rlmax,drl
 call read_error(ios,inline)

! test configurations
 call cmtskp(inline)
 read(5,*,iostat=ios) ncnf
 call read_error(ios,inline)

 do jj=2,ncnf+1
   call cmtskp(inline)
   read(5,*,iostat=ios) nvcnf(jj)
   call read_error(ios,inline)
   do ii=nc+1,nc+nvcnf(jj)
     call cmtskp(inline)
     read(5,*,iostat=ios) nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
     call read_error(ios,inline)
   end do
 end do

! end of reading input data

 ncnf=0
 nvcnf(1)=nv
 do ii=1,nc+nv
   nacnf(ii,1)=na(ii)
   lacnf(ii,1)=la(ii)
   facnf(ii,1)=fa(ii)
 end do

 do jj=2,ncnf+1
   do ii=1,nc
     nacnf(ii,jj)=na(ii)
     lacnf(ii,jj)=la(ii)
     facnf(ii,jj)=fa(ii)
   end do
 end do

!TEMPORARY FIX
 do l1=1,lmax+1
   if(nproj(l1)==0) nproj(l1)=1
 end do

!DEFAULT SUBSTITUTONS FOR META-INCOMPATABLE DATA
 if(lloc .ne. 4) then
   nproj(lloc+1)=1
   ep(lloc+1)=0.25
   rc(5)=rc(lloc+1)
   lloc=4
   lpopt=5
 end if

 call check_data(atsym,zz,fcfact,rcfact,epsh1,epsh2,depsh,rlmax,drl,fa,facnf, &
&                rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile)

 nrl=int((rlmax/drl)-0.5d0)+1

!PWSCF wants an even number of mesh pointe
!if(trim(psfile)=='upf') then
  if(mod(nrl,2)/=0) nrl=nrl+1
!end if

!amesh=1.012d0
!amesh=1.006d0
 amesh=1.003d0

 eta = 0.01

 if(iexc==5) then
   write(6,'(a/)') 'Metagga fumctional r2scan01'
   eta = 0.01
   call r2scaninit(eta)
   call Init_GlobalConstants()
 end if

 mxprj=5

 al=dlog(amesh)

 rr1=0.0005d0/zz
 rr1=dmin1(rr1,0.0005d0/10)
 mmax=dlog(45.0d0 /rr1)/al

!calculate zion for output
 zion=zz
 do ii=1,nc
  zion=zion-fa(ii)
 end do


 allocate(rr(mmax))
 allocate(rho(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax),up(mmax),uupsa(mmax,30))
 allocate(evkb(mxprj,4), cvgplt(2,7,mxprj,4),qq(mxprj,mxprj))
 allocate(vp(mmax,5),vfull(mmax),vkb(mmax,mxprj,4),pswf(mmax,mxprj,4))
 allocate(pspwf(mmax,mxprj,4),psppwf(mmax,mxprj,4))
 allocate(vwell(mmax))
 allocate(vpuns(mmax,5))
 allocate(vo(mmax),vxc(mmax))
 allocate(vh(mmax),vtau(mmax),voep(mmax),tau(mmax),tauc(mmax))
 allocate(taumod(mmax),taumodps(mmax),vtaumodps(mmax))
 allocate(rhomod(mmax,5))
 allocate(rhoae(mmax,nv),rhomodps(mmax),tautae(mmax))
 allocate(rhov(mmax),tauv(mmax))
 allocate(rhodiff(mmax),taudiff(mmax))
 allocate(vhps(mmax),vxcps(mmax),vtaups(mmax),rhops(mmax),taups(mmax))
 allocate(npa(mxprj,6))
 allocate(epa(mxprj,6),fpa(mxprj,6))
 allocate(uua(mmax,mxprj),upa(mmax,mxprj))
 allocate(u0(mmax),u0p(mmax))
 allocate(work(mmax))
 allocate(rholoc(mmax),tauloc(mmax))

 vp(:,:)=0.0d0
 vkb(:,:,:)=0.0d0
 epa(:,:)=0.0d0
 vfull(:)=0.0d0 ; rho(:)=0.0d0 ; rhoc(:)=0.0d0

 do ii=1,mmax
  rr(ii)=rr1*exp(al*(ii-1))
 end do

 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=dmax1(rcmax,rc(l1))
 end do
 do l1=1,lmax+1
   if(rc(l1)==0.0d0) then
     rc(l1)=rcmax
   end if
 end do
  
!find log mesh point nearest input rc
 rcmax=0.0d0
 irc(:)=0
 do l1=1,lmax+1
   irc(l1)=0
   do ii=2,mmax
     if(rr(ii)>rc(l1)) then
       irc(l1)=ii
       exit
     end if
   end do
   rcmax=dmax1(rcmax,rc(l1))
 end do

 do ii=2,mmax
   if(rr(ii)>rc(lloc+1)) then
     irc(lloc+1)=ii
     exit
   end if
 end do

 ircmin=mmax
 ircmax=0
 do l1=1,lmax+1
   if(irc(l1)<ircmin) ircmin=irc(l1)
   if(irc(l1)>ircmax) ircmax=irc(l1)
 end do
!
! full potential atom solution
!
   iexct=4
   it=0
   call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&              rr,vfull,zz,mmax,iexct,etot,ierr,srel,ircmin)

   eapbe(:)=ea(:)
!
   call Init_GlobalConstants()

  rho(:)=0.0d0 ; vh(:)=0.0d0 ; vxc(:)=0.0d0 ; tau(:)=0.d0 ;  vtau(:)=0.0d0
  rhoc(:)=0.0d0; tauc(:)=0.0d0

   call sratom_m(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&              rr,vfull,vh,vxc,tauc,tau,vtau,zz,mmax,iexc, &
&              etot,dltetot,ierr,srel,ircmin)

    rhov(:)=rho(:)-rhoc(:)
    tauv(:)=tau(:)-tauc(:)

 nproj(lloc+1)=0
 rc0(:)=rc(:)

! output printing (echos input data, with all-electron eigenvalues added)

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z   nc   nv     iexc    psfile'
 write(6,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',psfile
 write(6,'(a/a,a)') '#','#   n    l    f      MGGA eval (Ha)', &
&      '         PBE         delta'

 flush(6)
 do ii=1,nc+nv
   write(6,'(2i5,f8.2,1p,2d16.6,2d12.2)') na(ii),la(ii),fa(ii),ea(ii), &
&        eapbe(ii),ea(ii)-eapbe(ii)
 end do
 flush(6)

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc(l1),ep(l1),ncon(l1),&
&        nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

!DEFAULT SUBSTITUTONS FOR META-INCOMPATABLE DATA
!defaults used unless input icmod = 5
 if(icmod /= 5) then
    icmod=5
    fcfact=3.0d0
    rcfact=0.5d0
 end if
 write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&      '# icmod, fcfact, rcfact'
 write(6,'(i5,2f10.5)') icmod,fcfact,rcfact

 write(6,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
&      '# epsh1, epsh2, depsh'
 write(6,'(3f8.2)') epsh1,epsh2,depsh

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.4)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

 write(6,'(a,f12.6)') ' amesh = ',amesh

 write(6,'(//a)') 'Reference configufation results'
 write(6,'(a,i6)') '  iterations',it
 if(it .ge. 100) then
  write(6,'(a)') 'oncvpsp: ERROR all-electron reference atom not converged'
  stop
 end if
 write(6,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot

 cvgplt(:,:,:,:)=0.0d0

! iniitialize arrays for psedo wave function metagga input
 fpa(:,:)=0.0d0
 rhops(:)=0.0d0
 taups(:)=0.0d0
 vtaups(:)=0.0d0
 zval=0.0d0

 vp(:,:)=0.0d0
!
! loop to construct pseudopotentials for all angular momenta
!
 write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

!temporarily set this to 1 so that the pseudo wave function needed for the
!local potential will be generated.  Reset after run_vkb.
 nproj(lloc+1)=1

 do l1=1,lmax+1
   ll=l1-1
   uu(:)=0.0d0; qq(:,:)=0.0d0
   iprj=0

!get principal quantum number for the highest core state for this l
   npa(1,l1)=l1
   do kk=1,nc
     if(la(kk)==l1-1) npa(1,l1)=na(kk)+1
   end do !kk

!get all-electron bound states for projectors
   if(nv/=0) then
     do kk=nc+1,nc+nv
       if(la(kk)==l1-1) then
         iprj=iprj+1
         et=ea(kk)
         call lschfb_m(na(kk),la(kk),ierr,et, &
&                      rr,vfull,vtau,uu,up,zz,mmax,mch,srel,ircmin)
         if(ierr /= 0) then
            write(6,'(/a,3i4)') 'oncvpsp-530: lschfb convergence ERROR n,l,iter=', &
&           na(ii),la(ii),it
           stop
         end if
         epa(iprj,l1)=ea(kk)
         npa(iprj,l1)=na(kk)
         fpa(iprj,l1)=fa(kk)
         uua(:,iprj)=uu(:)
         upa(:,iprj)=up(:)
       end if !la(kk)==l1-1
       if(iprj==nproj(l1)) exit
     end do !kk
   end if !nv/=0

!get all-electron well states for projectors
!if there were no valence states, use ep from input data for 1st well state
!otherwise shift up by input debl

   if(iprj==0) epa(1,l1)=ep(l1)
   if(iprj<nproj(l1))then
     do kk=1,nproj(l1)-iprj
       iprj=iprj+1
       if(iprj>1 .and. debl(l1)<=0.0d0) then
         write(6,'(a,f8.3,a/a)') 'oncvpsp: ERROR debl =',debl, 'for l=', &
&              ' ERROR not allowed with 2 or more scattering states', &
&              'program will stop'
         stop
       end if
       if(iprj>1) then
         epa(iprj,l1)=epa(iprj-1,l1)+debl(l1)
         npa(iprj,l1)=npa(iprj-1,l1)+1
       end if
        call wellstate_m(npa(iprj,l1),ll,irc,epa(iprj,l1),rr, &
&                     vfull,vxc,vtau,uu,up,zz,mmax,mch,srel,ircmin)
       uua(:,iprj)=uu(:)
       upa(:,iprj)=up(:)
     end do !kk
   end if !iprj<nproj(l1)

!get all-electron overlap matrix
   do jj=1,nproj(l1)
     do ii=1,jj
       call fpovlp(uua(1,ii),uua(1,jj),irc(l1),ll,zz,qq(ii,jj),rr,srel,mmax)
       qq(jj,ii)=qq(ii,jj)
     end do
   end do
   
   do kk=1,6
     call run_optimize_m(epa(1,l1),ll,mmax,mxprj,rr,uua,qq, &
&                    irc(l1),qcut(l1),qmsbf(l1),ncon(l1),nbas(l1),nproj(l1), &
&                    pswf(1,1,l1),pspwf(1,1,l1),psppwf(1,1,l1), &
&                    vp(1,l1),vkb(1,1,l1),vfull,cvgplt(1,1,1,l1),nderr)

     if(nderr) then
       qcut(l1)=0.95d0*qcut(l1)
       write(6,'(a,i1,a,f10.5)') ' WARNING pspot_m: first ll=',l1-1, &
&            ' PSWF has node, reducing qcut to',qcut(l1)
     else
       exit
     end if
   end do !kk

   if(nderr) then
     write(6,'(a,i1,a)') 'ERROR qcut reduction by 0.74 for ll=',l1-1, &
&          ' failed to remove node; stopping'
     stop
   end if


 end do !l1

! accumulate pseudo wavefunctions rho and tau
 rhops(:)=0.0d0
 taups(:)=0.0d0
 do l1=1,lmax+1
   do iprj=1,nproj(l1)
     if(fpa(iprj,l1)>0.0d0) then

       rhops(:)=rhops(:)+fpa(iprj,l1)*(pswf(:,iprj,l1)/rr(:))**2

       ll=l1-1
       sls=ll*(ll+1)
       taups(:)=taups(:) + 0.5d0*fpa(iprj,l1)* &
                (((pspwf(:,iprj,l1)/al - pswf(:,iprj,l1))/rr(:)**2)**2 &
&                + (sls/rr(:)**2)*(pswf(:,iprj,l1)/rr(:))**2)
      end if
   end do
 end do

!Create  new models rho and tau.
!These are based on the differences between the AE and PS rhos and taus.
!First, a fraction rcfact of the PS maxima of rho and tau are added to 
!the diffs starting at the minimum rc times (1 - (r/rcmin)**2)**4
!Then the sums are analytically continued to r=0 

 rhopsmax=0.0d0
 taupsmax=0.0d0
 do ii=mmax,1,-1
   if(rhops(ii)>rhopsmax) rhopsmax=rhops(ii)
   if(taups(ii)>taupsmax) taupsmax=taups(ii)
 end do
 

 work(:)=0.0d0

 rhodiff(:)=0.0d0
 taudiff(:)=0.0d0
 do ii=1,mmax
   xx=(1.3d0*rr(ircmax)-rr(ii))/0.3d0*rr(ircmax)
   rhodiff(ii)= fcrossover(xx)*(rho(ii)-rhops(ii)) + fcrossover(-xx)*rhoc(ii)
   taudiff(ii)= fcrossover(xx)*(tau(ii)-taups(ii)) + fcrossover(-xx)*tauc(ii)
 end do

!interval-halving search for analytically extending rhodiff to origin
!so that rholoc(1)=fctfact*rhopsmax

 irxmax=0
 do ii=ircmin,1,-1
   if(irxmax==0 .and. rr(ii)<0.8d0*rr(ircmin)) irxmax=ii
   if(rr(ii)<0.1d0*rr(ircmin)) then
     irxmin=ii
     exit
   end if
 end do

 mrho0=fcfact*rhopsmax

 do ii=1,ircmin
   xx=rr(ii)/rr(ircmin)
   rhodiff(ii)=rhodiff(ii)+rcfact*rhopsmax*(1-xx**2)**4
 end do

 irx=ircmin

 do iter=1,20
   call rtloc_m(rr,rhodiff,rholoc,0.0d0,irx,mmax,icmod)
   if(rholoc(1)<mrho0) then
     irxmax=irx
   else
     irxmin=irx
   end if
     irx=(irxmax+irxmin)/2
 end do !iter

 rhomodps(:)=rholoc(:)+rhops(:)

 rhomod(:,:)=0.0d0
 rhomod(:,1)=rholoc(:)


 call derivs(mmax, rhomod(1,1), al, rr, rhomod(1,2), rhomod(1,3), work)

!interval-halving search for analytically extending taudiff to origin
!so that tauloc(1)=fctfact*taupsmax

 irxmax=0
 do ii=ircmin,1,-1
   if(irxmax==0 .and. rr(ii)<0.8d0*rr(ircmin)) irxmax=ii
   if(rr(ii)<0.1d0*rr(ircmin)) then
     irxmin=ii
     exit
   end if
 end do

 mtau0=fcfact*taupsmax

 do ii=1,ircmin
   xx=rr(ii)/rr(ircmin)
   taudiff(ii)=taudiff(ii)+rcfact*taupsmax*(1-xx**2)**4
 end do
 tauloc(:)=taudiff(:)

 irx=ircmin

 do iter=1,20

   call rtloc_m(rr,taudiff,tauloc,0.0d0,irx,mmax,icmod)
   if(tauloc(1)<mtau0) then
     irxmax=irx
   else
     irxmin=irx
   end if
     irx=(irxmax+irxmin)/2
 end do !iter

 taumodps(:)=tauloc(:)+taups(:)

!write(6,'(/a,f8.4)') 'Model tau extension to rr=0 begins at rr = ',rr(irx)
 
 work(:)=0.0d0
 vxcps(:)=0.0d0

 call vout_m(1,rhomodps,taumodps,work,vxcps,vtaumodps,zval,sf,eeel,eexc,etot, &
&                rr,mmax,iexc)

! construct Vanderbilt / Kleinman-Bylander projectors

 write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb_m(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,mxprj, &
&                      pswf,pspwf,psppwf,vfull,vp,vtaumodps,evkb,vkb,nlim)

!restore this to its proper value
 nproj(lloc+1)=0

 deallocate(uua,upa)

! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

! null charge and eigenvalue accumulators
 uupsa(:,:)=0.0d0
 ekvnl=0.0d0
 zval=0.0d0
 rhops(:)=0.0d0
 nodes(:)=0
 taups(:)=0.0d0
 irps=0

 write(6,'(a)') 'Psueodwavefunction solutions for first PS total energy terms'
 write(6,'(a/)') '  n  l       MGGA eval      <|KE,VLOC,VKB|>     delta'
 do kk=1,nv

   et=ea(nc+kk)
   ll=la(nc+kk)
   l1=ll+1
   call lschfb_m(na(nc+kk),ll,ierr,et, &
&                rr,vfull,vtau,uu,up,zz,mmax,mch,srel,ircmin)
   if(ierr /= 0) then
     write(6,'(/a,3i4)') 'oncvpsp-768: lschfb_m convergence ERROR n,l,iter=', &
&           na(ii),la(ii),it
     stop
   end if

   sls=ll*(ll+1)
   tautae(:)=tautae(:) + 0.5d0*fa(nc+kk)*(((up(:)/al - uu(:))/rr(:)**2)**2 &
&      + (sls/rr(:)**2)*(uu(:)/rr(:))**2)

   emax=0.75d0*et
   emin=1.25d0*et
   call lschvkbb_m(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax,ebar, &
&                rr,vp(1,lloc+1),vtaumodps,vkb(1,1,l1),evkb(1,l1), &
&                uu,up,mmax,mch)

   if(ierr/=0) then 
     write(6,'(a,3i4)') 'oncvpsp-752: lschvkbb_m ERROR',ll+nodes(l1)+1,ll,ierr
     flush(6)
     stop
   end if

  write(6,'(2i3,1p,2e18.6,e14.2)') ll+nodes(l1)+1,ll,et,ebar,ebar-et
!   write(6,'(2i5,f8.2,1p,2d16.6,2d12.2)') na(ii),la(ii),fa(ii),ea(ii), &
!&        et,ebar,ebar-et

! save valence pseudo wave functions for upfout
   do ii=1,mmax
     if(uu(ii)==0.0d0) exit
   
     uupsa(ii,kk)=uu(ii)
     rhops(ii)=rhops(ii)+fa(nc+kk)*(uu(ii)/rr(ii))**2

     sls=ll*(ll+1)
     taups(ii)=taups(ii) + 0.5d0*fa(nc+kk)*(((up(ii)/al - uu(ii))/rr(ii)**2)**2 &
&              + (sls/rr(ii)**2)*(uu(ii)/rr(ii))**2)
   end do

   ekvnl=ekvnl+fa(nc+kk)*ebar

   zval=zval+fa(nc+kk)
   nodes(l1)=nodes(l1)+1
   irps=max(irps,irc(l1))
 end do !kk


 rhomodps(:)=rhops(:)+rhomod(:,1)
 taumodps(:)=taups(:)+tauloc(:)
 taumod(:)=tauloc(:)

 vhps(:)=0.0d0 ;  vxcps(:)=0.0d0 ; vtaumodps(:)=0.0d0
 
 sf=0.0d0

 call vout_m(0,rhops,taumodps,vhps,vxcps,vtaumodps,zval,sf,eeel,eexc,etot, &
&                rr,mmax,iexc)
 work(:)=0.0d0
 call vout_m(1,rhomodps,taumodps,work,vxcps,vtaumodps,zval,sf,eeel,eexc,etot, &
&                rr,mmax,iexc)

! total energy output - note "ekvnl" is the sum of the kinetic, local
! potential, and non-local potental expectation values for the occupied
! valence states.

 epstot= ekvnl + eexc - 0.5d0*eeel
 write(6,'(/a)') 'Exchange-correlation and electron-electron Coulomb terms added'
 write(6,'(a,1p,e18.8/)') 'Pseudoatom total energy', epstot


 call run_diag_m(lmax,npa,epa,lloc,irc, vkb,evkb,nproj,rr, &
&                vfull,vtau,vp,vtaumodps,zz,mmax,mxprj,srel,ircmin)

 call run_ghosts_m(lmax,la,ea,nc,nv,lloc,irc,qmsbf, &
&                    vkb,evkb,nproj,rr,vp,vtaumodps,mmax,mxprj)
 flush 6



! unscreen local and semi-local potentials local vtau

 do l1=1,max(lmax+1,lloc+1)
   vpuns(:,l1)=vp(:,l1)-vhps(:)-vxcps(:)
 end do

!fix unscreening error due to greater range of all-electron charge
 do ii=mmax,1,-1
   if(rhops(ii)==0.0d0) then 
     do l1=1,max(lmax+1,lloc+1)
       vpuns(ii,l1)=-zion/rr(ii)
     end do
   else
     exit
   end if
 end do

 rho(:)=rhops(:)


 call run_plot_m(lmax,npa,epa,lloc,irc, &
&                vkb,evkb,nproj,rr,vfull,vp,vpuns, &
&                zz,mmax,mxprj,drl,nrl,rhops,rhoc,rhomod,srel,cvgplt, &
&                tauc,taups,taumodps,vtau,vtaumodps,ircmin)


 call run_phsft_m(lmax,lloc,nproj,epa,epsh1,epsh2,depsh,vkb,evkb, &
&               rr,vfull,vtau,vp,vtaumodps,zz,mmax,mxprj,irc,srel,ircmin)
 
 call gnu_script_m(epa,evkb,lmax,lloc,mxprj,nproj)


!Exit if no psfile output is wanted yet

 if(trim(psfile)=='none') stop

 if(trim(psfile)=='psp8') then

  rhomod(:,4)=taups(:)
  rhomod(:,5)=taumod(:)
  call linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             rhov,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile)
 end if

 if(trim(psfile)=='upf') then
  call upfout_m(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             taups,taumod, &
&             zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym,epstot, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile,uupsa,ea,srel)
 end if

 if(trim(psfile)=='psml') then

  irct=ircmax

   write(6,'(/a)',advance="no") '<?xml version="1.0" encoding="UTF-8" ?>'

   call get_ec_hints(cvgplt)

   call psmlout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&               taups,taumod, &
&               irct, srel, &
&               zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym,epstot, &
&               na,la,nc,nv, &
&               fa,epa, &
&               psfile,uupsa,ea, &
&               lpopt,ncnf,nvcnf,nacnf,lacnf,ncon,nbas, &
&               depsh,dvloc0,epsh1,epsh2,fcfact,rlmax,rcfact, &
&               ep,qcut,debl,facnf)

 end if

 stop
 end program oncvpsp_m


 subroutine cmtskp(inline)
! skips lines of standard input (file 5) whose first character is #
 implicit none

!In/Out variable
 integer :: inline

!Local variable
 character*1 tst

 tst='#'
 do while (tst=='#')
  read(5,*) tst
  inline=inline+1
 end do
 backspace(5)
 inline=inline-1

 return
 end subroutine cmtskp

 subroutine read_error(ios,inline)
! report data read error and stop
 implicit none

!Input variables
 integer :: ios,inline

 inline=inline+1
 if(ios/=0) then
  write(6,'(a,i4)') 'Read ERROR, input data file line',inline
  write(6,'(a)') 'Program will stop'
  stop
 end if

 return
 end subroutine read_error
