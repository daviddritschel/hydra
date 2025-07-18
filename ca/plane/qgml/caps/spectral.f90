module spectral

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! This module contains all subroutines related to spectral inversion
! and differentiation.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 !Import constants and parameters:
use constants
 !Import FFT module:
use sta2dfft

 !Declarations:
implicit none

 !Common arrays, constants:
double precision:: green(ng,ng,nz),rksq(ng,ng)
double precision:: diss(ng,ng),filt(ng,ng)
double precision:: bflo(ng,ng),bfhi(ng,ng)
double precision:: bety(ng,ng),qb(ng,ng)
double precision:: sfwind(ng,ng)
double precision:: srwfm

 !Vertical structure:
double precision:: hhat(nz),lambda(nz),kdsq(nz),kk0(nz),kkm(nz)
double precision:: vl2m(nz,nz),vm2l(nz,nz)

 !For 2D FFTs:
double precision:: hrkx(ng),hrky(ng),rk(ng)
double precision:: xtrig(2*ng),ytrig(2*ng)
integer:: xfactors(5),yfactors(5)

!==========================================================================!
! From main code: call init_spectral                to initialise          !
! then            call main_invert(qq,pp,uu,vv)     to perform inversion   !
!==========================================================================!

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================

subroutine init_spectral
! Initialises this module

implicit none

!Local variables:
double precision:: wkp(ng,ng),daf(ng),yg(ng)
double precision:: fac,rkmax,td
integer:: kx,ky,k,kc,ix,iy,iz,m

!----------------------------------------------------------------------
 !Set up 2D FFTs:
call init2dfft(ng,ng,twopi,twopi,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define wavenumbers and filtered wavenumbers:
rk(1)=zero
do k=1,ng/2-1
   rk(k+1)   =hrkx(2*k)
   rk(ng+1-k)=hrkx(2*k)
enddo
rk(ng/2+1)=hrkx(ng)

 !Define de-aliasing filter (2/3 rule):
rkmax=dble(ng/2)
daf(1)=one
do k=2,ng
   if (rk(k) < f23*rkmax) then
      daf(k)=one
   else
      daf(k)=zero
   endif
enddo

 !Define filter and squared horizontal wavenumber:
do ky=1,ng
   do kx=1,ng
      filt(kx,ky)=daf(kx)*daf(ky)
      rksq(kx,ky)=rk(kx)**2+rk(ky)**2
   enddo
enddo

 !Define Butterworth low-pass (F) & high-pass (1-F) filters:
fac=9.d0/rkmax**2
bflo=filt/(one+(fac*rksq)**2)
bfhi=filt*(one-bflo)

 !Define horizontal hyperdiffusion spectral operator:
fac=one/rkmax**2
diss=cdamp*(fac*rksq)**nnu

!--------------------------------------------------------------------
 !Define y grid lines for use below:
do iy=1,ng
   yg(iy)=gl*dble(iy-1)-pi
enddo

!--------------------------------------------------------------------
 !Read vertical structure and mode files (generated by vertical.f90):
open(60,file='vertical.asc',status='old')
do iz=1,nz
   read(60,*) hhat(iz),kdsq(iz)
enddo
close(60)
 !hhat = mean layer depth / total mean depth (sum(hhat) = 1).
 !kdsq = f^2/(b'*H) where f = Coriolis frequency, b' = buoyancy
 !       difference between layer iz and iz+1, H = total mean depth.
 !Note: kdsq(nz) is unused.

 !Ensure hhat sums to 1 - this is essential:
fac=one/sum(hhat)
hhat=fac*hhat

!-----------------------------------------------------------------
 !For computing relative vorticity and implementing Ekman drag:
do iz=1,nz-1
   kk0(iz)=kdsq(iz)/hhat(iz)
   kkm(iz+1)=kdsq(iz)/hhat(iz+1)
enddo

!-----------------------------------------------------------------
 !Read vertical eigenvalues and eigenmodes:
open(60,file='modes.asc',status='old')
do m=1,nz
   read(60,*) lambda(m)
enddo
 !Note, lambda(m) corresponds to lambda_m in the notes.
do m=1,nz
   do iz=1,nz
      read(60,*) vl2m(iz,m),vm2l(iz,m)
   enddo
enddo
close(60)
 !vl2m converts layer quantities to  mode quantities
 !vm2l converts  mode quantities to layer quantities

!--------------------------------------------------------------------
 !Define Green function for inverting Lap(psi)-lambda*psi = S for psi:
rksq(1,1)=one
green(:,:,1)=-one/rksq
rksq(1,1)=zero
green(1,1,1)=zero
 !The above is for vertical mode m = 1, for which lambda(m) = 0.
 !The 0 wavevector is removed; green cannot be used for this.

 !Define Green function for m > 1, for which lambda(m) > 0:
do m=2,nz
   green(:,:,m)=-one/(lambda(m)+rksq)
enddo

!------------------------------------------------------------------
 !If bathymetry eta_b is present, read in qb = f_0*eta_b/H_{nz}:
if (bath) then
   open(11,file='bath.r8',form='unformatted', &
        access='direct',status='old',recl=2*nhbytes)
   read(11,rec=1) td,qb
   close(11)
endif

!------------------------------------------------------------------
 !If wind-stress forcing is present, define spectral function to
 !add to PV tendency in the uppermost layer, i.e.
 !Dq_1/Dt = fwind*sin(y):
if (wind) then
   do ix=1,ng
      wkp(:,ix)=fwind*sin(yg)
   enddo
   !convert to spectral space:
   call ptospc(ng,ng,wkp,sfwind,xfactors,yfactors,xtrig,ytrig)
   !sfwind is added to (spectral) Dq_1/Dt in evolution.f90.
   sfwind(1,1)=zero
   ! The mean tendency is here zero; for a more general wind-stress
   ! forcing function, it may be desirable to keep the this zero to
   ! avoid a build up of mean PV in the uppermost layer, which will
   ! induce a growth in barotropic circulation in the domain.
endif
 !Note: we have to transform a sine in y function to a cosine
 !      series above to be compatable with the representation
 !      of PV as a cosine series in y.

!------------------------------------------------------------------
 !Define beta*y for use in PV inversion --- see main_invert below:
do ix=1,ng
   bety(:,ix)=beta*yg
enddo
 !This is written as a 2D array for simpler use in the main code.

 !Define maximum Rossby wave frequency (for time stepping):
if (beta > 0) then
   srwfm=f12*beta/rk(1)
else
   srwfm=small
endif

return
end subroutine init_spectral

!=======================================================================

subroutine main_invert(qq,pp,uu,vv)
 !Given the PV qq in all layers, this routine computes the
 !streamfunction pp and velocity field (uu,vv) in all layers.

 !Declarations:
implicit none

 !Passed arrays:
double precision:: qq(ng,ng,nz)
double precision:: pp(ng,ng,nz),uu(ng,ng,nz),vv(ng,ng,nz)

 !Local quantities:
double precision:: qm(ng,ng,nz),pm(ng,ng)
double precision:: wks(ng,ng)
integer:: iz,m

!------------------------------------------------------------------
 !Take away bathymetry contribution to bottom layer PV (iz = nz):
if (bath) qq(:,:,nz)=qq(:,:,nz)-qb

 !Project layer PV (in qq) onto vertical modes (as qm):
qm=zero
do m=1,nz
   do iz=1,nz
      qm(:,:,m)=qm(:,:,m)+vl2m(iz,m)*qq(:,:,iz)
   enddo
enddo

 !Remove beta*y (if it exists) from the barotropic mode PV:
if (beffect) qm(:,:,1)=qm(:,:,1)-bety
 !Note, sum_{iz} vl2m(iz,1) = 1 (see init/vertical.f90)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !For each vertical mode m, find the associated streamfunction pm:
do m=1,nz

    !Store mode PV in pm temporarily:
   pm=qm(:,:,m)

    !FFT mode PV to spectral space as wks temporarily:
   call ptospc(ng,ng,pm,wks,xfactors,yfactors,xtrig,ytrig)

    !Invert to get streamfunction - store in pp (spectral):
   pp(:,:,m)=green(:,:,m)*wks
    !See init_spectral above for definition of green.

enddo
 !End of loop over vertical modes.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 !Project modal streamfunction in pp onto vertical layers as qm:
qm=zero
do m=1,nz
   do iz=1,nz
      qm(:,:,iz)=qm(:,:,iz)+vm2l(iz,m)*pp(:,:,m)
   enddo
enddo

 !Compute velocity field spectrally:
do iz=1,nz
   call xderiv(ng,ng,hrkx,qm(:,:,iz),wks)
   call spctop(ng,ng,wks,pm,xfactors,yfactors,xtrig,ytrig)
   vv(:,:,iz)= pm !u in physical space

   call yderiv(ng,ng,hrky,qm(:,:,iz),wks)
   call spctop(ng,ng,wks,pm,xfactors,yfactors,xtrig,ytrig)
   uu(:,:,iz)=-pm !v in physical space

   wks=qm(:,:,iz)
   call spctop(ng,ng,wks,pm,xfactors,yfactors,xtrig,ytrig)
   pp(:,:,iz)= pm !psi in physical space
enddo

return
end subroutine main_invert

!=======================================================================

subroutine vorticity(qq,pp,zz)

! This routine computes the relative vorticity, given the current total
! PV field qq and streamfunction pp, i.e. following a call to inversion.

implicit none

! Passed arrays:
double precision:: qq(ng,ng,nz),pp(ng,ng,nz),zz(ng,ng,nz)

! Local variable:
integer:: iz

!-------------------------------------------------------------
! Top layer (iz = 1):
zz(:,:,1)=qq(:,:,1)-bety+kk0(1)*(pp(:,:,1)-pp(:,:,2))
! bety = beta*y here

! Intermediate layers:
do iz=2,nz-1
   zz(:,:,iz)=qq(:,:,iz)-bety+kkm(iz)*(pp(:,:,iz)-pp(:,:,iz-1))+ &
                              kk0(iz)*(pp(:,:,iz)-pp(:,:,iz+1))
enddo

! Bottom layer (iz = nz):
zz(:,:,nz)=qq(:,:,nz)-bety+kkm(nz)*(pp(:,:,nz)-pp(:,:,nz-1))
! Remove PV due to bathymetry if present:
if (bath) zz(:,:,nz)=zz(:,:,nz)-qb

! Note, kkm(iz) = f^2/(H_iz b'_{iz-1}) for iz = 2 to nz
! while kk0(iz) = f^2/(H_iz b'_{iz}) for iz = 1 to nz-1
! --- see init_spectral above.

return
end subroutine vorticity

!=======================================================================

subroutine gradient(fs,fx,fy)
! Computes the gradient of a spectral field fs, represented as a cosine
! series in y, and returns the gradient in physical space as (fx,fy).

implicit none

 !Passed arrays:
double precision:: fs(ng,ng,nz)
double precision:: fx(ng,ng,nz),fy(ng,ng,nz)

 !Local quantities:
double precision:: wks(ng,ng)
integer:: iz

!-------------------------------------------------------------------
 !Loop over layers:
do iz=1,nz
   call xderiv(ng,ng,hrkx,fs(:,:,iz),wks)
   call spctop(ng,ng,wks,fx(:,:,iz),xfactors,yfactors,xtrig,ytrig)

   call yderiv(ng,ng,hrky,fs(:,:,iz),wks)
   call spctop(ng,ng,wks,fy(:,:,iz),xfactors,yfactors,xtrig,ytrig)
enddo

return
end subroutine gradient

!=======================================================================

end module spectral
