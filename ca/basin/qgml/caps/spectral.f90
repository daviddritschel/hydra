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
double precision:: green(0:nx,0:ny,nz),rksq(0:nx,0:ny)
double precision:: decx(nxm1,nym1,nz),decy(nym1,nxm1,nz)
double precision:: diss(0:nx,0:ny),filt(0:nx,0:ny)
double precision:: bflo(0:nx,0:ny),bfhi(0:nx,0:ny)
double precision:: xh0(0:nx,nz),xh1(0:nx,nz)
double precision:: yh0(0:ny,nz),yh1(0:ny,nz)
double precision:: psibar(0:ny,0:nx),danorm(0:ny,0:nx)
double precision:: bety(0:ny,0:nx),qb(0:ny,0:nx)
double precision:: sfwind(0:nx,0:ny)
double precision:: srwfm

 !Vertical structure:
double precision:: hhat(nz),lambda(nz),kdsq(nz),kk0(nz),kkm(nz)
double precision:: vl2m(nz,nz),vm2l(nz,nz)

double precision:: rkx(0:nx),rky(0:ny)
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

!==========================================================================!
! From main code: call init_spectral                to initialise          !
! then            call main_invert(qq,pp,uu,vv)     to perform inversion   !
!==========================================================================!

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine init_spectral

 !Declarations:
implicit none

 !Local variables:
double precision:: wkp(0:ny,0:nx),skx(nx),sky(ny)
double precision:: dafx(0:nx),dafy(0:ny)
double precision:: scx,rkxmax,scy,rkymax
double precision:: td,fac,div,argm,argp
double precision:: rkmsi
integer:: ix,iy,iz,kx,ky,m,k

!--------------------------------------------------------------------
 !Set up 2D FFTs and other commonly-used spectral arrays:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,skx,sky)

 !Define x wavenumbers:
scx=pi/ellx
rkxmax=scx*dble(nx)
do kx=0,nx
   rkx(kx)=scx*dble(kx)
enddo

 !Define y wavenumbers:
scy=pi/elly
rkymax=scy*dble(ny)
do ky=0,ny
   rky(ky)=scy*dble(ky)
enddo

 !Define de-aliasing filter (2/3 rule):
dafx(0)=one
do kx=1,nx
   if (rkx(kx) < f23*rkxmax) then
      dafx(kx)=one
   else
      dafx(kx)=zero
   endif
enddo

dafy(0)=one
do ky=1,ny
   if (rky(ky) < f23*rkymax) then
      dafy(ky)=one
   else
      dafy(ky)=zero
   endif
enddo

do ky=0,ny
   do kx=0,nx
      filt(kx,ky)=dafx(kx)*dafy(ky)
   enddo
enddo

 !Define squared horizontal wavenumber:
do ky=0,ny
   do kx=0,nx
      rksq(kx,ky)=rkx(kx)**2+rky(ky)**2
   enddo
enddo

 !Define Butterworth low-pass (F) & high-pass (1-F) filters:
fac=18.d0/(rkxmax**2+rkymax**2)
bflo=filt/(one+(fac*rksq)**2)
bfhi=filt*(one-bflo)

 !Define horizontal hyperdiffusion spectral operator:
rkmsi=one/max(rkxmax**2,rkymax**2)
diss=cdamp*(rkmsi*rksq)**nnu

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
fac=1.d0/sum(hhat)
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
 !Define Green function for inverting Lap(psi)-lambda(m)*psi = S for psi:
rksq(0,0)=one
green(:,:,1)=-one/rksq
rksq(0,0)=zero
green(0,0,1)=zero
 !The above is for vertical mode m = 1, for which lambda(m) = 0.

 !Define Green function for m > 1, for which lambda(m) > 0:
do m=2,nz
   green(:,:,m)=-one/(lambda(m)+rksq)
enddo

!--------------------------------------------------------------------
 !Define functions for removing corner values of the streamfunction:
fac=one/dble(nx)
do ix=0,nx
   xh1(ix,1)=fac*dble(ix)
   xh0(ix,1)=one-xh1(ix,1)
enddo
fac=one/dble(ny)
do iy=0,ny
   yh1(iy,1)=fac*dble(iy)
   yh0(iy,1)=one-yh1(iy,1)
enddo
 !The above is for vertical mode m = 1 for which lambda(m) = 0.

 !Define corresponding functions for m > 1:
do m=2,nz
   fac=ellx*elly*sqrt(lambda(m)/(ellx**2+elly**2))
   div=one/sinh(fac)
   do ix=0,nx
      xh1(ix,m)=div*sinh(fac*xh1(ix,1))
      xh0(ix,m)=div*sinh(fac*xh0(ix,1))
   enddo
   do iy=0,ny
      yh1(iy,m)=div*sinh(fac*yh1(iy,1))
      yh0(iy,m)=div*sinh(fac*yh0(iy,1))
   enddo
enddo
 !Note: xh0 = 1 at x = x_min & xh0 = 0 at x = x_max
 !while xh1 = 0 at x = x_min & xh1 = 1 at x = x_max
 !The same holds replacing x by y.  The products xhj*yhk
 !for j = 0,1 and k = 0,1 are solutions of Helmholtz' equation
 !Lap(phi)-lambda(m)*phi = 0 for each mode m.  These products
 !vanish along two edges of the rectangular domain and equal
 !unity only on the opposite corner.  A linear combination
 !of xh0*yh0, xh0*yh1, xh1*yh0 and xh1*yh1 is used to remove
 !corner values of the streamfunction, ensuring that the
 !remaining variation along each edge can be represented as
 !a sine series in x or y.  See Dritschel, Dritschel & Carr,
 !JCP-X 17 (2023), 100129, available as full open access at
 !https://doi.org/10.1016/j.jcpx.2023.100129

!--------------------------------------------------------------------
 !Hyperbolic functions needed as solutions of Helmholtz's equation:
do m=1,nz
   do ky=1,nym1
      fac=sqrt(lambda(m)+rky(ky)**2)*ellx
      div=one/(one-exp(-two*fac))
      do ix=1,nxm1
         argm=fac*(one-xh1(ix,1))
         argp=fac*(one+xh1(ix,1))
         decx(ix,ky,m)=(exp(-argm)-exp(-argp))*div
      enddo
   enddo
   do kx=1,nxm1
      fac=sqrt(lambda(m)+rkx(kx)**2)*elly
      div=one/(one-exp(-two*fac))
      do iy=1,nym1
         argm=fac*(one-yh1(iy,1))
         argp=fac*(one+yh1(iy,1))
         decy(iy,kx,m)=(exp(-argm)-exp(-argp))*div
      enddo
   enddo
enddo
 !The above functions allow one to remove the sinusoidally-varying
 !part of the streamfunction along each edge of the domain, after
 !the corner values have been removed as discussed in the comments
 !just above. The reference above discusses how these functions are
 !used to obtain a final solution in which the streamfunction is
 !constant along the entire domain boundary, while along the PV
 !to vary arbitrarily.

!--------------------------------------------------------------------
 !Define part of streamfunction proportional to the mean vorticity
 !of the barotropic mode (m = 1):
do ix=0,nx
   do iy=0,ny
      psibar(iy,ix)=-f14*(ellx**2*xh0(ix,1)*xh1(ix,1)+ &
                          elly**2*yh0(iy,1)*yh1(iy,1))
   enddo
enddo

 !Define area weights needed to compute a horizontal average:
danorm(1:nym1,1:nxm1)=dsumi
danorm(1:nym1, 0)=f12*dsumi
danorm(1:nym1,nx)=f12*dsumi
danorm( 0,1:nxm1)=f12*dsumi
danorm(ny,1:nxm1)=f12*dsumi
danorm( 0, 0)=f14*dsumi
danorm( 0,nx)=f14*dsumi
danorm(ny, 0)=f14*dsumi
danorm(ny,nx)=f14*dsumi
 !The average of a field f is computed using sum(f*danorm).

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
 !Dq_1/Dt = fwind*sin(2*pi*(y-y_min)/(y_max-y_min)):
if (wind) then
   do ix=0,nx
      wkp(:,ix)=fwind*sin(twopi*yh1(:,1))
   enddo
   !yh1 = (y-y_min)/(y_max-y_min) is defined in init_spectral
   call ptospc_cc(nx,ny,wkp,sfwind,xfactors,yfactors,xtrig,ytrig)
   !sfwind is added to (spectral) Dq_1/Dt in evolution.f90.
   sfwind(0,0)=zero
   ! The mean tendency is here zero; for a more general wind-stress
   ! forcing function, it may be desirable to keep the this zero to
   ! avoid a build up of mean PV in the uppermost layer, which will
   ! induce a growth in barotropic circulation in the domain.
endif
 !Note: we have to transform a sine in y function to a cosine
 !      series above to be compatable with the representation
 !      of PV as a double cosine series.

!------------------------------------------------------------------
 !Define beta*y for use in PV inversion --- see main_invert below:
do ix=0,nx
   do iy=0,ny
      bety(iy,ix)=beta*(ymin+gly*dble(iy))
   enddo
enddo
 !This is written as a 2D array for simpler use in the main code

 !Define maximum Rossby wave frequency (for time stepping):
if (beta > 0) then
   srwfm=beta*max(rkx(1),rky(1))/(rkx(1)**2+rky(1)**2)
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
double precision:: qq(0:ny,0:nx,nz)
double precision:: pp(0:ny,0:nx,nz),uu(0:ny,0:nx,nz),vv(0:ny,0:nx,nz)

 !Local quantities:
double precision:: qm(0:ny,0:nx,nz),pm(0:ny,0:nx)
double precision:: wks(0:nx,0:ny)
double precision:: pbot(nx),ptop(nx),cppy(nym1,nx)
double precision:: plft(ny),prgt(ny),cppx(nxm1,ny)
double precision:: pmbar(nz),pmha(nz),ppha(nz)
double precision:: sw00,sw01,sw10,sw11
integer:: ix,iy,iz,kx,ky,m

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

 !Remove beta*y from the barotropic mode PV:
qm(:,:,1)=qm(:,:,1)-bety
 !Note, sum_{iz} vl2m(iz,1) = 1 (see init/vertical.f90)

!------------------------------------------------------------------
 !For each vertical mode m, find the associated streamfunction pm:
do m=1,nz

    !Store mean source term in pmbar for computing mean pm below:
   pmbar(m)=sum(qm(:,:,m)*danorm)
    !Here, danorm = dx * dy / (L_x * L_y) essentially.

   if (m == 1) then
       !For the barotropic mode, remove mean vorticity from qm(:,:,1):
      pm=qm(:,:,1)-pmbar(1)
   else
       !For all other vertical modes, store modal PV in 2D pm array:
      pm=qm(:,:,m)
   endif
 
    !FFT mode PV and invert to get uncorrected streamfunction pm:
   call ptospc_cc(nx,ny,pm,wks,xfactors,yfactors,xtrig,ytrig)
    !ptospc_cc allows the PV to vary generally along the boundaries.
   wks=green(:,:,m)*wks
    !See init_spectral above for definition of green.
   call spctop_cc(nx,ny,wks,pm,xfactors,yfactors,xtrig,ytrig)

    !For the barotropic mode, add part due to mean vorticity:
   if (m == 1) pm=pm+pmbar(1)*psibar
    !psibar is a quadratic function of x & y defined in init_spectral.

    !Remove homogeneous solutions of Helmholtz' equation so that pm
    !vanishes at the four corners of the x-y domain (for each mode m):
   sw00=pm(0,0)
   sw10=pm(ny,0)
   sw01=pm(0,nx)
   sw11=pm(ny,nx)
   do ix=0,nx
      do iy=0,ny
         pm(iy,ix)=pm(iy,ix) &
              -(sw00*xh0(ix,m)+sw01*xh1(ix,m))*yh0(iy,m) &
              -(sw10*xh0(ix,m)+sw11*xh1(ix,m))*yh1(iy,m)
      enddo
   enddo
    !See init_spectral above for definitions of xh0, xh1, yh0 & yh1.

    !Do a 1D sine transform of pm at y = ymin and ymax and obtain
    !the interior field (cppy) that must be subtracted to give pm = 0
    !at y = ymin and ymax:
   do ix=1,nxm1
      pbot(ix)=pm(0,ix)
      ptop(ix)=pm(ny,ix)
   enddo
   call dst(1,nx,pbot,xtrig,xfactors)
   call dst(1,nx,ptop,xtrig,xfactors)

    !Define the interior semi-spectral field:
   do kx=1,nxm1
      do iy=1,nym1
         cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx,m)+ptop(kx)*decy(iy,kx,m)
      enddo
   enddo
    !Invert using a sine transform for all interior y grid lines:
   call dst(nym1,nx,cppy,xtrig,xfactors)

    !Do a 1D sine transform of pm at x = xmin and xmax and obtain
    !the interior field (cppx) that must be subtracted to give pm = 0
    !at x = xmin and xmax:
   do iy=1,nym1
      plft(iy)=pm(iy,0)
      prgt(iy)=pm(iy,nx)
   enddo
   call dst(1,ny,plft,ytrig,yfactors)
   call dst(1,ny,prgt,ytrig,yfactors)

    !Define the interior semi-spectral field:
   do ky=1,nym1
      do ix=1,nxm1
         cppx(ix,ky)=plft(ky)*decx(nx-ix,ky,m)+prgt(ky)*decx(ix,ky,m)
      enddo
   enddo
    !Invert using a sine transform for all interior x grid lines:
   call dst(nxm1,ny,cppx,ytrig,yfactors)
    !See init_spectral above for definitions of decx and decy.

    !Remove cppx and cppy to obtain the final streamfunction pm
    !which, by construction, vanishes along the edge of the domain:
   pm(:,0 )=zero
   pm(:,nx)=zero
   pm(0, :)=zero
   pm(ny,:)=zero
   do ix=1,nxm1
      do iy=1,nym1
         pm(iy,ix)=pm(iy,ix)-cppx(ix,iy)-cppy(iy,ix)
      enddo
   enddo

    !Store streamfunction for this mode in the 3D array qm:
   qm(:,:,m)=pm
enddo

!------------------------------------------------------------------
 !Project modal streamfunction in qm onto vertical layers as pp:
pp=zero
do m=1,nz
   do iz=1,nz
      pp(:,:,iz)=pp(:,:,iz)+vm2l(iz,m)*qm(:,:,m)
   enddo
enddo

!------------------------------------------------------------------
 !Compute the velocity field (uu,vv) from pp:
call velocity(pp,uu,vv)

!------------------------------------------------------------------
 !Commpute right-hand side of system A pavg = zavg - qbar:
 !Domain mean relative vorticity is the circulation / domain area:
do iz=1,nz
   ppha(iz)=(gly*sum(vv(:,nx,iz)-vv(:,0,iz)) &
            -glx*sum(uu(ny,:,iz)-uu(0,:,iz)))/domarea &
                -sum(qq(:,:,iz)*danorm)
enddo

 !Solve A Psi_avg = zeta_avg - q_avg after projection onto
 !vertical modes; since A is singular (the barotropic mode
 !has a zero Rossby deformation wavenumber), we can only do
 !this for modes m > 1. Without loss of generality, we take
 !the barotropic layer-mean streamfunction to be zero:
pmha(1)=zero
do m=2,nz
   pmha(m)=sum(vl2m(:,m)*ppha)/lambda(m)
enddo

 !Project pmha back to layers as ppha:
ppha=zero
do m=2,nz
   ppha=ppha+vm2l(:,m)*pmha(m)
enddo

 !Restore correct mean streamfunction values:
do iz=1,nz
   ppha(iz)=ppha(iz)-sum(pp(:,:,iz)*danorm)
   pp(:,:,iz)=pp(:,:,iz)+ppha(iz)
enddo

 !Restore original PV in the bottom layer (iz = nz):
if (bath) qq(:,:,nz)=qq(:,:,nz)+qb

return
end subroutine main_invert

!=======================================================================

subroutine velocity(pp,uu,vv)
! Computes the velocity components uu & vv from the streamfunction
! pp via uu = -d(pp)/dy and vv = d(pp)/dx (done in spectral space).

! pp, uu & vv are all in physical space and include the domain edges.

implicit none

 !Passed arrays:
double precision:: pp(0:ny,0:nx,nz),uu(0:ny,0:nx,nz),vv(0:ny,0:nx,nz)

 !Local quantities:
double precision:: ppi(ny,nx),pps(nx,ny)
double precision:: ppx(0:nx,ny),vvi(ny,0:nx)
double precision:: ppy(nx,0:ny),uui(0:ny,nx)
integer:: ix,iy,iz

!-------------------------------------------------------------------
 !Loop over layers:
do iz=1,nz
    !Copy non-zero interior values of pp into 2D array ppi:
   ppi=pp(1:ny,1:nx,iz)

    !Transform ppi to spectral space:
   call ptospc_ss(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)

    !Compute d(ppi)/dx = ppx spectrally:
   call xderiv_ss(nx,ny,rkx(1:nx),pps,ppx)

    !Transform ppx back to physical space as vvi:
   call spctop_cs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)

    !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
   do ix=0,nx
      vv(0,ix,iz)=zero
      do iy=1,nym1
         vv(iy,ix,iz)=vvi(iy,ix)
      enddo
      vv(ny,ix,iz)=zero
   enddo

    !Compute d(ppi)/dy = ppy spectrally:
   call yderiv_ss(nx,ny,rky(1:ny),pps,ppy)

    !Transform ppy back to physical space as uui:
   call spctop_sc(nx,ny,ppy,uui,xfactors,yfactors,xtrig,ytrig)

    !Copy -uui into uu and add on zero edge values at ix = 0 & nx:
   uu(:, 0,iz)=zero
   do ix=1,nxm1
      uu(:,ix,iz)=-uui(:,ix)
   enddo
   uu(:,nx,iz)=zero
enddo

return
end subroutine velocity

!=======================================================================

subroutine vorticity(qq,pp,zz)

! This routine computes the relative vorticity, given the current total
! PV field qq and streamfunction pp, i.e. following a call to inversion.

implicit none

! Passed arrays:
double precision:: qq(0:ny,0:nx,nz),pp(0:ny,0:nx,nz),zz(0:ny,0:nx,nz)

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
! series in x & y, and returns the gradient in physical space as (fx,fy).

implicit none

 !Passed arrays:
double precision:: fs(0:nx,0:ny,nz)
double precision:: fx(0:ny,0:nx,nz),fy(0:ny,0:nx,nz)

 !Local quantities:
double precision:: dsx(nx,0:ny),dpx(0:ny,nx)
double precision:: dsy(0:nx,ny),dpy(ny,0:nx)
integer:: iz

!-------------------------------------------------------------------
 !Loop over layers:
do iz=1,nz
    !Compute df/dx spectrally:
   call xderiv_cc(nx,ny,rkx(1:nx),fs(:,:,iz),dsx)
    !Convert to physical space:
   call spctop_sc(nx,ny,dsx,dpx,xfactors,yfactors,xtrig,ytrig)
    !Store in fx and insert zero boundary values at xmin & xmax:
   fx(:,0,iz)=zero
   fx(:,1:nxm1,iz)=dpx(:,1:nxm1)
   fx(:,nx,iz)=zero

    !Compute df/dy spectrally:
   call yderiv_cc(nx,ny,rky(1:ny),fs(:,:,iz),dsy)
    !Convert to physical space:
   call spctop_cs(nx,ny,dsy,dpy,xfactors,yfactors,xtrig,ytrig)
    !Store in fy and insert zero boundary values at ymin & ymax:
   fy(0,:,iz)=zero
   fy(1:nym1,:,iz)=dpy(1:nym1,:)
   fy(ny,:,iz)=zero
enddo

return
end subroutine gradient

!=======================================================================

end module spectral
