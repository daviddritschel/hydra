module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: rksq(ng,ng),filt(ng,ng),c2g2(ng,ng),helm(ng,ng)
double precision:: rlap(ng,ng),adop(ng,ng),pope(ng,ng),opak(ng,ng)
double precision:: simp(ng,ng),pdis(ng,ng),diss(ng,ng),rdis(ng,ng)

 !For 2D FFTs:
double precision:: hrkx(ng),hrky(ng),rk(ng)
double precision:: xtrig(2*ng),ytrig(2*ng)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:ng),alk(ng)
integer:: kmag(ng,ng),kmax,kmaxred


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine init_spectral
! Initialises this module

implicit none

!Local variables:
double precision:: rkmax,rks,snorm
double precision:: anu,rkfsq,fsq
integer:: kx,ky,k

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

!-----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
rkmax=dble(ng/2)
kmax=nint(rkmax*sqrt(two))
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,ng
  do kx=1,ng
    k=nint(sqrt(rk(kx)**2+rk(ky)**2))
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(ng*ng) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(ng*ng)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=ng/2

!-----------------------------------------------------------------------
 !Define a variety of spectral operators:

 !Hyperviscosity coefficient (Dritschel, Gottwald & Oliver, JFM (2017)):
anu=cdamp*cof/rkmax**(2*nnu)
 !Assumes Burger number = 1.

 !Used for de-aliasing filter below:
rkfsq=(dble(ng)/3.d0)**2

 !Squared Coriolis frequency:
fsq=cof**2

do ky=1,ng
  do kx=1,ng
    rks=rk(kx)**2+rk(ky)**2
     !1-(H^2/3)*grad^2
    pope(kx,ky)=one+hbsq3*rks
     !Spectral c^2*grad^2 - f^2 operator:
    opak(kx,ky)=-(fsq+csq*rks)
     !Hyperviscous operator:
    diss(kx,ky)=anu*rks**nnu
     !De-aliasing filter:
    if (rks .gt. rkfsq) then
      filt(kx,ky)=zero
      rksq(kx,ky)=zero
      adop(kx,ky)=zero
      c2g2(kx,ky)=zero
      rlap(kx,ky)=zero
      helm(kx,ky)=zero
      rdis(kx,ky)=zero
      pdis(kx,ky)=zero
    else
      filt(kx,ky)=one
       !Squared wavenumber (-Laplace operator):
      rksq(kx,ky)=rks
       !Inverse of pope operator above (for gamma_t iteration):
      adop(kx,ky)=one/pope(kx,ky)
       !c^2*grad^2:
      c2g2(kx,ky)=-csq*rks
       !grad^{-2}:
      rlap(kx,ky)=-one/(rks+1.d-12)
       !(c^2*grad^2 - f^2)^{-1}:
      helm(kx,ky)=one/opak(kx,ky)
       !R operator in paper:
      rdis(kx,ky)=dt2i+diss(kx,ky)
       !P*R operator product:
      pdis(kx,ky)=pope(kx,ky)*rdis(kx,ky)
    endif
     !Semi-implicit operator for delta and gamma_l:
    simp(kx,ky)=one/(pope(kx,ky)*rdis(kx,ky)**2-opak(kx,ky))
     !Redfine damping operator for use in q_l evolution:
    diss(kx,ky)=two/(one+dt2*diss(kx,ky))
     !Multiply P = 1-(H^2/3)*grad^2 by 4/Delta{t} for use in time stepping:
    pope(kx,ky)=pope(kx,ky)*dt4i
  enddo
enddo

 !Ensure mean potentials and height anomaly remain zero:
rlap(1,1)=zero

return 
end subroutine

!======================================================================
subroutine main_invert(qs,ds,gs,hh,uu,vv,zz)
! Given the linearised PV qs, divergence ds and acceleration divergence gs
! (all in spectral space), this routine computes the dimensionless depth 
! anomaly hh and the velocity field (uu,vv) in physical space.  It also 
! returns the relative vorticity (zz) in physical space.

! Note: we assume zero momentum <(1+hh)(u,v)> = 0.  The mean flow is
! determined from this condition.

implicit none

 !Passed variables:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng) !Spectral
double precision:: hh(ng,ng),uu(ng,ng),vv(ng,ng) !Physical
double precision:: zz(ng,ng)                     !Physical

 !Local variables:
double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hh, uu & vv from ql, dd & gg.  The energy error is computed from 
 !      <(u-u0)^2+(v-v0)^2+c^2*(h-h0)^2>/<u^2+v^2+c^2*h^2>
 !      where <:> means a domain average and (u0,v0,h0) is the previous
 !      guess for (u,v,h).

 !Physical work arrays:
double precision:: bb(ng,ng),dd(ng,ng),gg(ng,ng)
double precision:: htot(ng,ng),hx(ng,ng),hy(ng,ng),wkp(ng,ng)
double precision:: badd(ng,ng)

 !Spectral work arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: uds(ng,ng),vds(ng,ng),fhs(ng,ng)

 !Other constants:
double precision:: uio,vio
double precision:: dhrms,durms,enorm

!-------------------------------------------------------
 !Define total dimensionless height:
htot=one+hh

 !Define spectral divergent velocity (never changes in iteration):
wkc=rlap*ds
 !This solves Lap(wkc) = dd in spectral space
call xderiv(ng,ng,hrkx,wkc,uds)
call yderiv(ng,ng,hrky,wkc,vds)

 !Obtain a physical space copy of ds & gs:
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
wka=gs
call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)

 !Obtain gg-2*dd^2 (filtered) for use below:
wkp=dd**2
call dealias(wkp)
badd=gg-two*wkp

 !Obtain fixed part of h inversion, f*q_l - gamma in spectral space:
fhs=cof*qs-gs

!-------------------------------------------------------
 !Iteratively solve for hh, uu & vv:

 !Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
   !Compute Green-Naghdi term B = (H^2/3)*(1+h)*(gamma-2*delta^2+2*J(u,v)):
  call jacob(uu,vv,wkp)
   !Filter J(u,v) (in wkp) before multiplying it by (1+h):
  call dealias(wkp)
   !Complete definition of B (= bb) then de-alias:
  bb=hbsq3*htot*(badd+two*wkp)
  call dealias(bb)

   !Compute terms involving B:
  wkp=hh
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
   !wka is hh in spectral space
  call xderiv(ng,ng,hrkx,wka,wkd)
  call spctop(ng,ng,wkd,hx,xfactors,yfactors,xtrig,ytrig)
   !hx is dh/dx in physical space
  call yderiv(ng,ng,hrky,wka,wkd)
  call spctop(ng,ng,wkd,hy,xfactors,yfactors,xtrig,ytrig)
   !hy is dh/dy in physical space

   !Compute div(hx*B,hy*B) (to be put into wka, in spectral space):
  hx=bb*hx
  hy=bb*hy
  call divs(hx,hy,wka)

   !Compute Laplace((1+h)*B) and add div(hx*B,hy*B) just computed:
  wkp=htot*bb
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
   !wkc is the spectral version of (1+h)*B

   !Invert [c^2*grad^2-f^2)]h = f*q_l - gamma
   !                          + Laplace((1+h)*B) + div(hx*B,hy*B):
  wka=helm*(fhs+wka-rksq*wkc) !helm includes the de-aliasing filter
   !fhs:  q_l + f*h (spectral)
   !rksq: -grad^2   (spectral)
   !helm: (c^2*grad^2-f^2)^{-1} (spectral)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
   !wkp: corrected de-aliased height field (to be hh below)

   !Compute rms error in hh:
  dhrms=sum((hh-wkp)**2)
   !Re-assign hh & htot:
  hh=wkp
  htot=one+hh

   !Store spectral copy of hh in wkb:
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)

   !Compute relative vorticity wkb = q_l + f*h (spectral):
  wkb=qs+cof*wkb

   !Solve Lap(wka) = wkb spectrally to define streamfunction (wka):
  wka=rlap*wkb

   !Bring back vorticity to physical space as zz:
  call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)

   !Compute streamfunction derivatives in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

   !New velocity components in spectral space, written in (wkb,wkd):
  wkb=uds-wkb  !uds is the fixed divergent part of uu
  wkd=vds+wkd  !vds is the fixed divergent part of vv

   !Convert to physical space as (hx,hy):
  call spctop(ng,ng,wkb,hx,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,hy,xfactors,yfactors,xtrig,ytrig)

   !Add mean flow (uio,vio):
  uio=-sum(hh*hx)*dsumi
  vio=-sum(hh*hy)*dsumi
  hx=hx+uio
  hy=hy+vio

   !Compute rms error in uu & vv:
  durms=sum((uu-hx)**2+(vv-hy)**2)

   !Re-assign velocity components:
  uu=hx
  vv=hy

   !Compute overall error:
  enorm=sqrt((durms+csq*dhrms)/sum(uu**2+vv**2+csq*hh**2))
enddo
 !Passing this, we have converged.

!------------------------------------------------------------------

return
end subroutine

!=================================================================
subroutine dealias(aa)
! De-aliases the array physical space array aa.  Returns aa in
! physical space.

implicit none

 !Passed array:
double precision:: aa(ng,ng)   !Physical

 !Work array:
double precision:: wka(ng,ng)  !Spectral

!---------------------------------------------------------
call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=================================================================
subroutine gradient(ff,ffx,ffy)
! Computes the gradient ffx = dF/dx & ffy = dF/dy of a field F.
! *** ff is in spectral space whereas (ffx,ffy) are in physical space

implicit none

 !Passed arrays:
double precision:: ff(ng,ng)             !Spectral
double precision:: ffx(ng,ng),ffy(ng,ng) !Physical

 !Local array:
double precision:: vtmp(ng,ng)           !Spectral

 !Get derivatives of F:
call xderiv(ng,ng,hrkx,ff,vtmp)
call spctop(ng,ng,vtmp,ffx,xfactors,yfactors,xtrig,ytrig)

call yderiv(ng,ng,hrky,ff,vtmp)
call spctop(ng,ng,vtmp,ffy,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=================================================================
subroutine jacob(aa,bb,cc)
! Computes the Jacobian of aa and bb and returns it in cc.
! All passed variables are in physical space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng),cc(ng,ng)           !Physical

 !Work arrays:
double precision:: ax(ng,ng),ay(ng,ng),bx(ng,ng),by(ng,ng) !Physical
double precision:: wka(ng,ng),wkb(ng,ng)                   !Spectral

!---------------------------------------------------------
cc=aa
call ptospc(ng,ng,cc,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Get derivatives of aa:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,ax,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,ay,xfactors,yfactors,xtrig,ytrig)

cc=bb
call ptospc(ng,ng,cc,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Get derivatives of bb:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,bx,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,by,xfactors,yfactors,xtrig,ytrig)

cc=ax*by-ay*bx

return
end subroutine

!=================================================================
subroutine divs(aa,bb,cs)
! Computes the divergence of (aa,bb) and returns it in cs.
! Both aa and bb in physical space but cs is in spectral space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng)   !Physical
double precision:: cs(ng,ng)             !Spectral

 !Work arrays:
double precision:: wkp(ng,ng)            !Physical
double precision:: wka(ng,ng),wkb(ng,ng) !Spectral

!---------------------------------------------------------
wkp=aa
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call xderiv(ng,ng,hrkx,wka,wkb)

wkp=bb
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,cs)

cs=wkb+cs

return
end subroutine

!===================================================================

subroutine spec1d(ss,spec)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.

implicit none

 !Passed variables:
double precision:: ss(ng,ng),spec(0:ng)

 !Local variables:
integer:: kx,ky,k

!--------------------------------------------------------
do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(1,1)
spec(k)=spec(k)+f14*ss(1,1)**2

 !y-independent mode:
do kx=2,ng
  k=kmag(kx,1)
  spec(k)=spec(k)+f12*ss(kx,1)**2
enddo

 !x-independent mode:
do ky=2,ng
  k=kmag(1,ky)
  spec(k)=spec(k)+f12*ss(1,ky)**2
enddo

 !All other modes:
do ky=2,ng
  do kx=2,ng
    k=kmag(kx,ky)
    spec(k)=spec(k)+ss(kx,ky)**2
  enddo
enddo

return
end subroutine

!===================================================================

end module     
