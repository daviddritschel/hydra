module evolution

! Module containing subroutines to evolve PV contours and fields.

use common

implicit none

double precision:: qt(ng,nt),qc(ng,nt),qd(ng,nt) !Various PVs

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine advect

! Main subroutine for advection

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops
!      to rebuild the contours in a separate memory space.
integer:: ireg,itime,jtime
logical:: ggen

!-----------------------------------------------------------------------
 !Re-initialise the PV fields following recontouring:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Save data periodically:
  itime=nint(t/dt)
  jtime=itime/ngsave
  if (ngsave*jtime .eq. itime) then
    call inversion(1) !1 indicates a contour->grid conversion must be done
    call savegrid(jtime+1)
    ggen=.false.
  else
    ggen=.true.
  endif
   !ggen is used to indicate if calling inversion is needed in advance below
  jtime=itime/ncsave
  if (ncsave*jtime .eq. itime) call savecont(jtime+1)

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
       !Prepare PV residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

    twist=twist-twistmax
  endif

   !Advect PV from time t to t + dt:
  call advance(ggen)
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Possibly save final data:
itime=nint(t/dt)
jtime=itime/ngsave
if (ngsave*jtime .eq. itime) then
  call inversion(1)
  call savegrid(jtime+1)
endif
jtime=itime/ncsave
if (ncsave*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration (module congen).

implicit none
integer:: i

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,n,npt

!------------------------------------------------------------
 !Re-instate qs as the PV anomaly (here qs is on the grid):
do i=1,nt
  qs(:,i)=qs(:,i)-cof
enddo

 !Convert qs to semi-spectral space:
call forfft(ng,nt,qs,trig,factors) 

 !Convert PV contours (x,y,z) to gridded PV anomaly as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(ng,nt,qc,trig,factors) 

 !Define (semi-spectral) residual PV qd = (1-F)*(qs-qc):
qd=qs-qc
call hipass(qd)

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current PV anomaly field is stored in qs (semi-spectral), and
! the full and residual PV fields needed in congen.f90 is stored in qq
! and qr (physical).

implicit none

! Local variables:
double precision:: qa(ng,nt),aqa
integer:: i

!------------------------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !          Reset  qs = F*(qs-qc) + qc + qd
 !            and  qd = (1-F)*(qs-qc)
 !Here F & (1-F) are low & high pass filter functions (see spectral.f90)
call reset

 !Obtain full PV anomaly (qs) on the grid:
call revfft(ng,nt,qs,trig,factors)

 !Add f to define full PV:
do i=1,nt
  qs(:,i)=qs(:,i)+cof
enddo

 !Correct average PV by enforcing zero average vorticity:
qa=qs*(one+hh)
aqa=average(qa)
qs=qs-aqa

 !Obtain PV anomaly due to contours (qc) on the grid:
call revfft(ng,nt,qc,trig,factors)

 !Define residual PV (qr) needed for recontouring:
do i=1,nt
  qr(:,i)=qs(:,i)-qc(:,i)-cof
enddo

 !Note: We are leaving this module after this; qd will be redefined
 !      upon re-entry in subroutine init.

return
end subroutine prepare

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for the semi-spectral fields qs & qd).

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer,parameter:: niter=2

 !Semi-spectral PV fields needed in time stepping:
double precision:: qsi(ng,nt),sqs(ng,nt)
double precision:: qdi(ng,nt),qdm(ng,nt),sqd(ng,nt)
 !Contour positions needed in time stepping:
double precision:: xi(npt),yi(npt),zi(npt)
 !Contour velocities:
double precision:: u(npt),v(npt),w(npt)
 !Other local quantities:
double precision:: xtmp,ytmp,ztmp,fac
integer:: i,iter

!------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call inversion(0)
 !If ggen is false, inversion was called previously at this time level.
 !A zero argument above means qc does not need to be regenerated.

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset  qs = F*(qs-qc) + qc + qd
 !            and  qd = (1-F)*(qs-qc)
 !Here F & (1-F) are low & high pass filter functions (see spectral.f90)
call reset

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all PV contours and fields:

 !Contours:
call velint(uu,vv,u,v,w)
 !Here (u,v,w) stands for the Cartesian velocity on the contours
 !at time level n, i.e. u(x^n,t^n)
do i=1,npt
   !Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xtmp=x(i)+dt2*u(i)
  ytmp=y(i)+dt2*v(i)
  ztmp=z(i)+dt2*w(i)
   !Ensure point remains on unit sphere:
  fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
  xi(i)=fac*xtmp
  yi(i)=fac*ytmp
  zi(i)=fac*ztmp
   !Preliminary guess for x^{n+1}:
  xtmp=x(i)+dt*u(i)
  ytmp=y(i)+dt*v(i)
  ztmp=z(i)+dt*w(i)
   !Ensure point remains on unit sphere:
  fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
  x(i)=fac*xtmp
  y(i)=fac*ytmp
  z(i)=fac*ztmp
enddo

 !Calculate the source terms (sqs,sqd) for semi-spectral PV (qs,qd):
call pvsource(sqs,sqd)

 !Update PV fields:
qdi=qd
qdm=qd+dt4*sqd
qd=diss*(qdm+dt4*sqd)-qdi
qsi=qs+dt2*sqs
qs=qs+dt*sqs

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion(1)

   !Interpolate gridded velocity (uu,vv) at contour nodes as (u,v):
  call velint(uu,vv,u,v,w)
  do i=1,npt
    xtmp=xi(i)+dt2*u(i)
    ytmp=yi(i)+dt2*v(i)
    ztmp=zi(i)+dt2*w(i)
    fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
    x(i)=fac*xtmp
    y(i)=fac*ytmp
    z(i)=fac*ztmp
  enddo
   !Now (x,y,z) contain a new guess for x^{n+1}.

   !Calculate the source terms (sqs,sqd) for PV (qs,qd):
  call pvsource(sqs,sqd)

   !Update PV fields:
  qd=diss*(qdm+dt4*sqd)-qdi
  qs=qsi+dt2*sqs
enddo

 !Apply latitudinal hyperviscous damping to qd:
call latdamp(qd)
 !The longitudinal part is incorporated above in diss.

 !Advance time:
t=t+dt

return
end subroutine advance

!=======================================================================

subroutine pvsource(sqs,sqd)

! Gets the nonlinear source terms (sqs,sqd) for the PV (qs,qd) in
! semi-spectral space.

implicit none

 !Passed variables:
double precision:: sqs(ng,nt),sqd(ng,nt)

 !Local variables:
double precision:: wkp(ng,nt),wkq(ng,nt)
integer:: i

!---------------------------------------------------------------
 !qd source --- only NL advection term is needed:
call deriv(ng,nt,rk,qd,wkp)
call revfft(ng,nt,wkp,trig,factors) ! wkp = d(qd)/d(lon)
sqd=qd
call revfft(ng,nt,sqd,trig,factors) 
call latder(sqd,wkq)                ! wkq = d(qd)/d(lat)
do i=1,nt
  sqd(:,i)=-uu(:,i)*clati*wkp(:,i)-vv(:,i)*wkq(:,i)
enddo
 !Convert to semi-spectral space and de-alias:
call dealias(sqd)

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:
call deriv(ng,nt,rk,qs,wkp)
call revfft(ng,nt,wkp,trig,factors) ! wkp = d(qs)/d(lon)
sqs=qs
call revfft(ng,nt,sqs,trig,factors) 
call latder(sqs,wkq)                ! wkq = d(qs)/d(lat)
do i=1,nt
  sqs(:,i)=-uu(:,i)*clati*wkp(:,i)-vv(:,i)*(wkq(:,i)+bet) !include beta*v
enddo
 !Convert to semi-spectral space and de-alias:
call dealias(sqs)

return
end subroutine pvsource

!=======================================================================

subroutine dgsource(sds,sgs)

! Gets the source terms (sds,sgs) for delta and gamma (ds,gs) in
! semi-spectral space.

! Note that (sds,sgs) only include the nonlinear terms; the linear
! terms are used to carry out the balanced flow inversion.
  
implicit none

 ! Passed variables:
double precision:: sds(ng,nt),sgs(ng,nt)

 ! Local variables:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: wke(ng,nt),wkf(ng,nt),wkp(ng,nt),wkq(ng,nt)
double precision:: wku(ng,nt),wkv(ng,nt)
double precision:: dd(ng,nt),rhs(ng),avgval
integer:: i

!-----------------------------------------------------------------
 !Compute A = z*U + dV/d(lon) -> wka & B = z*V - dU/d(lon) -> wkb,
 !where U = u/r & V = v/r, while z = sin(lat) and r = cos(lat):
do i=1,nt
  wku(:,i)=clati*uu(:,i)
  wkv(:,i)=clati*vv(:,i)
  wka(:,i)=wkv(:,i)
  wkb(:,i)=wku(:,i)
enddo
 !*** Do not re-use wku ***

 !Convert wka & wkb to semi-spectral space:
call forfft(ng,nt,wka,trig,factors) 
call forfft(ng,nt,wkb,trig,factors) 

 !Compute longitudinal derivatives:
call deriv(ng,nt,rk,wka,wkc)
call deriv(ng,nt,rk,wkb,wkd)

 !Recover physical fields:
call revfft(ng,nt,wkc,trig,factors) 
call revfft(ng,nt,wkd,trig,factors) 

 !Complete definition of A & B and define wkf = f*zeta & wkb = 2*f*beta*v:
do i=1,nt
  wka(:,i)=slat*wku(:,i)+wkc(:,i) !slat = z = sin(lat)
  wkb(:,i)=slat*wkv(:,i)-wkd(:,i)
  wkf(:,i)=cof*zz(:,i)            !cof = f
  wkd(:,i)=dfb*vv(:,i)            !dfb = 2*f*beta
enddo
 !*** Do not re-use wka, wkb, wkf or wkd ***

 !Get physical space velocity divergence -> dd:
dd=ds
call revfft(ng,nt,dd,trig,factors) 
 !*** Do not re-use dd ***

 !Get physical space derivatives of divergence:
call deriv(ng,nt,rk,ds,wkv)
call revfft(ng,nt,wkv,trig,factors) ! wkv = d(delta)/d(lon)
call latder(dd,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the velocity divergence tendency (sds):
wkp=uu**2+vv**2 ! *** Do not re-use wkp ***
sds=two*(wka*(zz-wka)-wkb*(wkb+dd))-wkp-dd**2-wku*wkv-vv*wkc
 !wku = u/r on rhs above; *** wka & wkb now safe to re-use below ***

 !Convert to semi-spectral space and de-alias:
call dealias(sds)

 !Compute the nonlinear part of the acceleration divergence tendency (sgg):

 !First compute div((u,v)*Z) (immediately below wkf = Z = f*zeta):
wka=wkf
call forfft(ng,nt,wka,trig,factors)
call deriv(ng,nt,rk,wka,wkv)
call revfft(ng,nt,wkv,trig,factors) ! wkv = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u,v)*Z) + 2*f*beta*v in sgs temporarily:
sgs=dd*wkf+wku*wkv+vv*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs)
 !*** wkf can now be re-used ***

 !Define B = c^2*h + (u^2 + v^2)/2 -> wkb:
wkb=csq*hh+f12*wkp
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(ng,nt,rk,wkb,wkf)

 !Compute div((u,v)*h) -> wke:
wka=hh
call forfft(ng,nt,wka,trig,factors) 
call deriv(ng,nt,rk,wka,wkv)
call revfft(ng,nt,wkv,trig,factors) ! wkv = dh/d(lon)
call latder(hh,wkc)                 ! wkc = dh/d(lat)
wke=dd*hh+wku*wkv+vv*wkc
 !Note wku = u/r on rhs above

 !Compute Laplacian of div((u,v)*h) in wke after de-aliasing:
call dealias(wke)     ! Now wke is in semi-spectral space
call laplace(wke,wka) ! wka is returned in semi-spectral space

 !Complete calculation of nonlinear part of gamma tendency:
sgs=csq*wka+fpole*wkf-sgs
 !Note: each part has been separately de-aliased.

!-----------------------------------------------------
 !Remove global mean values of sds & sgs:
rhs=sds(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sds(:,1)=sds(:,1)-avgval

rhs=sgs(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sgs(:,1)=sgs(:,1)-avgval

end subroutine dgsource

!=======================================================================

subroutine balance

! Solves the equations resulting from setting delta_t = gamma_t = 0
! to determine delta and gamma; also determines h, u & v.  Effectively,
! this is PV inversion.

implicit none

! Error in energy norm below which solution is said to be converged:
double precision,parameter:: tole=1.d-9

! Local variables:
double precision:: hho(ng,nt),uuo(ng,nt),vvo(ng,nt)
double precision:: dso(ng,nt),gso(ng,nt)
double precision:: sds(ng,nt),sgs(ng,nt)
double precision:: berr

!-----------------------------------------------------------------------
! Initialise (use previous field values as first guesses):
dso=ds
gso=gs
hho=hh
uuo=uu
vvo=vv

! Iterate to find balance:
do
  ! Invert PV, delta and gamma to obtain the velocity field:
  call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
  ! Note: qt, ds & gs are in spectral space while 
  !       hh, uu, vv, qq and zz are in physical space.

  ! Compute ***relative*** energy norm error:
  berr=average((uu-uuo)**2+(vv-vvo)**2+csq*(hh-hho)**2)/ &
       average(uu**2+vv**2+csq*hh**2)

  ! Leave this loop if error is small enough:
  if (berr .lt. tole) exit
  
  ! Otherwise, copy new fields into old ones for the next iteration:
  uuo=uu
  vvo=vv
  hho=hh
  dso=ds
  gso=gs

  ! Compute NL source terms for delta and gamma:
  call dgsource(sds,sgs)

  ! Update divergence and acceleration divergence:
  call helminv(sgs,ds,gs)
  ds=-csqi*ds
  gs=-sds
enddo

return
end subroutine balance

!=======================================================================

subroutine inversion(iopt)

! Finds the gridded dimensionless height anomaly (hh) and velocity
! field (uu,vv) from the PV contours and the (semi-spectral)
! divergence (ds) and acceleration divergence (gs).  The latter
! are determined by the balance conditions delta_t = gamma_t = 0.  

! if iopt = 1, qc (semi-spectral) is found from the PV contours.
  
implicit none
integer:: iopt

!------------------------------------------------------------
if (iopt .eq. 1) then
   !Convert PV contours (x,y,z) to gridded PV anomaly as qc:
  call con2grid(qc)

   !Convert qc to semi-spectral space:
  call forfft(ng,nt,qc,trig,factors) 
endif

 !Combine fields to update qt with full (semi-spectral) field,
 !qt = F*(qs-qc) + qc + qd, where F is a low pass filter:
qt=qs-qc
call lopass(qt)
qt=qt+qc+qd

 !Find balanced divergence and acceleration divergence:
call balance
 !This also returns the updated dimensionless height anomaly hh,
 !the velocity field (uu,vv), the gridded PV anomaly (qq) and
 !relative vorticity (zz).

return
end subroutine inversion

!=======================================================================

subroutine reset
! Resets qs & qd at the beginning of each time step, or
! just before recontouring at the end of a run.

implicit none

 !Axillary PV array:
double precision:: qa(ng,nt)

!-------------------------------------------------------------------
 !Reset  qs = F*(qs-qc) + qc + qd  and  qd = (1-F)*(qs-qc)
 !where F & (1-F) are low & high pass filter functions
 !(see spectral.f90)

 !Convert PV contours (x,y,z) to gridded PV anomaly as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(ng,nt,qc,trig,factors) 
qa=qs-qc
call lopass(qa)
qs=qa+qc+qd
qd=qs-qc
call hipass(qd)

return
end subroutine

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wkp(ng,nt)
double precision:: zmax,frmax,hmax,hmin
double precision:: zrms,drms

!----------------------------------------------------------------------
 !Increment the integral of |zeta|_max:
zmax=maxval(abs(zz))
twist=twist+dt*zmax

 !Compute diagnostics:
frmax=sqrt(csqi*maxval((uu**2+vv**2)/(one+hh)))
call getrms(zz,zrms)
wkp=ds
call revfft(ng,nt,wkp,trig,factors)
call getrms(wkp,drms)
hmax=maxval(hh)
hmin=minval(hh)

 !Record |zeta|_max/f_pole, Fr_max, h_min, h_max, zeta_rms, delta_rms:
write(16,'(1x,f12.5,6(1x,f14.8))') t,zmax/fpole,frmax,hmin,hmax,zrms,drms

return
end subroutine diagnose

!=======================================================================

subroutine savegrid(igrids)

! Saves PV, energy and various spectra at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: wkp(ng,nt)
double precision:: spec(ng)
double precision:: ekin,epot,etot,angm
integer:: i,m

!-----------------------------------------------------------------
 !Compute energy components and total as well as angular momentum:
wkp=(one+hh)*(uu**2+vv**2)
ekin=twopi*average(wkp)
wkp=hh**2
epot=twopi*csq*average(wkp)
etot=ekin+epot
do i=1,nt
  wkp(:,i)=clat*((one+hh(:,i))*uu(:,i)+omega*clat*hh(:,i))
enddo
angm=fourpi*average(wkp)

 !Write energies & angular momentum to ecomp.asc:
write(15,'(f13.6,4(1x,f16.9))') t,ekin,epot,etot,angm
write(*,'(a,f13.6,a,f14.7,a,f14.7)') ' t = ',t,'     E_tot = ',etot, &
                                     '     A = ',angm

!-----------------------------------------------------------------
 !Compute 1d longitudinal power spectra for various fields:
wkp=zz
call forfft(ng,nt,wkp,trig,factors)
call longspec(wkp,spec)
spec=log10(spec+1.d-32)
write(51,'(f13.6,1x,i5)') t,ng
do m=1,ng
  write(51,'(2(1x,f12.8))') alm(m),spec(m)
enddo

call longspec(ds,spec)
spec=log10(spec+1.d-32)
write(52,'(f13.6,1x,i5)') t,ng
do m=1,ng
  write(52,'(2(1x,f12.8))') alm(m),spec(m)
enddo

call longspec(gs,spec)
spec=log10(spec+1.d-32)
write(53,'(f13.6,1x,i5)') t,ng
do m=1,ng
  write(53,'(2(1x,f12.8))') alm(m),spec(m)
enddo

wkp=hh
call forfft(ng,nt,wkp,trig,factors)
call longspec(wkp,spec)
spec=log10(spec+1.d-32)
write(54,'(f13.6,1x,i5)') t,ng
do m=1,ng
  write(54,'(2(1x,f12.8))') alm(m),spec(m)
enddo

!-----------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(31,rec=igrids) real(t),real(qq)
wkp=ds
call revfft(ng,nt,wkp,trig,factors)
write(32,rec=igrids) real(t),real(wkp)
wkp=gs
call revfft(ng,nt,wkp,trig,factors)
write(33,rec=igrids) real(t),real(wkp)
write(34,rec=igrids) real(t),real(hh)
write(35,rec=igrids) real(t),real(zz)

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: qa(ng,nt)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the contours subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') n,npt,t,dq

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
qa=qt-qc
call revfft(ng,nt,qa,trig,factors)
write(83,rec=irec) real(t),real(qa)

 !Save PV contours:
open(81,file='contours/qqindex'//pind,form='unformatted',status='replace')
write(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='replace')
write(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
