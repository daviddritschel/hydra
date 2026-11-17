!######################################################################
!  Finds balanced fields given only the total PV in qq_init.r8 or
!  in evolution/qq.r4, using first-order delta-gamma balance in which
!  delta_t and gamma_t are set to zero to diagnose all fields.
!  Here, delta is the velocity divergence while gamma is the
!  acceleration divergence.

!  This routine can be used for both initialisation and post-processing.

!  Revised 14 November 2026 by D G Dritschel @ Kyoto
!######################################################################

program dgbal

 !Import spectral module, constants and parameters:
use spectral

implicit none

 !Global commons (to this program):

 !Physical fields (PV, delta, gamma, thickness and relative vorticity):
double precision:: qq(ng,nt)
double precision:: dd(ng,nt)
double precision:: gg(ng,nt)
double precision:: hh(ng,nt)
double precision:: zz(ng,nt)

 !Physical velocity field:
double precision:: uu(ng,nt)
double precision:: vv(ng,nt)

 !Spectral fields (delta and gamma):
double precision:: ds(ng,nt)
double precision:: gs(ng,nt)

 !Bottom topography and associated spectral fields (if used):
double precision,allocatable,dimension(:,:):: bb,tsds,tsgs

 !Other quantities:
double precision:: ekin,epot,etot,angm,t
real:: qqr4(ng,nt),tr4
integer:: iopt,loop,iread,i

 !---------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !If fixed topography is present, allocate memory and read in its form:
if (topogr) then
  allocate(bb(ng,nt),tsds(ng,nt),tsgs(ng,nt))   
  open(11,file='bb_init.r8',form='unformatted',&
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,bb
  close(11)
  !Store its contribution to the gamma and delta source terms:
  tsgs=csq*bb
  call forfft(ng,nt,tsgs,trig,factors)
  call laplace(tsgs,tsds)
  call revfft(ng,nt,tsds,trig,factors)
endif

write(*,*) ' Choose one of the following options:'
write(*,*)
write(*,*) ' (0) initialise from data in qq_init.r8'
write(*,*) ' (1) post-process from data in evolution/qq.r4'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt
write(*,*)

if (iopt .eq. 0) then
   !Only find balance at t = 0 to initialise; read total PV:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,qq
  close(11)
  
   !Find balanced fields:
  call balance

   !Save initial (total) PV:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,qq
  close(11)
   !Note, the PV is modified by an overall constant offset due to the
   !mathematical requirement that the domain integral of the relative
   !vorticity be zero: see main_invert in spectral.f90.

   !Save initial velocity divergence:
  open(11,file='dd_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,dd
  close(11)

   !Save initial acceleration divergence:
  open(11,file='gg_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,gg
  close(11)

   !Save initial dimensionless height anomaly:
  open(11,file='hh_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,hh
  close(11)

   !Save initial relative vorticity:
  open(11,file='zz_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,zz
  close(11)

  write(*,*)
  write(*,*) ' The balanced fields are available in xx_init.r8 where'
  write(*,*) ' xx = qq, dd, gg, hh or zz.'
  write(*,*)
  
else
   !Process all times using only the (total) PV in evolution/qq.r4:
  open(30,file='evolution/qq.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)

  open(31,file='evolution/bqq.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(32,file='evolution/bdd.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(33,file='evolution/bhh.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(34,file='evolution/bgg.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(35,file='evolution/bzz.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(51,file='evolution/becomp.asc',status='replace')

   !Read data at all available times and balance:
  loop=0
  do
    loop=loop+1
    iread=0
     !Read the PV anomaly from qq.r4:
    read(30,rec=loop,iostat=iread) tr4,qqr4
    if (iread .ne. 0) exit
    qq=dble(qqr4)
    t=dble(tr4)

     !Balance using only the PV anomaly at time t = t:
    write(*,'(a,f9.2)') ' *** Processing t = ',t
    call balance

     !Write data:
    write(31,rec=loop) tr4,real(qq)
    write(32,rec=loop) tr4,real(dd)
    write(33,rec=loop) tr4,real(hh)
    write(34,rec=loop) tr4,real(gg)
    write(35,rec=loop) tr4,real(zz)

     !Write balanced energy components and angular momentum:
    write(51,'(f13.6,4(1x,f16.9))') t,ekin,epot,etot,angm
  enddo

  close(30)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(51)

  write(*,*)
  write(*,*) ' The balanced fields are available in evolution/bxx.r4'
  write(*,*) ' where xx = qq, dd, hh, gg & zz.'
  write(*,*)
  write(*,*) ' Also, becomp.asc contains the energy components.'
  write(*,*)

endif


 !Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine balance
 !Given the *total* PV in qq, this routine finds all the other
 !fields by requiring delta_t = gamma_t = 0.
  
implicit none

 !Local parameters and variables:

double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      d, g, h, u & v from q (for j = 1 & 2). The energy norm is
 !      <du^2+dv^2+c^2*dh^2> where <:> means a horizontal domain average
 !      and (du,dv,dh) is the difference from the previous guess.

 !rlx:  relaxation parameter to help convergence (rlx = 1 does not work):
double precision,parameter:: rlx=0.5d0, rlxc=one-rlx

 !Divergent part of the velocity field:
double precision:: ud(ng,nt)
double precision:: vd(ng,nt)

 !Previous guess needed to compute the energy norm:
double precision:: hhpre(ng,nt)
double precision:: uupre(ng,nt)
double precision:: vvpre(ng,nt)

 !Used to store source for the height inversion:
double precision:: sh(ng,nt)

 !Spectral nonlinear delta dn gamma tendencies
double precision:: sds(ng,nt)
double precision:: sgs(ng,nt)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)

 !Other quantities:
double precision:: enorm
integer:: i,m

 !----------------------------------------------------------------------
 !Start with zero divergence and acceleration divergence:
ds=zero
gs=zero

 !Initialise previous guesses for h, u & v:
hh=zero
uu=zero
vv=zero

hhpre=zero
uupre=zero
vvpre=zero

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Iteratively solve for h, u & v (initial guess for d = g = 0):

 !Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
   !Correct average PV by enforcing zero average vorticity:
  wkc=one+hh
  wka=qq*wkc
  qq=qq-average(wka)
   !Recall q = (z + f)/(1 + h)

   !---------------------------------------------------------------
   !Define initial "source" sh needed to invert hh:
  do i=1,nt
    sh(:,i)=bet*uu(:,i)-cof*(qq(:,i)-cof)*wkc(:,i)
  enddo
   !These use the fact that ds = gs = 0 in this first iterative loop.
   !Convert source to semi-spectral space and de-alias:
  call dealias(sh)

   !Solve for hh:
  sh=-csqi*sh
  call helminv(sh,hh,wkb)
   !Solves f^2*h - c^2*Lap(h) = S_h

   !Return hh to physical space:
  call revfft(ng,nt,hh,trig,factors) 

   !---------------------------------------------------------------
   !Correct average PV by enforcing zero average vorticity:
  wkc=one+hh
  wka=qq*wkc
  qq=qq-average(wka)

   !Compute relative vorticity (now guaranteed to have zero average):
  do i=1,nt
    zz(:,i)=qq(:,i)*wkc(:,i)-cof
  enddo

   !Compute non-divergent velocity field (as here ds = 0):
  
   !Convert the relative vorticity to semi-spectral space as wkd and de-alias:
  wkd=zz
  call dealias(wkd)

   !Invert Laplace's operator on relative vorticity (wkd):
  call laplinv(wkd,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(ng,nt,rk,wka,wkc)

   !Define velocity in semi-spectral space:
  do m=1,nt
    uu(:,m)=-wkb(:,m)
    vv(:,m)=clati*wkc(:,m)
  enddo

   !Return to physical space:
  call revfft(ng,nt,uu,trig,factors)
  call revfft(ng,nt,vv,trig,factors)

   !Compute relative error (measured by the energy norm):
  wka=(uu-uupre)**2+(vv-vvpre)**2+csq*(hh-hhpre)**2
  wkb=uu**2+vv**2+csq*hh**2
  enorm=average(wka)/(average(wkb)+small)
   !Small = 1.d-12 is added to avoid dividing by zero for a rest state

  if (enorm > one) then
     write(*,*) ' *** Diverging ***'
     stop
  endif

   !Apply relaxation:
  hh=rlx*hh+rlxc*hhpre
  uu=rlx*uu+rlxc*uupre
  vv=rlx*vv+rlxc*vvpre
  
   !Store current guess to recompute error in the next iteration:
  hhpre=hh
  uupre=uu
  vvpre=vv

enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Passing this, we have converged.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Iterate to find all the balanced fields, including delta & gamma:
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 !Energy norm error (must be > tole to start):
enorm=f12
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do while (enorm .gt. tole)

   !Compute nonlinear source terms S_d & S_g for delta & gamma:
  call source(sds,sgs)

   !Solve for delta:
  sgs=-csqi*sgs 
  call helminv(sgs,ds,wkb) !wkb is not used 
  !Solves f^2*d - c^2*Lap(d) = S_g,

   !gamma is just the negative of the divergence source term S_d:
  gs=-sds

   !-----------------------------------------------------
   !Compute divergent part of the velocity field (ud,vd):

   !Invert Laplace's operator on the divergence (ds):
  call laplinv(ds,wka,vd)
   !Here the divergence potential is wka, and its latitudinal
   !derivative is vd, the divergent meridional velocity
   !(all in semi-spectral space)

   !Compute divergent zonal velocity ud:
  call deriv(ng,nt,rk,wka,ud)
  do m=1,nt
    ud(:,m)=clati*ud(:,m)
  enddo
   !clati = 1/r, where r = cos(lat)

   !--------------------------------------------------------
   !Compute non-divergent velocity field and add to (ud,vd):

   !First correct average PV by enforcing zero average vorticity:
  wkc=one+hh
  wka=qq*wkc
  qq=qq-average(wka)

   !Compute relative vorticity (now guaranteed to have zero average):
  do i=1,nt
    zz(:,i)=qq(:,i)*wkc(:,i)-cof
  enddo

   !Convert the relative vorticity to semi-spectral space as wkd and de-alias:
  wkd=zz
  call dealias(wkd)
   !Preserve wkd until it is used to define sh below

   !Invert Laplace's operator on relative vorticity:
  call laplinv(wkd,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(ng,nt,rk,wka,wkc)

   !Define hh in spectral space as sh for use in the loop below:
  sh=hh
  call forfft(ng,nt,sh,trig,factors)

  do m=1,nt
    uu(:,m)=ud(:,m)-wkb(:,m)
    vv(:,m)=vd(:,m)+clati*wkc(:,m)
    sh(:,m)=gs(:,m)+bet*uu(:,m)-cof*(wkd(:,m)-cof*sh(:,m))
  enddo
   !All the terms in this loop are semi-spectral.

  !-------------------------------------------------------------------
   !Solve for hh:
  sh=-csqi*sh
  call helminv(sh,hh,wkb)
   !Solves f^2*h - c^2*Lap(h) = S_h

   !Return h, u & v to physical space:
  call revfft(ng,nt,hh,trig,factors) 
  call revfft(ng,nt,uu,trig,factors) 
  call revfft(ng,nt,vv,trig,factors) 

   !-----------------------------------------------------
   !Compute relative error (measured by the energy norm):
  wka=(uu-uupre)**2+(vv-vvpre)**2+csq*(hh-hhpre)**2
  wkb=uu**2+vv**2+csq*hh**2
  enorm=average(wka)/(average(wkb)+small)
   !Small = 1.d-12 is added to avoid dividing by zero for a rest state

  if (enorm > one) then
     write(*,*) ' *** Diverging ***'
     stop
  endif

   !Apply relaxation:
  hh=rlx*hh+rlxc*hhpre
  uu=rlx*uu+rlxc*uupre
  vv=rlx*vv+rlxc*vvpre

   !Store current guess to recompute error in the next iteration:
  hhpre=hh
  uupre=uu
  vvpre=vv

enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Passing this, we have converged.

 !Return delta and gamma to physical space:
dd=ds
call revfft(ng,nt,dd,trig,factors)

gg=gs
call revfft(ng,nt,gg,trig,factors)

 !Compute energy components, divided by rho*H where H = H_1+H_2:

 !Kinetic energy:
wkc=(one+hh)*(uu**2+vv**2)
ekin=twopi*average(wkc)
 !This is the area integral of the kinetic energy density

 !Potential energy:
if (topogr) then
   wkc=hh*(two*wkc+hh)
else
   wkc=hh**2
endif
epot=twopi*csq*average(wkc)

 !Total energy (sum of the above):
etot=ekin+epot

 !Compute angular momentum:
do i=1,nt
  wka(:,i)=clat*((one+hh(:,i))*uu(:,i)+omega*clat*hh(:,i))
enddo
 !wka: angular momentum density
 !Omega = planetary rotation rate (= 2*pi, see constants.f90)
angm=fourpi*average(wka)

return
end subroutine balance

!=======================================================================

subroutine source(sds,sgs)

! Gets the nonlinear parts of the source terms (sds,sgs) for divergence
! and acceleration divergence (ds,gs) --- all in semi-spectral space.

! The semi-spectral fields ds & gs are spectrally truncated.
! Note, h, u, v & z obtained by main_invert before calling this routine
! are also spectrally truncated. These variables are globally available.

implicit none

 !Passed variables:
double precision:: sds(ng,nt),sgs(ng,nt)

 !Local variables:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: wke(ng,nt),wkf(ng,nt),wkg(ng,nt),wkh(ng,nt)
double precision:: om(ng,nt),usq(ng,nt),dd(ng,nt)
double precision:: rhs(ng),avgval
integer:: i

!--------------------------------------------------------------
 !First define Omega = u/r (where r = cos(lat)):
do i=1,nt
  om(:,i)=clati*uu(:,i)
enddo
 !*** Do not re-use om below ***

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain delta source (sds):

 !Compute A = z*u/r + (1/r)*dv/d(lon) -> wka
 !    and B = z*v/r - (1/r)*du/d(lon) -> wkb:
 !where z = sin(lat) and r = cos(lat):

 !Convert vv (= v) & om (= u/r) to semi-spectral space as wka & wkb:
wka=vv
call forfft(ng,nt,wka,trig,factors) 
wkb=om
call forfft(ng,nt,wkb,trig,factors) 

 !Compute longitudinal derivatives of wka & wkb:
call deriv(ng,nt,rk,wka,wkc)
call deriv(ng,nt,rk,wkb,wkd)

 !Recover physical fields (inverse FFT):
call revfft(ng,nt,wkc,trig,factors) 
call revfft(ng,nt,wkd,trig,factors) 

 !Complete definition of A & B:
do i=1,nt
  wka(:,i)=slat*om(:,i)+clati*wkc(:,i)
  wkb(:,i)=tlat*vv(:,i)-wkd(:,i)
enddo
 !tlat = tan(lat) above; *** Do not re-use wka or wkb ***

 !Get physical space velocity divergence delta -> dd:
dd=ds
call revfft(ng,nt,dd,trig,factors) 
 !*** Do not re-use dd ***

 !Get physical space derivatives of delta:
call deriv(ng,nt,rk,ds,wkg)
call revfft(ng,nt,wkg,trig,factors) ! wkg = d(delta)/d(lon)
call latder(dd,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the delta tendency:
usq=uu**2+vv**2 ! u^2 + v^2;  *** Do not re-use usq ***
sds=two*(wka*(zz-wka)-wkb*(wkb+dd))-usq-dd**2-om*wkg-vv*wkc
 !om = u/r on rhs above; *** wka & wkb now safe to re-use below ***

 !Add effect of fixed topography (if present):
if (topogr) sds=sds-tsds

 !Convert to semi-spectral space and de-alias:
call dealias(sds)

!---------------------------------------------------------------
 !gamma source:

 !First obtain c^2*Lap(div(h*(u,v))) -> wkg:
wka=hh
call forfft(ng,nt,wka,trig,factors) 
call deriv(ng,nt,rk,wka,wkd)
call revfft(ng,nt,wkd,trig,factors) ! wkd = dh/d(lon)
call latder(hh,wkc)                 ! wkc = dh/d(lat)
wke=dd*hh+om*wkd+vv*wkc
 !Note om = u/r on rhs above

 !If thermal damping is present, add its contribution:
if (thermal) wke=wke+rth*hh

 !Compute Laplacian of wke = div(h*(u,v)) after de-aliasing:
call dealias(wke) !Also converts wke to semi-spectal space
call laplace(wke,wkg)
wkg=csq*wkg
 !wkg = c^2*Lap(div(h*(u,v))) in semi-spectral space

 !Define wkf = f*zeta & wkd = 2*f*beta*v:
do i=1,nt
  wkf(:,i)=cof*zz(:,i)            !cof = f
  wkd(:,i)=dfb*vv(:,i)            !dfb = 2*f*beta
enddo

 !Compute div((u,v)*Z) where Z = f*zeta:
wka=wkf
call forfft(ng,nt,wka,trig,factors)
call deriv(ng,nt,rk,wka,wke)
call revfft(ng,nt,wke,trig,factors) ! wke = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u,v)*Z) + 2*f*beta*v in sgs temporarily:
sgs=dd*wkf+om*wke+vv*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs)

 !Define B = c^2*h + (u^2 + v^2)/2 -> wkb:
wkb=csq*hh+f12*usq
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Add effect of fixed topography (if present):
if (topogr) wkb=wkb+tsgs
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(ng,nt,rk,wkb,wkf)

 !Complete calculation of (de-aliased) nonlinear part of gamma tendency:
sgs=wkg+fpole*wkf-sgs
 !wkg is defined above as the Laplacian of c^2*div(h*(u,v))

!-----------------------------------------------------
 !Remove global mean values of sds & sgs:
rhs=sds(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sds(:,1)=sds(:,1)-avgval

rhs=sgs(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sgs(:,1)=sgs(:,1)-avgval

return
end subroutine source

 ! End main program
end program dgbal
!=======================================================================
