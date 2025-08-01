module constants

 !Include all modifiable parameters for use below:
use parameters

 !Integer kind parameters used in contours.f90 and 
 !congen.f90: dbleint represents up to 10^18
 !            halfint represents up to 127
integer,parameter:: dbleint=selected_int_kind(16)
integer,parameter:: halfint=selected_int_kind(-1)

integer,parameter:: ngm1=ng-1, ngp1=ng+1, ngp2=ng+2
integer,parameter:: ntm2=nt-2, ntm1=nt-1, ntp1=nt+1, ntp2=nt+2

integer,parameter:: ngridp=ng*nt,nbytes=4*(ngridp+1)

integer,parameter:: kda1=ngp1-ng/3, kda2=ngp1+ng/3
!      Wavenumbers between kda1 and kda2 are removed for de-aliasing;
!      this is done in longitude and in latitude (using great circles)

integer,parameter:: mgf=4, ngf=ng*mgf, ntf=nt*mgf
!      mgf:      Fine conversion grid ratio
!      ntf,ngf:  number of grid boxes in the x & y directions (fine grid)
integer,parameter:: mgu=16,ntu=mgu*nt, ngu=mgu*ng
!      mgu:      Ultra-fine conversion grid ratio (used in module congen.f90)
!      ntu,ngu:  number of grid boxes in the x & y directions (ultra-fine grid)

 !Maximum number of contour node points:
integer,parameter:: npm=80*ngridp
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npm/4

double precision,parameter:: zero=0.d0,  one=1.d0, two=2.d0,  three=3.d0
double precision,parameter:: four=4.d0, five=5.d0, six=6.d0, twelve=12.d0
double precision,parameter:: f12=one/two,  f13=one/three,  f23=two/three
double precision,parameter:: f14=one/four, f32=three/two,  f43=four/three
double precision,parameter:: f53=five/three, f56=five/six, f74=7.d0/four
double precision,parameter:: f16=one/six, f112=one/twelve, f1112=11.d0/twelve
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision,parameter:: hpi=pi/two, twopi=two*pi
double precision,parameter:: small=1.d-12, small3=small*small*small

double precision,parameter:: hlx=pi, hlxi=one/(pi+small), oms=one-small
! The x range is  -pi  <= x <=  pi   (longitude)
! The y range is -pi/2 <= y <= pi/2  (latitude)
double precision,parameter:: glxu=twopi/dble(ntu), glyu=pi/dble(ngu)

 !Basic constants:
double precision,parameter:: dl=twopi/dble(nt), dli=dble(nt)/(twopi+small)
double precision,parameter:: hpidl=(pi+dl)/two
double precision,parameter:: fpole=two*omega, csq=cgw**2, csqi=one/csq
double precision,parameter:: tfin=tsim*dble(nperiod)
double precision,parameter:: kappai=one/kappa, aa0i=one/aa0

 !Constants used in contour node placement (Fontane & Dritschel JCP 2009):
double precision,parameter:: amu=0.2d0, ell=6.25d0*dl
double precision,parameter:: elf=one/ell**2, densf=one/(amu*sqrt(ell))
double precision,parameter:: dm=amu**2*ell/four, dm2=dm**2, dmsq=four*dm2
double precision,parameter:: dmi=two/dm, d4small=small*dl**4

 !Logicals for stochastic and thermal forcing:
logical,parameter:: stoch=(esr .gt. zero), thermal=(alpha .gt. zero)

end module
