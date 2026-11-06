program random_pv
!-----------------------------------------------------------------
!    Generates a random phased baroclinic (mode 2) PV distribution
!    with spectrum S(k) = c * k^{2p-3} * exp[-(p-1)*(k/k_0)^2],
!    for integer p > 1.  Here, k_0 = <k^2 S>/<S> is the spectral
!    centroid.
!-----------------------------------------------------------------

use spectral

implicit none

 !Maximum x & y wavenumbers:
integer,parameter:: nw=ng/2

double precision:: qs(ng,ng),qa(ng,ng)
double precision:: qeddy,fac,s,phix,phiy
double precision:: p,p1,p2,rk0,amp,cx,sx,cy,sy

integer:: kx,ky,kxc,kyc,k,i,iy
integer, dimension(:), allocatable :: seed

!-----------------------------------------------------------------------
 !Initialise random number generator:
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=iseed
do i=1,iseed
  call random_seed(put=seed)
enddo

!-----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!-----------------------------------------------------------------------
write(*,*) ' The baroclinic (mode 2) PV has the spectrum'
write(*,*) '  S(k) = c k^{2p-3} * exp[-(p-1)*(k/k_0)^2].'
write(*,*)
write(*,*) ' Enter p > 1 & k_0:'
read(*,*) p,rk0

write(*,*) ' Enter the maximum baroclinic PV magnitude / f:'
read(*,*) qeddy
qeddy=cof*qeddy

! Generate spectral amplitude sqrt{k^{2(p-2)} * exp[-(p-1)*(k/k_0)^2]} :
fac=one/rk0**2
p1=p-one
p2=p-two
do ky=1,nw+1
  do kx=1,nw+1
    s=fac*(rk(kx)**2+rk(ky)**2)
    qs(kx,ky)=sqrt(fac*s**p2*exp(-p1*s))
  enddo
enddo

! Apply to generate full spectrum:
do ky=2,nw
  kyc=ng+2-ky
  do kx=2,nw
    kxc=ng+2-kx
    phix=twopi*rand(0)-pi
    phiy=twopi*rand(0)-pi
    cx=cos(phix)
    sx=sin(phix)
    cy=cos(phiy)
    sy=sin(phiy)
    amp=qs(kx,ky)
    qs(kx ,ky )=amp*cx*cy
    qs(kxc,ky )=amp*sx*cy
    qs(kx, kyc)=amp*cx*sy
    qs(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nw
  kxc=ng+2-kx
  phix=twopi*rand(0)-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nw
  kyc=ng+2-ky
  phiy=twopi*rand(0)-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

ky=nw+1
do kx=2,nw
  kxc=ng+2-kx
  phix=twopi*rand(0)-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=nw+1
do ky=2,nw
  kyc=ng+2-ky
  phiy=twopi*rand(0)-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

qs(1,1)=zero
qs(nw+1,nw+1)=zero

! Transform to physical space:
call spctop(ng,ng,qs,qa,xfactors,yfactors,xtrig,ytrig)

! Work out max/min values and scale:
fac=qeddy/max(maxval(qa),abs(minval(qa)))
qa=fac*qa

!-------------------------------------------------------------
! Write PV anomaly in each layer:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,-mrat2*qa
write(11,rec=2) zero, mrat1*qa
close(11)
! Here, mrat1 = m_1/(m_1+m_2) and mrat2 = m_2/(m_1+m_2) where
! m_j is the mass of layer j --- see constants.f90.

! Write zero fields of divergence and acceleration divergence:
qa=zero
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qa
write(11,rec=2) zero,qa
close(11)

open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qa
write(11,rec=2) zero,qa
close(11)

end program random_pv
