program ranpv
!-----------------------------------------------------------------
!    Generates a random phased PV distribution with an energy
!    spectrum E(k) = c k^{2p-3} * exp[-(p-1)*(k/k_0)^2], p > 1.
!    Here, k_0 is the enstrophy-energy centroid, the ratio
!    sqrt(enstrophy/energy), for all (integer) p > 1.
!-----------------------------------------------------------------

use spectral

implicit double precision(a-h,o-z)
double precision:: qs(ng,ng),qa(ng,ng)

 !Initialise inversion constants and arrays:
call init_spectral

nw=ng/2

write(*,*) ' We assume initial energy spectrum of the form'
write(*,*) '   E(k) = c k^{2p-3} * exp[-(p-1)*(k/k_0)^2]'
write(*,*)
write(*,*) ' Enter p > 1 & k_0:'
read(*,*) pow,ak0
! The potential enstrophy spectrum = (k_d^2 + k^2)*E(k)
! where k_d is the Rossby deformation wavenumber.

! Set the max abs QG PV anomaly = 4*pi:'
qeddy=4.d0*pi

write(*,*) ' Enter an integer seed for the random # generator:'
read(*,*) ngen
do i=1,ngen
  uni=rand(0)
enddo

! Generate potential enstrophy spectrum / k (actually, its square root)
! sqrt{k^{2(p-1)} * exp[-(p-1)*(k/k_0)^2]} :
efac=one/ak0**2
p1=pow-one
p2=pow-two
do ky=1,nw+1
  do kx=1,nw+1
    rks=rk(kx)**2+rk(ky)**2
    s=efac*rks
    qs(kx,ky)=sqrt(efac*(rks+kdsq)*s**p2*exp(-p1*s))
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

! Work out max/min values and total pot. enstrophy:
ens=zero
qamin=qa(1,1)
qamax=qamin
do ix=1,ng
  do iy=1,ng
    ens=ens+qa(iy,ix)**2
    qamin=min(qamin,qa(iy,ix))
    qamax=max(qamax,qa(iy,ix))
  enddo
enddo
ens=ens/(two*dble(ng*ng))

! Renormalise PV:
fmult=qeddy/max(abs(qamax),abs(qamin))
do ix=1,ng
  do iy=1,ng
    qa(iy,ix)=fmult*qa(iy,ix)
  enddo
enddo

! Work out max/min values and total pot. enstrophy:
ens=ens*fmult
qamin=qamin*fmult
qamax=qamax*fmult

! Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qa
close(11)

write(*,'(a,f12.5)') ' rms PV = ',sqrt(2.d0*ens)
write(*,'(a,f12.7,a,f11.7)') ' min PV = ',qamin,'  &  max PV = ',qamax

end program
