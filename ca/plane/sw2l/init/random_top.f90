program random_top
! |---------------------------------------------------|
! |      This routine sets up random topography.      |
! |      Optionally, one can add A*sin(y).            |
! |---------------------------------------------------|

use spectral

implicit none

 !Maximum x & y wavenumbers:
integer,parameter:: nw=ng/2

 !Common arrays, constants:
double precision:: ss(ng,ng),bb(ng,ng)
double precision:: fac,s,phix,phiy,rms,p,p1,p2,rk0
double precision:: amp,cx,sx,cy,sy,tamp

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

!----------------------------------------------------------------------
 !Initialise FFTs:
call init_spectral

!-----------------------------------------------------
 !Compute random field with a given variance spectrum:
write(*,*) ' We consider a topography spectrum of the form'
write(*,*) ' c k^{2p+1} exp[-p(k/k_0)^2]. Enter p > 0 and k_0:'
read(*,*) p,rk0
write(*,*) ' Enter the rms value of the field (determines c):'
read(*,*) rms
write(*,*) ' This field is superposed on top of sinusoidal'
write(*,*) ' topography in y with amplitude A. Enter A:'
read(*,*) tamp

! Generate spectral amplitude sqrt{k^{2(p-2)} * exp[-(p-1)*(k/k_0)^2]} :
fac=one/rk0**2
p1=p-one
p2=p-two
do ky=1,nw+1
  do kx=1,nw+1
    s=fac*(rk(kx)**2+rk(ky)**2)
    ss(kx,ky)=sqrt(fac*s**p2*exp(-p1*s))
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
    amp=ss(kx,ky)
    ss(kx ,ky )=amp*cx*cy
    ss(kxc,ky )=amp*sx*cy
    ss(kx, kyc)=amp*cx*sy
    ss(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nw
  kxc=ng+2-kx
  phix=twopi*rand(0)-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cx
  ss(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nw
  kyc=ng+2-ky
  phiy=twopi*rand(0)-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cy
  ss(kx, kyc)=amp*sy
enddo

ky=nw+1
do kx=2,nw
  kxc=ng+2-kx
  phix=twopi*rand(0)-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cx
  ss(kxc,ky )=amp*sx
enddo

kx=nw+1
do ky=2,nw
  kyc=ng+2-ky
  phiy=twopi*rand(0)-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cy
  ss(kx, kyc)=amp*sy
enddo

ss(1,1)=zero
ss(nw+1,nw+1)=zero

! Transform to physical space:
call spctop(ng,ng,ss,bb,xfactors,yfactors,xtrig,ytrig)

! Work out rms:
fac=sqrt(dsumi*sum(bb**2))

! Renormalise field:
fac=rms/fac
bb=fac*bb

! Add uniform sloping part in y:
do iy=1,ng
  bb(iy,:)=bb(iy,:)+tamp*sin(gl*dble(iy-1)-pi)
enddo

! Write data to a file:
open(11,file='bt_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,bb
close(11)

end program random_top
