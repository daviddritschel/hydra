program bcjet
! Initialises a baroclinic (mode 2) jet with a Gaussian shape PV profile.
! Use balinit to balance the initial flow using zero first derivatives
! of divergence and acceleration divergence.

use constants

implicit none

double precision:: qa(ng,ng),yc(ng)
double precision:: robc,wid,amp2,amp3,fac,x,y

integer:: i,ix,iy

write(*,*) ' Consider a baroclinic jet with PV = f*A*exp(-((y-y_c)/w)^2)'
write(*,*) ' where f is the Coriolis frequency and'
write(*,*) ' y_c(x) = A_2*sin(2x)+A_3*sin(3x).'
write(*,*)
write(*,*) ' Enter the baroclinic Rossby number A and the width w:'
read(*,*) robc,wid
write(*,*) ' Enter A_2 & A_3:'
read(*,*) amp2,amp3

do ix=1,ng
   x=gl*dble(ix-1)-pi
   yc(ix)=amp2*sin(two*x)+amp3*sin(three*x)
enddo

fac=one/wid**2
robc=cof*robc
do ix=1,ng
   do iy=1,ng
      y=gl*dble(iy-1)-pi
      qa(iy,ix)=robc*exp(-fac*(y-yc(ix))**2)
   enddo
enddo

!-------------------------------------------------------------
! Write PV anomaly in each layer:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
fac=-vm12/(vm11*vm22-vm12*vm21)
write(11,rec=1) zero,fac*qa
fac=-fac*vm11/vm12
write(11,rec=2) zero,fac*qa
close(11)

!-------------------------------------------------------------
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

end program bcjet
