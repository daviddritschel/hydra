program ellipole

!===============================================================
!  Sets up a PV anomaly in the shape of an elliptical Gaussian
!  centred on the north pole.

!  Developed to test the initial balancing procedure.
!===============================================================

 !Import constants are parameters:
use constants

 !Declarations:
implicit none

double precision:: qq(ng,nt),qp(ng)
double precision:: amp,asp,aspi,b,bsq
double precision:: lon,lat,x,y,z
integer:: i,j

!------------------------------------------------------------
 !Define planetary vorticity (qp):
do j=1,ng
  qp(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

 !Define the PV anomaly:
write(*,*) ' The PV anomaly has the form'
write(*,*) '   A*exp(-(x^2*lambda+y^2/lambda)/(4*b^2*(1-z)).'
write(*,*)
write(*,*) ' Enter A, lambda and b:'
read(*,*) amp,asp,b

aspi=one/asp
bsq=b**2

do i=1,nt
   lon=dl*dble(i-1)-pi
   do j=1,ng
      lat=(dble(j)-f12)*dl-hpi
      x=cos(lat)*cos(lon)
      y=cos(lat)*sin(lon)
      z=sin(lat)
      qq(j,i)=amp*exp(-(asp*x**2+aspi*y**2)/(bsq*(one+z)))
   enddo
enddo

 !Add planetary vorticity to the anomaly in qq:
do i=1,nt
   qq(:,i)=qq(:,i)+qp
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program ellipole
