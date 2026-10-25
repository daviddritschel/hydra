program rest

!=============================================================
!  Sets up PV = planetary vorticity
!=============================================================

 !Import constants are parameters:
use constants

 !Declarations:
implicit none

double precision:: qq(ng,nt),qp(ng)
integer:: i,j

!------------------------------------------------------------
 !Define PV for a rest state:
do j=1,ng
  qp(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

do i=1,nt
   qq(:,i)=qq(:,i)+qp
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program rest
