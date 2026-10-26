program ranhh

!=====================================================
! Sets up a random dimensionless height anomaly field
! using a function in module randomise.f90

! Both zeta and delta = 0; gamma is obtained after by
! a call to inigamma. See flow-setup script.  
!=====================================================

 !Import necessary modules:
use randomise

 !Declarations:
implicit none

double precision:: hh(ng,nt),qp(ng),hmax
integer:: ordbeg,ordend,i,j

!------------------------------------------------------------
write(*,*) ' Enter the range in order, n_beg to n_end:'
read(*,*) ordbeg,ordend

write(*,*) ' Enter the maximum absolute amplitude of the dimensionless:'
write(*,*) ' height anomaly:'
read(*,*) hmax

!Generate random height anomaly field:
call generate_random_field(hh, hmax, ordbeg, ordend)

!Write data:
open(20,file='hh_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

!Define resting-state PV, 2*Omega*sin(latitude)/(1+hh):
do j=1,ng
   qp(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

!Re-use hh array for total PV:
do i=1,nt
   hh(:,i)=qp/(one+hh(:,i))
enddo

!Write data:
open(20,file='qq_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

!Write zero divergence:
hh=zero
open(20,file='dd_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

end program ranhh
