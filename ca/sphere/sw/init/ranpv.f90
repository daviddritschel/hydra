program ranpv

!=====================================================
! Sets up a random PV anomaly field using the forcing
! function in module force.f90.
!=====================================================

 !Import necessary modules:
use randomise

 !Declarations:
implicit none

double precision:: qq(ng,nt),qp(ng),qmax
integer:: ordbeg,ordend,i,j

!------------------------------------------------------------
write(*,*) ' Enter the range in order, n_beg to n_end:'
read(*,*) ordbeg,ordend

write(*,*) ' Enter the maximum absolute amplitude of the PV anomaly:'
read(*,*) qmax

!Generate random PV anomaly field:
call generate_random_field(qq, qmax, ordbeg, ordend)

!Add resting-state PV, 2*Omega*sin(latitude):
do j=1,ng
  qp(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

do i=1,nt
   qq(:,i)=qq(:,i)+qp
enddo

!Write total initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program ranpv
