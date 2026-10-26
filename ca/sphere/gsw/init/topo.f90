program topo

!=============================================================
!  Sets up a random bottom topography field
!  using a function in module randomise.f90
!=============================================================

 !Import necessary modules:
use randomise

 !Declarations:
implicit none

double precision:: bb(ng,nt),bmax
integer:: ordbeg,ordend

!------------------------------------------------------------
write(*,*) ' Enter the range in order, n_beg to n_end:'
read(*,*) ordbeg,ordend

write(*,*) ' Enter the maximum absolute amplitude of the topography'
write(*,*) ' scaled on the mean layer depth:'
read(*,*) bmax

 !Generate random topography scaled on the mean depth:
call generate_random_field(bb, bmax, ordbeg, ordend)

 !Write data:
open(11,file='bb_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,bb
close(11)

end program topo
