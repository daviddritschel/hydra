program random_huv

!=========================================================
!  Sets up random h, u & v fields to test the inversion
!  procedure. Uses a function in module randomise.f90
!=========================================================

 !Import necessary modules:
use spectral
use randomise

 !Declarations:
implicit none

double precision:: h1(ng,nt), h2(ng,nt)
double precision:: d1(ng,nt), d2(ng,nt)
double precision:: g1(ng,nt), g2(ng,nt)
double precision:: z1(ng,nt), z2(ng,nt)
double precision:: q1(ng,nt), q2(ng,nt)
double precision:: u1(ng,nt), u2(ng,nt)

double precision:: wka(ng,nt), wkb(ng,nt)
double precision:: qmax
integer:: ordbeg, ordend, i, j, m

!------------------------------------------------------
! Initialise inversion constants and arrays:
call init_spectral

!------------------------------------------------------
! Generate random field for h_1:
qmax = 0.2d0
ordbeg = 1
ordend = 8
call generate_random_field(h1, qmax, ordbeg, ordend)

! Generate random field for h_2:
qmax = 0.1d0
ordbeg = 3
ordend = 12
call generate_random_field(h2, qmax, ordbeg, ordend)

!------------------------------------------------------
! Generate random field for d_1:
qmax = 0.4d0
ordbeg = 2
ordend = 10
call generate_random_field(d1, qmax, ordbeg, ordend)

! Generate random field for d_2:
qmax = 0.8d0
ordbeg = 3
ordend = 8
call generate_random_field(d2, qmax, ordbeg, ordend)

!------------------------------------------------------
! Generate random field for z_1:
qmax = 5.0d0
ordbeg = 4
ordend = 6
call generate_random_field(z1, qmax, ordbeg, ordend)

! Ensure zero domain average:
qmax = average(z1)
z1 = z1 - qmax

! Generate random field for z_2:
qmax = 3.0d0
ordbeg = 1
ordend = 12
call generate_random_field(z2, qmax, ordbeg, ordend)

! Ensure zero domain average:
qmax = average(z2)
z2 = z2 - qmax

!------------------------------------------------------
! Compute PV in each layer:
do i=1,nt
   q1(:,i) = (z1(:,i)+cof)/(one+h1(:,i))
   q2(:,i) = (z2(:,i)+cof)/(one+h2(:,i))
enddo

! Write PV field:
open(11,file = 'qq_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,q1
write(11,rec=2) zero,q2
close(11)

! Write h:
open(11,file = 'hh0_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,h1
write(11,rec=2) zero,h2
close(11)

! Write delta:
open(11,file = 'dd0_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,d1
write(11,rec=2) zero,d2
close(11)

! Write zeta:
open(11,file = 'zz0_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,z1
write(11,rec=2) zero,z2
close(11)

!------------------------------------------------------
! Compute u & gamma in each layer:

! First convert the height field in each layer to semi-spectral space:
call forfft(ng,nt,h1,trig,factors) 
call forfft(ng,nt,h2,trig,factors) 

! ==> Layer 1 <==

! Convert the relative vorticity to semi-spectral space:
call forfft(ng,nt,z1,trig,factors) 

! Invert Laplace's operator on relative vorticity:
call laplinv(z1,wka,u1)
! Here the streamfunction is wka while d(wka)/dlat = u_1.

! Flip sign of u_1 (in semi-spectral space):
u1 = -u1

! Get pressure term:
wka = csq1*h1+denrat*csq2*h2
call laplace(wka,wkb)

! Add things up to define gamma_1 (in semi-spectral space):
do m=1,nt
  g1(:,m)=cof*z1(:,m)-bet*u1(:,m)-wkb(:,m)
enddo

! Convert u_1 & gamma_1 to physical space for output below:
call revfft(ng,nt,u1,trig,factors) 
call revfft(ng,nt,g1,trig,factors) 

! ==> Layer 2 <==

! Convert the relative vorticity to semi-spectral space:
call forfft(ng,nt,z2,trig,factors) 

! Invert Laplace's operator on relative vorticity:
call laplinv(z2,wka,u2)
! Here the streamfunction is wka while d(wka)/dlat = u_2.

! Flip sign of u_2 (in semi-spectral space):
u2 = -u2

! Get pressure term:
wka = csq1*h1+csq2*h2
call laplace(wka,wkb)

! Add things up to define gamma_2 (in semi-spectral space):
do m=1,nt
  g2(:,m)=cof*z2(:,m)-bet*u2(:,m)-wkb(:,m)
enddo

! Convert u_2 & gamma_2 to physical space for output below:
call revfft(ng,nt,u2,trig,factors) 
call revfft(ng,nt,g2,trig,factors) 

!------------------------------------------------------
! Write u:
open(11,file = 'uu0_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,u1
write(11,rec=2) zero,u2
close(11)

! Write gamma:
open(11,file = 'gg0_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,g1
write(11,rec=2) zero,g2
close(11)

end program random_huv
