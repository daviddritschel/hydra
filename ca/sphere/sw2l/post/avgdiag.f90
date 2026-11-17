program average

! Computes the time-averaged values of the diagnostics in the file
! evolution/rossby-froude.asc from t_c <= t <= t_max where t_c is specified
! and t_max is the maximum time found in the data.

use parameters

implicit none

double precision:: tc, t, afac
double precision::  romax, frmax, rorms
double precision:: aromax,afrmax,arorms
integer:: iread, na

!---------------------------------------------------------------------
! Open data file:
open(15,file='evolution/rossby-froude.asc',status='old')

write(*,*) 'Time from which to compute the time averages?'
read(*,*) tc
write(*,*)

! Initialise:
aromax=0.d0
afrmax=0.d0
arorms=0.d0
na=0
  
! Read data and process:
do
  iread=0
   !Read Ro_max, Fr_max & Ro_rms:
  read(15,*,iostat=iread) t,romax,frmax,rorms
  if (iread .ne. 0) exit 

  if (t+1.d-6 >= tc) then
    ! Accumulate time average:
    aromax=aromax+romax
    afrmax=afrmax+frmax
    arorms=arorms+zrms
    na=na+1
  endif

enddo

! Finalise averages:
afac=1.d0/dble(na)
aromax=afac*aromax
afrmax=afac*afrmax
arorms=afac*arorms

! Close files:
close(51)
close(61)

write(*,'(a,f12.5)') ' The final time found is ',t
write(*,*)
write(*,*) ' The time averaged diagnotics are as follows:'
write(*,*)
write(*,*) '     Ro_max      Fr_max      Ro_rms'
write(*,'(3(1x,f11.6))') aromax,afrmax,arorms

! Record info:
open(16,file='evolution/avgdiag.asc',status='replace')
write(16,'(2(1x,f10.3),3(1x,f11.6))') tc,t,aromax,afrmax,arorms

! End main program
end program average
