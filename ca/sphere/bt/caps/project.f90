program project
!  -------------------------------------------------------------------------
!  |   Projects the vorticity on the spherical harmonics of order 2:       |
!  |                                                                       |
!  |   a0*Y_2^0 + a1*Y_2^1-conj(a1)*Y_2^{-1} + a2*Y_2^2+conj(a2)*Y_2^{-2}  |
!  |                                                                       |
!  |   where "conj" means complex conjugate.                               |
!  |                                                                       |
!  |   Output is to the formatted file "coeffs.asc" which lists time vs    |
!  |   the coefficients a0, Re(a1), Im(a1), Re(a2) & Im(a2).               |
!  -------------------------------------------------------------------------

use constants
implicit none

real:: tr4,zzr4(ng,nt)
double precision:: zz(ng,nt)
double precision:: clon(nt),slon(nt)
double precision:: clat(ng),slat(ng),rdt(ng),y20(ng)
double precision:: y21c(ng,nt),y21s(ng,nt)
double precision:: y22c(ng,nt),y22s(ng,nt)
double precision:: c20,c21,c22,w1,w2
double precision:: a0,a1r,a1i,a2r,a2i
double precision:: aspsqm1,rlon,rlat,x,y
integer:: i,j,loop,iread

!---------------------------------------------------------------
 !Define various fixed arrays:
do i=1,nt
  rlon=dble(i-1)*dl-pi
  clon(i)=cos(rlon)
  slon(i)=sin(rlon)
enddo

do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

 !Define dmu/dtheta = rho*tau:
aspsqm1=asp**2-one
do j=1,ng
  rdt(j)=clat(j)*sqrt(one+aspsqm1*clat(j)**2)
enddo

 !Spherical harmonics of order 2 (real forms):
c20=f14*sqrt(five/pi)
do j=1,ng
  y20(j)=c20*(three*slat(j)**2-one)
enddo
c21=sqrt(15.d0/pi)
c22=f12*c21
do i=1,nt
  do j=1,ng
    x=clat(j)*clon(i)
    y=clat(j)*slon(i)
    y21c(j,i)=-c21*x*slat(j)
    y21s(j,i)=c21*y*slat(j)
    y22c(j,i)=c22*(x**2-y**2)
    y22s(j,i)=-c21*x*y
  enddo
enddo

 !Normalisation coefficients used below in projection:
c20=f1112*(y20(1)**2*rdt(1)+y20(ng)**2*rdt(ng))
do j=2,ngm1
  c20=c20+y20(j)**2*rdt(j)
enddo
c20=one/(c20*dble(nt))

c21=zero
c22=zero
do i=1,nt
  c21=c21+f1112*(y21c(1,i)**2*rdt(1)+y21c(ng,i)**2*rdt(ng))
  c22=c22+f1112*(y22c(1,i)**2*rdt(1)+y22c(ng,i)**2*rdt(ng))
  do j=2,ngm1
    c21=c21+y21c(j,i)**2*rdt(j)
    c22=c22+y22c(j,i)**2*rdt(j)
  enddo
enddo
c21=sqrt(two)/c21
c22=sqrt(two)/c22

!---------------------------------------------------------------
 !Open input data file:
open(44,file='zz.r4',form='unformatted', & 
    & access='direct',status='old',recl=nbytes)

 !Open output file:
open(22,file='coeffs.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(44,rec=loop,iostat=iread) tr4,zzr4
  if (iread .ne. 0) exit 

   !Convert to double precision:
  zz=dble(zzr4)

   !Project onto spherical harmonics:
  a0=zero
  a1r=zero
  a1i=zero
  a2r=zero
  a2i=zero
  do i=1,nt
    w1=zz(1,i)*rdt(1)
    w2=zz(ng,i)*rdt(ng)
    a0=a0+f1112*(y20(1)*w1+y20(ng)*w2)
    a1r=a1r+f1112*(y21c(1,i)*w1+y21c(ng,i)*w2)
    a1i=a1i+f1112*(y21s(1,i)*w1+y21s(ng,i)*w2)
    a2r=a2r+f1112*(y22c(1,i)*w1+y22c(ng,i)*w2)
    a2i=a2i+f1112*(y22s(1,i)*w1+y22s(ng,i)*w2)
    do j=2,ngm1
      w1=zz(j,i)*rdt(j)
      a0=a0+y20(j)*w1
      a1r=a1r+y21c(j,i)*w1
      a1i=a1i+y21s(j,i)*w1
      a2r=a2r+y22c(j,i)*w1
      a2i=a2i+y22s(j,i)*w1
    enddo
  enddo
  a0=a0*c20
  a1r=a1r*c21
  a1i=a1i*c21
  a2r=a2r*c22
  a2i=a2i*c22

  write(22,'(f12.5,5(1x,f12.8))') tr4,a0,a1r,a1i,a2r,a2i
  write( *,'(a,f12.5)') ' Processed t = ',tr4

enddo

 !Close all files:
close(44)
close(22)

write(*,*)
write(*,*) ' t vs a0, Re(a1), Im(a1), Re(a2) & Im(a2) is listed in coeffs.asc;'
write(*,*) ' to view the results type'
write(*,*)
write(*,*) ' plotcol coeffs.asc'
write(*,*)

end program
