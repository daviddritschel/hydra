program ini_ortho
!  ---------------------------------------------------------------------
!  |   Creates an orthographic projection of the data in ff_init.r8    |
!  |   ff is any selected field in the main job directory.             |
!  |   Writes ffopi.r4 which can be imaged using odv.                  |
!  ---------------------------------------------------------------------

use constants

implicit none

integer,parameter:: nfield=5
real:: rlatc,rlonc
integer:: ny,nz,nbprec
integer:: iopt
character(len=2):: dname(nfield),field(nfield)
character(len=10):: infile
character(len=8):: outfile
character(len=5):: viewfile

data dname /'qq','zz','dd','gg','hh'/
data field /'qq','zz','dd','gg','hh'/

!---------------------------------------------------------------------
write(*,*) ' Which field do you wish to project?'
write(*,*) '   (1)  PV (q),'
write(*,*) '   (2)  relative vorticity (zeta),'
write(*,*) '   (3)  velocity divergence (delta),'
write(*,*) '   (4)  acceleration divergence (gamma), or'
write(*,*) '   (5)  dimensionless height anomaly (h)?'
write(*,*) ' Option?'
read(*,*) iopt
if (iopt*(nfield+1-iopt) .le. 0) then
  write(*,*) ' *** Not a valid option; exiting!'
  stop
endif

infile=dname(iopt)//'_init.r8'
outfile=field(iopt)//'opi.r4'
viewfile=field(iopt)//'opi'

write(*,*) ' Latitude & longitude of the direction of view (degrees)?'
read(*,*) rlatc,rlonc

!write(*,*) ' Width (in pixels) of output image?'
!read(*,*) ny
ny=ng !need this to be compatible with the odv viewing script [could use ny=nt]
nz=ny

 !Define output record length (in bytes):
nbprec=4*(ny*nz+1)

 !Open input data file:
open(11,file=infile,form='unformatted',access='direct', &
                   status='old',recl=2*nbytes)

 !Open output data file:
open(22,file='evolution/'//outfile,form='unformatted',access='direct', &
                                 status='replace',recl=nbprec)

 !Project data orthographically:
call opgrid(rlatc,rlonc,ny,nz)

 !Close files:
close(11)
close(22)

write(*,*)
write(*,*) ' *** View the results using'
write(*,*)
write(*,*) ' odv '//viewfile
write(*,*)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!============================================================================

subroutine opgrid(rlatc,rlonc,ny,nz)
! Projects lat-lon data in zz.r4 orthographically on an ellipsoid

implicit none

integer:: ny,nz

integer:: i,j
integer:: iy,iz,ic,ip1,jp1

double precision:: tr,qr(ng,nt)
real:: qq(0:ngp1,nt),qqop(nz,ny),yg(ny),zg(nz)
real:: qqbg,t,dy,dz,det
real:: rlatc,clatc,slatc,rlonc,clonc,slonc
real:: xp,yp,zp,xm,xt,yt,zt,ri,rj,aa,bb,cc,dd,qqt

!------------------------------------------------------------------
!write(*,*) ' Background field level (for points off the surface)?'
!read(*,*) qqbg
qqbg=0.0

!------------------------------------------------------------------
 !Define points (yg,zg) in cross-sectional plane of view:
dy=2.0/float(ny)
do iy=1,ny
  yg(iy)=dy*(float(iy)-0.5)-1.0
enddo
dz=dy
do iz=1,nz
  zg(iz)=dz*(float(iz-nz/2)-0.5)
enddo

clatc=cos(rlatc*pi/180.)
slatc=sin(rlatc*pi/180.)

clonc=cos(rlonc*pi/180.)
slonc=sin(rlonc*pi/180.)

!------------------------------------------------------------------
read(11,rec=1) tr,qr

t=real(tr)
qq(1:ng,:)=real(qr)

 !Copy latitudes adjacent to poles (j = 1 and ng) with a pi
 !shift in longitude to simplify interpolation below:
do i=1,ng
  ic=i+ng
  qq(0,i)=qq(1,ic)
  qq(0,ic)=qq(1,i)
  qq(ngp1,i)=qq(ng,ic)
  qq(ngp1,ic)=qq(ng,i)
enddo

 !Do interpolation in chosen perspective:
do iy=1,ny
  yp=yg(iy)
  do iz=1,nz
    zp=zg(iz)
    det=1.-yp**2-zp**2
    if (det .gt. zero) then
      xp=sqrt(det)
      xm=xp*clatc-zp*slatc
      zt=zp*clatc+xp*slatc
      xm=xp*clatc-zp*slatc
      yt=yp*clonc+xm*slonc
      xt=xm*clonc-yp*slonc

       !Find longitude & latitude then bi-linearly interpolate qq:
      ri=dli*(pi+atan2(yt,xt))
      i=1+int(ri)
      ip1=1+mod(i,nt)
      bb=float(i)-ri
      aa=one-bb

      rj=dli*(hpidl+asin(zt))
      j=int(rj)
      jp1=j+1
      cc=rj-float(j)
      dd=one-cc

      qqop(iz,iy)=bb*(dd*qq(j,i)  +cc*qq(jp1,i))   &
                 +aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
    else
      qqop(iz,iy)=qqbg
    endif
  enddo
enddo

write(22,rec=1) qqop

end subroutine opgrid

!=========================================================

end program ini_ortho
