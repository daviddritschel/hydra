module congen

! Module containing subroutines for rebuilding PV contours given either
! the PV field and no contours, or the PV residual field together with
! contours.  This creates new contours in either case.

! Revised by D G Dritschel on 31 July 2020 to use any number of layers.

use common

implicit none

 !Array for storing the PV field interpolated to the ultra-fine grid: 
double precision:: qa(ngu+1,ngu)

 !Temporary arrays for contour storage:
double precision:: xa(npm),ya(npm)
integer:: nexta(npm),inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta
 !Note: contours are built level by level to keep storage to a minimum.


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=====================================================================

subroutine recontour(qg)
! Main routine for recontouring (from Dritschel & Ambaum, 1996, QJRMS).

! qg           : a gridded field added to that due to contours
! xq(i),yq(i)  : location of node i in the domain
! nextq(i)     : index of the node following node i 
! layq(j)      : layer (integer) containing contour j
! indq(j)      : field level (integer) of contour j
! npq(j)       : number of nodes on contour j
! i1q(j)       : beginning node index on contour j
! i2q(j)       : ending node index on contour j
! jl1q(iz)     : beginning contour index in layer iz
! jl2q(iz)     : ending contour index in layer iz
! il1q(iz)     : beginning node index in layer iz
! il2q(iz)     : ending node index in layer iz
! nq           : number of contours
! nptq         : total number of nodes

implicit none

 !Passed array:
double precision:: qg(ng,ng,nz)

 !Local variables:
integer:: ix,ixf,ix0,ix1
integer:: iy,iyf,iy0,iy1
integer:: i,j,iz
logical:: create

 !-----------------------------------------------------------------
 !Counters for total number of nodes and contours:
npta=0
na=0

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !              Begin a major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,nz
   !Obtain ultra-fine grid field qa in layer iz:
  create=jl2q(iz) > 0
  if (create) then
     !This means that there are contours in this layer; convert them 
     !to gridded values (in the array qa):
    call con2ugrid(iz)

     !Bi-linear interpolate the residual qg to the fine grid and add to qa:
    do ix=1,ngu
      ixf=igfw(ix)
      ix0=ig0w(ix)
      ix1=ig1w(ix)

      do iy=1,ngu
        iyf=igfw(iy)
        iy0=ig0w(iy)
        iy1=ig1w(iy)

        qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qg(iy0,ix0,iz) &
                           +w10(iyf,ixf)*qg(iy1,ix0,iz) &
                           +w01(iyf,ixf)*qg(iy0,ix1,iz) &
                           +w11(iyf,ixf)*qg(iy1,ix1,iz)

      enddo
    enddo

    ! Compute new contour interval
    qjump(iz)=(maxval(qa)-minval(qa))/dble(ncontq)

  else
    !There are no contours (the usual situation at t = 0):


    ! Compute new contour interval
    qjump(iz)=(maxval(qq(:,:,iz))-minval(qq(:,:,iz)))/dble(ncontq)
    create=qjump(iz) > small

    if (create) then 

       !Interpolate qg (here, the full field) to the fine grid as qa:
      do ix=1,ngu
        ixf=igfw(ix)
        ix0=ig0w(ix)
        ix1=ig1w(ix)

        do iy=1,ngu
          iyf=igfw(iy)
          iy0=ig0w(iy)
          iy1=ig1w(iy)

          qa(iy,ix)=w00(iyf,ixf)*qg(iy0,ix0,iz) &
                   +w10(iyf,ixf)*qg(iy1,ix0,iz) &
                   +w01(iyf,ixf)*qg(iy0,ix1,iz) &
                   +w11(iyf,ixf)*qg(iy1,ix1,iz)
        enddo
      enddo

    else
      qjump(iz)=zero
    endif
  endif

   !See if there are contours to create in this layer:
  if (create) then

     !Reset starting contour and node indices in each layer:
    jl1q(iz)=na+1
    il1q(iz)=npta+1

     !Generate new contours (xa,ya) from qa array:
    call ugrid2con(iz)

     !Reset ending contour and node indices in each layer:
    jl2q(iz)=na
    il2q(iz)=npta
  endif

enddo

 !Copy contour arrays and indices back to those in the main code:
do i=1,npta
  xq(i)=xa(i)
  yq(i)=ya(i)
  nextq(i)=nexta(i)
enddo

do j=1,na
  i1q(j)=i1a(j)
  i2q(j)=i2a(j)
  npq(j)=npa(j)
  indq(j)=inda(j)
enddo

do iz=1,nz
  if (jl2q(iz) > 0) then
    do j=jl1q(iz),jl2q(iz)
      layq(j)=iz
    enddo
  endif
enddo

nq=na
nptq=npta

return
end subroutine recontour

!==========================================================================

subroutine ugrid2con(iz)
! Generates contours (xa,ya) from the gridded field qa for the levels
! +/-qjump/2, +/-3*qjump/2, ....

implicit none

 !Passed index:
integer:: iz

 !Local parameters and arrays:
integer,parameter:: ncrm=3*nplm/4
 !ncrm:  maximum number of contour crossings of a single contour level
 !nplm:  maximum number of nodes in any contour level
 
integer(kind=dbleint),parameter:: &
     ngsq=int(ngu,kind=dbleint)*int(ngu,kind=dbleint)
integer(kind=dbleint):: k,kob,kaa,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(ngu+1),qdy(ngu+1)
double precision:: xd(nprm),yd(nprm)
double precision:: qoff,dq,dqi,qtmp,xgt,ygt,xx,yy
integer:: isx(ngu+1),isy(ngu+1),icre(nm),icrtab(ngsq,2)
integer:: levbeg,levend,lev,levt,noc,icrn
integer:: ncr,i,ix,iy,icr,inc,ibeg,iend,npd
integer(kind=halfint):: noctab(ngsq)
logical:: free(ncrm),keep

 !----------------------------------------
 !Initialise various constants used below:
dq=qjump(iz)
dqi=one/dq
qoff=dq*dble(nlevm)
 !qoff: should be a large integer multiple of the contour interval, dq.  
 !The multiple should exceed the maximum expected number of contour levels.

 !--------------------------------------------------
 !First get the beginning and ending contour levels:
levbeg=int((qoff+minval(qa(1:ngu,:)))*dqi+f12)+1
levend=int((qoff+maxval(qa(1:ngu,:)))*dqi+f12)

if (levbeg <= levend) then
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Loop over contour levels and process:
do lev=levbeg,levend
   !Integer index giving contour level:
  levt=lev-nlevm+(lev-1)/nlevm-1

   !Counter for total number of grid line crossings:
  ncr=0

   !Contour level being sought:
  qtmp=(dble(lev)-f12)*dq-qoff

   !Below, kib = grid box into which the contour (containing ncr) is going
   !       kob =   "   "  out of "    "     "         "       "    " coming
   !      [kob -> ncr -> kib:  ncr lies at the boundary between kob & kib]

   !*** grid boxes are numbered 1 (lower left) to ngu*ngu (upper right) ***

   !Initialise number of crossings per box:
  noctab=0

   !-----------------------------------------------------------
   !Find x grid line crossings first:
  do ix=1,ngu
    xgt=xgu(ix)

    qdy(1:ngu)=qa(1:ngu,ix)-qtmp
    isy(1:ngu)=sign(one,qdy(1:ngu))
    qdy(ngu+1)=qdy(1)
    isy(ngu+1)=isy(1)

    do iy=1,ngu
      if (isy(iy) /= isy(iy+1)) then
        ncr=ncr+1
        inc=(1-isy(iy))/2
        kaa=(iy-1)*ngu+1
        kib(ncr)=kaa+ibg(ix+inc)
        kob=kaa+ibg(ix+1-inc)
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        xcr(ncr)=xgt
        yy=ygu(iy)-glu*qdy(iy)/(qdy(iy+1)-qdy(iy))
        ycr(ncr)=oms*(yy-twopi*dble(int(yy*pinv)))
      endif
    enddo

  enddo

   !----------------------------------------------------------
   !Find y grid line crossings next:
  do iy=1,ngu
    ygt=xgu(iy)

    do ix=1,ngu
      qdx(ix)=qa(iy,ix)-qtmp
      isx(ix)=sign(one,qdx(ix))
    enddo
    qdx(ngu+1)=qdx(1)
    isx(ngu+1)=isx(1)

    do ix=1,ngu
      if (isx(ix) /= isx(ix+1)) then
        ncr=ncr+1
        inc=(1-isx(ix))/2
        kib(ncr)=ngu*ibg(iy+1-inc)+ix
        kob=ngu*ibg(iy+inc)+ix
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        ycr(ncr)=ygt
        xx=xgu(ix)-glu*qdx(ix)/(qdx(ix+1)-qdx(ix))
        xcr(ncr)=oms*(xx-twopi*dble(int(xx*pinv)))
      endif
    enddo

  enddo

   !----------------------------------------------------------------
   !Now re-build contours:
  free(1:ncr)=.true.

  do icr=1,ncr
    if (free(icr)) then
       !A new contour (indexed na) starts here:
      na=na+1
      inda(na)=levt
      ibeg=npta+1
      i1a(na)=ibeg

       !First point on the contour:
      npd=1
      xd(1)=xcr(icr)
      yd(1)=ycr(icr)

       !Find remaining points on the contour:
      k=kib(icr)
       !k is the box the contour is entering
      noc=noctab(k)
       !Use last crossing (noc) in this box (k) as the next node:
      icrn=icrtab(k,noc)
       !icrn gives the next point after icr (icrn is leaving box k)
      do while (icrn /= icr)
        noctab(k)=noc-1
         !noctab is usually zero now except for boxes with a
         !maximum possible 2 crossings
        npd=npd+1
        xd(npd)=xcr(icrn)
        yd(npd)=ycr(icrn)
        free(icrn)=.false.
        k=kib(icrn)
        noc=noctab(k)
        icrn=icrtab(k,noc)
      enddo

       !Re-distribute nodes on this contour 3 times to reduce complexity:
      keep=.false.
      do
        call renode(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
         !Delete contour if deemed too small (see renode in contours.f90):
        if (npa(na) == 0) exit
        call renode(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
         !Delete contour if deemed too small:
        if (npd == 0) exit
        call renode(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
         !Delete contour if deemed too small:
        if (npa(na) == 0) exit
         !Contour is big enough to keep:
        keep=.true.
        exit
      enddo

      if (keep) then 
        npta=npta+npa(na)
        iend=ibeg+npa(na)-1
        i2a(na)=iend
        do i=ibeg,iend-1
          nexta(i)=i+1
        enddo
        nexta(iend)=ibeg
      else
        na=na-1
      endif

      free(icr)=.false.
    endif
  enddo

enddo
 !End of loop over contour levels
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
endif

return
end subroutine ugrid2con

!=======================================================================

subroutine con2ugrid(iz)
! Contour -> ultra-fine grid conversion in layer iz.  

implicit none

 !Local parameters and arrays:
double precision:: qjx(ngu)
double precision:: dx(nptq),dy(nptq)
double precision:: dq,sdq,xx,yy,cc,p,px0,qavg0,qadd
integer:: ixc(nptq),ngc(nptq)
integer:: iz,i,ia,ix,ixdif,ixbeg,iy,jump,ncr
logical:: crossx(nptq)

!-------------------------------------------------------------------
! PV jump in this layer:
dq=qjump(iz)

 !Initialise interior x grid line crossing information and fill the
 !PV (q) jump array along the lower boundary (iy = 1 or y = -pi):
do i=il1q(iz),il2q(iz)
  ixc(i)=int(glui*(xq(i)+pi))
enddo

qjx=zero

do i=il1q(iz),il2q(iz)
  ia=nextq(i)
  xx=xq(ia)-xq(i)
  dx(i)=xx-twopi*dble(int(xx*pinv))
  yy=yq(ia)-yq(i)
  dy(i)=yy-twopi*dble(int(yy*pinv))
  ixdif=ixc(ia)-ixc(i)
  ngc(i)=ixdif-ngu*((2*ixdif)/ngu)
  crossx(i)=(ngc(i) /= 0)
  if (abs(yy) > pi) then
     !The contour segment (i,ia) crosses y = -pi; find x location:
    cc=sign(one,dy(i))
    p=-(yq(i)-pi*cc)/(dy(i)+small)
    ix=1+int(glui*(mod(thrpi+xq(i)+p*dx(i),twopi)))
    qjx(ix)=qjx(ix)-dq*cc
     !Note: qjx gives the jump in q going from ix to ix+1
  endif
enddo

 !Sum q jumps to obtain the gridded field qa along lower boundary:
qa(1,1)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=1,ngu-1
  qa(1,ix+1)=qa(1,ix)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
qa(2:ngu+1,:)=zero

 !Determine x grid line crossings and accumulate q jumps:
do i=il1q(iz),il2q(iz)
  if (crossx(i)) then
    jump=sign(1,ngc(i))
    ixbeg=ixc(i)+(1+jump)/2+ngu
    sdq=dq*dble(jump)
    ncr=0
    do while (ncr /= ngc(i)) 
      ix=1+mod(ixbeg+ncr,ngu)
      xx=xgu(ix)-xq(i)
      px0=(xx-twopi*dble(int(xx*pinv)))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      yy=yq(i)+px0*dy(i)
      iy=2+int(glui*(yy-twopi*dble(int(yy*pinv))+pi))
       !Increment q jump between the grid lines iy-1 & iy:
      qa(iy,ix)=qa(iy,ix)+sdq
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=1,ngu
  do iy=2,ngu
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

 !Remove domain average:
qavg0=sum(qa(1:ngu,:))/dble(ngu*ngu)
qadd=qavg(iz)-qavg0
qa(1:ngu,:)=qa(1:ngu,:)+qadd

return
end subroutine con2ugrid

!==========================================================================

 !Main end module
end module congen
