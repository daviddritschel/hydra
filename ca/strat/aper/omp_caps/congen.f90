module congen
! Converts contours (xq,yq) to gridded values on an ultra-fine grid 
! of dimensions mgu*nx x mgu*ny (in a closed rectangular domain), 
! optionally adds a residual field (interpolated to the ultra-fine 
! grid), and then creates new contours.

! Open contours originating and terminating in a boundary added by
! D G Dritschel on 18 June 2012 @ Moscow

use common
use generic

implicit none

double precision:: qa(0:nyu+1,0:nxu)
double precision:: xa(npm),ya(npm)
integer:: inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta

contains

!=====================================================================

subroutine recontour(qq,xq,yq,dq,qavg,nextq,indq,npq,i1q,i2q,nq,nptq,iopt)
! Main routine for recontouring (from D & Ambaum, 1996, QJRMS)

! qq           : a gridded field added to that due to contours if iopt = 1
! xq(i),yq(i)  : location of node i in the domain
! dq           : contour interval
! qavg         : average value of field (computed if nptq=0)
! nextq(i)     : index of the node following node i 
!                *** this must be zero for an endpoint on the boundary ***
! indq(j)      : field level (integer) of contour j
! npq(j)       : number of nodes on contour j
! i1q(j)       : beginning node index on contour j
! i2q(j)       : ending node index on contour j
! nq           : number of contours
! nptq         : total number of nodes
! iopt         : if 1, always combine a residual gridded field with
!                that due to any contours for recontouring;
!                if 0, and if nq = 0, qq is assumed to contain the 
!                full field to be contoured (this is normal at t = 0).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nx)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm),i2q(nm)

 !Local quantities:
logical:: contours,residual

 !-------------------------------------------------------------------
 !First check if there are any contours to convert to gridded values:
contours=nq .gt. 0

 !See if there is a residual field to add:
residual=(.not. contours) .or. (iopt .eq. 1)

 !-----------------------------------------------------------------
 !Counters for total number of nodes and contours:
npta=0
na=0

 !Obtain fine grid field qa:
if (contours) then
   !Convert contours to gridded values (in the array qa):
  call con2ugrid(xq,yq,dq,qavg,nextq,nptq)

  if (residual) then
     !Bi-linear interpolate the residual qq to the fine grid and add to qa:
    do ix=0,nxu
      ixf=ixfw(ix)
      ix0=ix0w(ix)
      ix1=ix1w(ix)

      do iy=0,nyu
        iyf=iyfw(iy)
        iy0=iy0w(iy)
        iy1=iy1w(iy)

        qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0) &
                         & +w10(iyf,ixf)*qq(iy1,ix0) &
                         & +w01(iyf,ixf)*qq(iy0,ix1) &
                         & +w11(iyf,ixf)*qq(iy1,ix1)

      enddo
    enddo
  endif

else

   !Check if field requires contouring by computing l1 norm of qq:
  call l1norm(qq,qql1)
  if (qql1 .lt. small) then
    qavg=zero
    return
  endif

   !Compute average value of the field (qavg):
  call average(qq,qavg)

   !No contours: interpolate qq (which here contains the full field)
   !to the fine grid as qa:
  do ix=0,nxu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=0,nyu
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=w00(iyf,ixf)*qq(iy0,ix0)+w10(iyf,ixf)*qq(iy1,ix0) &
             & +w01(iyf,ixf)*qq(iy0,ix1)+w11(iyf,ixf)*qq(iy1,ix1)
 
    enddo
  enddo

endif

 !Generate new contours (xa,ya) from qa array:
call ugrid2con(dq,nextq)

 !Copy arrays back to those in the argument of the subroutine:
do i=1,npta
  xq(i)=xa(i)
  yq(i)=ya(i)
enddo

do j=1,na
  i1q(j)=i1a(j)
  i2q(j)=i2a(j)
  npq(j)=npa(j)
  indq(j)=inda(j)
enddo

nq=na
nptq=npta

return
end subroutine

!==========================================================================

subroutine ugrid2con(dq,nextq)
! Generates contours (xa,ya) from the gridded field qa for the levels
! +/-dq/2, +/-3*dq/2, ....

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
integer:: nextq(npm)

 !Local parameters and arrays:
integer,parameter:: ncrm=3*nplm/4
 !ncrm:  max number of contour crossings of a single contour level
 !nplm:  max number of nodes in any contour level
 
integer,parameter:: nxny=nxu*nyu, koff=nxu*(nyu-1)
 
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(0:nxu),qdy(0:nyu)
double precision:: xd(nprm),yd(nprm)
integer:: isx(0:nxu),isy(0:nyu)
integer:: kib(ncrm),icre(nm)
integer:: icrtab(nxny,2)
integer*1:: noctab(nxny)
logical:: free(ncrm),keep

 !initialise constants and arrays:
dqi=one/dq
qoff=dq*dble(nlevm)
 !qoff: should be a large integer multiple of the contour interval, dq.  
 !The multiple should exceed the maximum expected number of contour levels.

 !--------------------------------------------------------
 !First get the beginning and ending contour levels:
qamax=qa(0,0)
qamin=qa(0,0)
do ix=0,nxu
  do iy=0,nyu
    qamax=max(qamax,qa(iy,ix))
    qamin=min(qamin,qa(iy,ix))
  enddo
enddo

levbeg=int((qoff+qamin)*dqi+f12)+1
levend=int((qoff+qamax)*dqi+f12)

if (levbeg .le. levend) then
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Loop over contour levels and process:
do lev=levbeg,levend
 !Integer index giving contour level:
indq=lev-nlevm+(lev-1)/nlevm-1

 !Counter for total number of grid line crossings:
ncr=0

 !Counter for total number of open contours originating in an edge:
npe=0

 !Contour level being sought:
qtmp=(dble(lev)-f12)*dq-qoff

 !Below, kib = grid box into which the contour (containing ncr) is going
 !       kob =   "   "  out of "    "     "         "       "    " coming
 !      [kob -> ncr -> kib:  ncr lies at the boundary between kob & kib]

 !   *** grid boxes are numbered 1 (lower left) to nxu*nyu (upper right) ***

!$omp parallel
 !Initialise number of crossings per box:
!$omp do
do k=1,nxny
  noctab(k)=0
enddo
!$omp enddo
 !-----------------------------------------------------------
 !Find x grid line crossings first (edge values are special):
 !Left edge:
ix=0
xgt=xgu(ix)

!$omp do
do iy=0,nyu
  qdy(iy)=qa(iy,ix)-qtmp
  isy(iy)=sign(one,qdy(iy))
enddo
!$omp enddo

!$omp do private(ncrloc,kob)
do iy=0,nyu-1
  if (isy(iy) .ne. isy(iy+1)) then
!$omp critical
    ncr=ncr+1
    ncrloc=ncr
!$omp end critical
    if (isy(iy) .lt. 0) then
       !A contour comes out of the boundary at this point:
      kib(ncrloc)=iy*nxu+1
!$omp critical
      npe=npe+1
      icre(npe)=ncrloc
!$omp end critical
    else
       !A contour goes into the boundary at this point:
      kib(ncrloc)=0
      kob=iy*nxu+1
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncrloc
    endif
    xcr(ncrloc)=xmin
    ycr(ncrloc)=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
  endif
enddo
!$omp enddo

 !Right edge:
ix=nxu
xgt=xgu(ix)

!$omp do
do iy=0,nyu
  qdy(iy)=qa(iy,ix)-qtmp
  isy(iy)=sign(one,qdy(iy))
enddo
!$omp enddo

!$omp do private(ncrloc,kob)
do iy=0,nyu-1
  if (isy(iy) .ne. isy(iy+1)) then
!$omp critical
    ncr=ncr+1
    ncrloc=ncr
!$omp end critical
    if (isy(iy) .gt. 0) then
       !A contour comes out of the boundary at this point:
      kib(ncrloc)=(iy+1)*nxu
!$omp critical
      npe=npe+1
      icre(npe)=ncrloc
!$omp end critical
    else
       !A contour goes into the boundary at this point:
      kib(ncrloc)=0
      kob=(iy+1)*nxu
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncrloc
    endif
    xcr(ncrloc)=xmax
    ycr(ncrloc)=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
  endif
enddo
!$omp enddo

 !Interior x grid lines:
!$omp do private(xgt,qdy,isy,ncrloc,inc,kaa,kob)
do ix=1,nxu-1
  xgt=xgu(ix)

  do iy=0,nyu
    qdy(iy)=qa(iy,ix)-qtmp
    isy(iy)=sign(one,qdy(iy))
  enddo

  do iy=0,nyu-1
    if (isy(iy) .ne. isy(iy+1)) then
!$omp critical
      ncr=ncr+1
      ncrloc=ncr
!$omp end critical
      inc=(1-isy(iy))/2
      kaa=iy*nxu+ix
      kib(ncrloc)=kaa+inc
      kob=kaa+1-inc
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncrloc
      xcr(ncrloc)=xgt
      ycr(ncrloc)=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
    endif
  enddo

enddo
!$omp enddo

 !----------------------------------------------------------
 !Find y grid line crossings next (edge values are special):
 !Bottom edge:
iy=0
ygt=ygu(iy)

!$omp do
do ix=0,nxu
  qdx(ix)=qa(iy,ix)-qtmp
  isx(ix)=sign(one,qdx(ix))
enddo
!$omp enddo

!$omp do private(ncrloc,kob)
do ix=0,nxu-1
  if (isx(ix) .ne. isx(ix+1)) then
!$omp critical
    ncr=ncr+1
    ncrloc=ncr
!$omp end critical
    if (isx(ix) .gt. 0) then
       !A contour comes out of the boundary at this point:
      kib(ncrloc)=ix+1
!$omp critical
      npe=npe+1
      icre(npe)=ncrloc
!$omp end critical
    else
       !A contour goes into the boundary at this point:
      kib(ncrloc)=0
      kob=ix+1
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncrloc
    endif
    ycr(ncrloc)=ymin
    xcr(ncrloc)=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
  endif
enddo
!$omp enddo

 !Top edge:
iy=nyu
ygt=ygu(iy)

!$omp do
do ix=0,nxu
  qdx(ix)=qa(iy,ix)-qtmp
  isx(ix)=sign(one,qdx(ix))
enddo
!$omp enddo

!$omp do private(ncrloc,kob)
do ix=0,nxu-1
  if (isx(ix) .ne. isx(ix+1)) then
!$omp critical
    ncr=ncr+1
    ncrloc=ncr
!$omp end critical
    if (isx(ix) .lt. 0) then
       !A contour comes out of the boundary at this point:
      kib(ncrloc)=koff+ix+1
!$omp critical
      npe=npe+1
      icre(npe)=ncrloc
!$omp end critical
    else
       !A contour goes into the boundary at this point:
      kib(ncrloc)=0
      kob=koff+ix+1
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncrloc
    endif
    ycr(ncrloc)=ymax
    xcr(ncrloc)=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
  endif
enddo
!$omp enddo
 !koff = nxu*(nyu-1) above

 !Interior y = constant grid lines:
!$omp do private(ygt,qdx,isx,ncrloc,inc,kaa,kob)
do iy=1,nyu-1
  ygt=ygu(iy)

  do ix=0,nxu
    qdx(ix)=qa(iy,ix)-qtmp
    isx(ix)=sign(one,qdx(ix))
  enddo

  do ix=0,nxu-1
    if (isx(ix) .ne. isx(ix+1)) then
!$omp critical
      ncr=ncr+1
      ncrloc=ncr
!$omp end critical
      inc=(1-isx(ix))/2
      kaa=(iy-1)*nxu+ix+1
      kib(ncrloc)=kaa+(1-inc)*nxu
      kob=kaa+inc*nxu
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncrloc
      ycr(ncrloc)=ygt
      xcr(ncrloc)=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
    endif
  enddo

enddo
!$omp enddo

 !----------------------------------------------------------------
 !Now re-build contours:
!$omp do
do icr=1,ncr
  free(icr)=.true.
enddo
!$omp enddo
!$omp end parallel

 !First deal with any open contours attached to boundaries:
if (npe .gt. 0) then
  do ie=1,npe
     !A new contour (indexed na) starts here:
    na=na+1
    inda(na)=indq
    ibeg=npta+1
    i1a(na)=ibeg

     !The starting node on the contour (coming out of a boundary):
    icr=icre(ie)

     !First point on the contour:
    npd=1
    xd(1)=xcr(icr)
    yd(1)=ycr(icr)

     !Find remaining points on the contour:
    k=kib(icr)
     !k is the box the contour is entering (0 if going into a boundary)
    do while (k .ne. 0)
      noc=noctab(k)
       !Use last crossing in this box (noc) as the next node:
      icrn=icrtab(k,noc)
       !icrn gives the next point after icr (icrn is leaving box k)
      noctab(k)=noc-1
       !noctab is usually zero now except for boxes with a
       !maximum possible 2 crossings
      npd=npd+1
       !Coordinates of new node:
      xd(npd)=xcr(icrn)
      yd(npd)=ycr(icrn)
      free(icrn)=.false.
      k=kib(icrn)
    enddo

     !Re-distribute nodes on this contour 3 times to reduce complexity:
    keep=.false.
    do
      call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
       !Delete contour if deemed too small (see renode_open):
      if (npa(na) .eq. 0) exit
      call renode_open(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
       !Delete contour if deemed too small (see renode_open):
      if (npd .eq. 0) exit
      call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
       !Delete contour if deemed too small (see renode_open):
      if (npa(na) .eq. 0) exit
       !Contour is big enough to keep:
      keep=.true.
      exit
    enddo
       
    if (keep) then
      npta=npta+npa(na)
      iend=ibeg+npa(na)-1
      i2a(na)=iend
      do i=ibeg,iend-1
        nextq(i)=i+1
      enddo
      nextq(iend)=0
    else
      na=na-1
    endif

    free(icr)=.false.
  enddo
endif

 !Next deal with remaining closed contours:
do icr=1,ncr
  if (free(icr)) then
     !A new contour (indexed na) starts here:
    na=na+1
    inda(na)=indq
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
    do while (icrn .ne. icr)
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
      call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
       !Delete contour if deemed too small (see renode_closed):
      if (npa(na) .eq. 0) exit
      call renode_closed(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
       !Delete contour if deemed too small (see renode_closed):
      if (npd .eq. 0) exit
      call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
       !Delete contour if deemed too small (see renode_closed):
      if (npa(na) .eq. 0) exit
       !Contour is big enough to keep:
      keep=.true.
      exit
    enddo

    if (keep) then 
      npta=npta+npa(na)
      iend=ibeg+npa(na)-1
      i2a(na)=iend
      do i=ibeg,iend-1
        nextq(i)=i+1
      enddo
      nextq(iend)=ibeg
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
end subroutine

!=======================================================================

subroutine con2ugrid(xq,yq,dq,qavg,nextq,nptq)
! Contour -> grid conversion.  The contours are represented by
! nodes (xq(i),yq(i)), i = 1, ..., nptq, where nextq(i) gives the
! index of the node following i, dq is the jump in q across all
! contours, and qavg is average value of the field.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: xq(npm),yq(npm)
integer:: nextq(npm)

 !Local parameters and arrays:
double precision:: qjx(0:nxu+1),qbot(0:nxu)
double precision:: dx(nptq),dy(nptq)
integer:: ixc(nptq),nxc(nptq)
logical:: crossx(nptq)

!----------------------------------------------------------------
!$omp parallel

 !Initialise interior x grid line crossing information and fill the
 !q jump array along lower boundary:
!$omp do
do i=1,nptq
  ixc(i)=int(one+dxxui*(xq(i)-xbeg))
enddo
!$omp enddo
 !Here xbeg is very slightly larger than xmin so that a point on
 !the left edge has ixc = 0, but one with xq just greater than xbeg
 !has ixc = 1.  Similarly, a point on the right edge has ixc = nxu+1.
 !Note: dxxui = dble(nxu)/((xmax-xmin)*(1-small)) where small = 1.d-12
 !(see casl.f90 initialisation)

!$omp do
do ix=0,nxu+1
  qjx(ix)=zero
enddo
!$omp enddo

!$omp do private(ia,py0,ix)
do i=1,nptq
  ia=nextq(i)
  if (ia .gt. 0) then
     !A node with ia = 0 terminates a contour at a boundary
    dx(i)=xq(ia)-xq(i)
    dy(i)=yq(ia)-yq(i)
    nxc(i)=ixc(ia)-ixc(i)
    crossx(i)=(nxc(i) .ne. 0)
    if ((yq(ia)-ybeg)*(ybeg-yq(i)) .gt. zero) then
       !The contour segment (i,ia) crosses y = ybeg; find x location:
      py0=(ybeg-yq(i))/dy(i)
      ix=int(one+dxxui*(xq(i)+py0*dx(i)-xbeg))
!$omp atomic
      qjx(ix)=qjx(ix)-dq*sign(one,dy(i))
       !Note: qjx gives the jump going from ix-1 to ix
    endif
  else
     !Here, there is no segment (i,next(i)) to consider:
    crossx(i)=.false.
  endif
enddo
!$omp enddo
 !Above, ybeg is very slightly greater than ymin to detect boundary crossings

!$omp single
 !Sum q jumps to obtain the gridded q along lower boundary:
qbot(0)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=1,nxu
  qbot(ix)=qbot(ix-1)+qjx(ix)
enddo
!$omp end single

!----------------------------------------------------------------
 !Initialise interior q jump array:
!$omp do
do ix=0,nxu
  do iy=0,nyu+1
    qa(iy,ix)=zero
  enddo
enddo
!$omp enddo

 !Determine x grid line crossings and accumulate q jumps:
!$omp do private(jump,ixbeg,sdq,ncr,ix,px0,iy)
do i=1,nptq
  if (crossx(i)) then
    jump=sign(1,nxc(i))
    ixbeg=ixc(i)+(jump-1)/2
    sdq=dq*sign(one,dx(i))
    ncr=0
    do while (ncr .ne. nxc(i)) 
      ix=ixbeg+ncr
      px0=(xxu(ix)-xq(i))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      iy=int(one+dyyui*(yq(i)+px0*dy(i)-ybeg))
       !Increment q jump between the grid lines iy-1 & iy:
!$omp atomic
      qa(iy,ix)=qa(iy,ix)+sdq
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo
!$omp enddo
!$omp end parallel

 !Get q values by sweeping through y:
do ix=0,nxu
  qa(0,ix)=qbot(ix)
  do iy=1,nyu
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

 !Restore average (use qjx as temp array):
do ix=0,nxu
  qjx(ix)=f12*(qa(0,ix)+qa(nyu,ix))
  do iy=1,nyu-1
    qjx(ix)=qjx(ix)+qa(iy,ix)
  enddo
enddo

qavg0=f12*(qjx(0)+qjx(nxu))
do ix=1,nxu-1
  qavg0=qavg0+qjx(ix)
enddo
qavg0=qavg0/dble(nxu*nyu)

qadd=qavg-qavg0
do ix=0,nxu
  do iy=0,nyu
    qa(iy,ix)=qa(iy,ix)+qadd
  enddo
enddo

return
end subroutine

!==========================================================================

 !Main end module
end module
