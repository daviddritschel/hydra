module congen

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! This module contains subroutines for rebuilding PV contours in a
! multi-layer periodic channel, given either the PV field and no
! contours, or the PV residual field together with contours.
! This creates new contours in either case.

! Revised by D G Dritschel on 31 July 2020 to use any number of layers.

! Open contours originating and terminating in a boundary added by
! D G Dritschel on 18 June 2012 @ Moscow
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Import contour advection module:
use contours

implicit none

 !Array for storing the PV field interpolated to the ultra-fine grid: 
double precision:: qa(0:nyup1,0:nxum1)

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

subroutine recontour(qq)
! Main routine for recontouring (from D & Ambaum, 1996, QJRMS)

! qq           : a gridded field added to that due to contours
! xq(i),yq(i)  : location of node i in the domain
! nextq(i)     : index of the node following node i 
!                *** this must be zero for an endpoint on the boundary ***
! layq(j)      : layer (integer) containing contour j
! indq(j)      : field level (integer) of contour j
! npq(j)       : number of nodes on contour j
! i1q(j)       : beginning node index on contour j
! i2q(j)       : ending node index on contour j
! qjump(iz)    : PV jump across all contours in layer iz
! qavg(iz)     : Average PV in layer iz
! jl1q(iz)     : beginning contour index in layer iz
! jl2q(iz)     : ending contour index in layer iz
! il1q(iz)     : beginning node index in layer iz
! il2q(iz)     : ending node index in layer iz
! nq           : number of contours
! nptq         : total number of nodes

implicit none

 !Passed array:
double precision:: qq(0:ny,0:nxm1,nz)

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

      !Bi-linear interpolate the residual qq to the fine grid and add to qa:
      do ix=0,nxum1
         ixf=ixfw(ix)
         ix0=ix0w(ix)
         ix1=ix1w(ix)

         do iy=0,nyu
            iyf=iyfw(iy)
            iy0=iy0w(iy)
            iy1=iy1w(iy)

            qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0,iz) &
                               +w10(iyf,ixf)*qq(iy1,ix0,iz) &
                               +w01(iyf,ixf)*qq(iy0,ix1,iz) &
                               +w11(iyf,ixf)*qq(iy1,ix1,iz)
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

         !Interpolate qq (which here contains the full field)
         !to the fine grid as qa:
         do ix=0,nxum1
            ixf=ixfw(ix)
            ix0=ix0w(ix)
            ix1=ix1w(ix)

            do iy=0,nyu
               iyf=iyfw(iy)
               iy0=iy0w(iy)
               iy1=iy1w(iy)

               qa(iy,ix)=w00(iyf,ixf)*qq(iy0,ix0,iz) &
                        +w10(iyf,ixf)*qq(iy1,ix0,iz) &
                        +w01(iyf,ixf)*qq(iy0,ix1,iz) &
                        +w11(iyf,ixf)*qq(iy1,ix1,iz)
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
 !Ends loop over layers

 !Copy arrays back to those in the argument of the subroutine:
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
! +/-qjump(iz)/2, +/-3*qjump(iz)/2, ....

implicit none

 !Passed index:
integer:: iz

 !Local parameters and arrays:
integer,parameter:: ncrm=3*nplm/4
 !ncrm:  max number of contour crossings of a single contour level
 !nplm:  max number of nodes in any contour level
 
integer(kind=dbleint),parameter:: &
     nxny=int(nxu,kind=dbleint)*int(nyu,kind=dbleint)
integer(kind=dbleint),parameter:: &
     koff=int(nxu,kind=dbleint)*int(nyu-1,kind=dbleint)
integer(kind=dbleint):: k,kob,kaa,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(0:nxu),qdy(0:nyu)
double precision:: xd(nprm),yd(nprm)
double precision:: qoff,dq,dqi,qtmp,xgt,ygt,xx
integer:: isx(0:nxu),isy(0:nyu),icre(nm),icrtab(nxny,2)
integer:: levbeg,levend,lev,levt,noc,icrn
integer:: ncr,i,ix,iy,icr,inc,ibeg,iend,npd,npe,ie
integer(kind=halfint):: noctab(nxny)
logical:: free(ncrm),keep

 !----------------------------------------
 !Initialise various constants used below:
dq=qjump(iz)
dqi=one/dq
qoff=dq*dble(nlevm)
 !qoff: should be a large integer multiple of the contour interval, dq.  
 !The multiple should exceed the maximum expected number of contour levels.

 !--------------------------------------------------------
 !First get the beginning and ending contour levels:
levbeg=int((qoff+minval(qa(0:nyu,:)))*dqi+f12)+1
levend=int((qoff+maxval(qa(0:nyu,:)))*dqi+f12)

if (levbeg <= levend) then
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !Loop over contour levels and process:
   do lev=levbeg,levend
      !Integer index giving contour level:
      levt=lev-nlevm+(lev-1)/nlevm-1

      !Counter for total number of grid line crossings:
      ncr=0

      !Counter for total number of open contours originating in an edge:
      npe=0

      !Contour level being sought:
      qtmp=(dble(lev)-f12)*dq-qoff

      !Below, kib = grid box into which the contour (containing ncr) is going
      !       kob =   "   "  out of "    "     "         "       "    " coming
      !      [kob -> ncr -> kib:  ncr lies at the boundary between kob & kib]

      !   *** grid boxes are numbered 1 (lower left) to nxu*nyu (upper right).

      !Initialise number of crossings per box:
      noctab=0

      !-----------------------------------------------------------
      !Find x grid line crossings first:
      do ix=0,nxum1
         xgt=xgu(ix)

         do iy=0,nyu
            qdy(iy)=qa(iy,ix)-qtmp
            isy(iy)=sign(one,qdy(iy))
         enddo

         do iy=0,nyu-1
            if (isy(iy) /= isy(iy+1)) then
               ncr=ncr+1
               inc=(1-isy(iy))/2
               kaa=iy*nxu
               kib(ncr)=kaa+ibx(ix+inc)
               kob=kaa+ibx(ix+1-inc)
               noctab(kob)=noctab(kob)+1
               icrtab(kob,noctab(kob))=ncr
               xcr(ncr)=xgt
               ycr(ncr)=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
            endif
         enddo

      enddo

      !----------------------------------------------------------
      !Find y grid line crossings next (edge values are special):
      !Bottom edge:
      iy=0
      ygt=ygu(iy)

      do ix=0,nxum1
         qdx(ix)=qa(iy,ix)-qtmp
         isx(ix)=sign(one,qdx(ix))
      enddo
      qdx(nxu)=qdx(0)
      isx(nxu)=isx(0)

      do ix=0,nxum1
         if (isx(ix) /= isx(ix+1)) then
            ncr=ncr+1
            if (isx(ix) > 0) then
               !A contour comes out of the boundary at this point:
               kib(ncr)=ix+1
               npe=npe+1
               icre(npe)=ncr
            else
               !A contour goes into the boundary at this point:
               kib(ncr)=0
               kob=ix+1
               noctab(kob)=noctab(kob)+1
               icrtab(kob,noctab(kob))=ncr
            endif
            ycr(ncr)=ymin
            xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
            xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
         endif
      enddo

      !Top edge:
      iy=nyu
      ygt=ygu(iy)

      do ix=0,nxum1
         qdx(ix)=qa(iy,ix)-qtmp
         isx(ix)=sign(one,qdx(ix))
      enddo
      qdx(nxu)=qdx(0)
      isx(nxu)=isx(0)

      do ix=0,nxum1
         if (isx(ix) /= isx(ix+1)) then
            ncr=ncr+1
            if (isx(ix) < 0) then
               !A contour comes out of the boundary at this point:
               kib(ncr)=koff+ix+1
               npe=npe+1
               icre(npe)=ncr
            else
               !A contour goes into the boundary at this point:
               kib(ncr)=0
               kob=koff+ix+1
               noctab(kob)=noctab(kob)+1
               icrtab(kob,noctab(kob))=ncr
            endif
            ycr(ncr)=ymax
            xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
            xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
         endif
      enddo
      !koff = nxu*(nyu-1) above

      !Interior y = constant grid lines:
      do iy=1,nyu-1
         ygt=ygu(iy)

         do ix=0,nxum1
            qdx(ix)=qa(iy,ix)-qtmp
            isx(ix)=sign(one,qdx(ix))
         enddo
         qdx(nxu)=qdx(0)
         isx(nxu)=isx(0)

         do ix=0,nxum1
            if (isx(ix) /= isx(ix+1)) then
               ncr=ncr+1
               inc=(1-isx(ix))/2
               kaa=(iy-1)*nxu+ix+1
               kib(ncr)=kaa+(1-inc)*nxu
               kob=kaa+inc*nxu
               noctab(kob)=noctab(kob)+1
               icrtab(kob,noctab(kob))=ncr
               ycr(ncr)=ygt
               xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
               xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
            endif
         enddo

      enddo

      !----------------------------------------------------------------
      !Now re-build contours:
      do icr=1,ncr
         free(icr)=.true.
      enddo

      !First deal with any open contours attached to boundaries:
      if (npe > 0) then
         do ie=1,npe
            !A new contour (indexed na) starts here:
            na=na+1
            inda(na)=levt
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
            do while (k /= 0)
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
               if (npa(na) == 0) exit
               call renode_open(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
               !Delete contour if deemed too small (see renode_open):
               if (npd == 0) exit
               call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
               !Delete contour if deemed too small (see renode_open):
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
               nexta(iend)=0
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
               call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
               !Delete contour if deemed too small (see renode_closed):
               if (npa(na) == 0) exit
               call renode_closed(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
               !Delete contour if deemed too small (see renode_closed):
               if (npd == 0) exit
               call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
               !Delete contour if deemed too small (see renode_closed):
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
! Contour -> grid conversion in layer iz.  The contours are represented
! by nodes (xq(i),yq(i)), i = 1, ..., nptq, where nextq(i) gives the
! index of the node following i, dq is the jump in q across all contours,
! and qavg(iz) is average value of the field.

implicit none

 !Passed index:
integer:: iz

 !Local parameters and arrays:
double precision:: qjx(0:nxum1),qbot(0:nxum1)
double precision:: dx(nptq),dy(nptq)
double precision:: dq,sdq,px0,py0,qavg0,qadd,xx
integer:: ixc(nptq),nxc(nptq)
integer:: i,ia,ix,ixdif,ixbeg,iy,jump,ncr
logical:: crossx(nptq)

!------------------------------------------------------------------
! PV jump in this layer:
dq=qjump(iz)

! Initialise interior x grid line crossing information and fill the
! q jump array along lower boundary:
do i=il1q(iz),il2q(iz)
   ixc(i)=1+int(glxui*(xq(i)-xmin))
enddo

do ix=0,nxum1
   qjx(ix)=zero
enddo

dq=qjump(iz)

do i=il1q(iz),il2q(iz)
   ia=nextq(i)
   if (ia > 0) then
     !A node with ia = 0 terminates a contour at a boundary
      xx=xq(ia)-xq(i)
      dx(i)=xx-ellx*dble(int(xx*hlxi))
      dy(i)=yq(ia)-yq(i)
      ixdif=ixc(ia)-ixc(i)
      nxc(i)=ixdif-nxu*((2*ixdif)/nxu)
      crossx(i)=(nxc(i) /= 0)
      if ((yq(ia)-ybeg)*(ybeg-yq(i)) > zero) then
         !The contour segment (i,ia) crosses y = ybeg; find x location:
         py0=(ybeg-yq(i))/dy(i)
         xx=xq(i)+py0*dx(i)
         xx=oms*(xx-ellx*dble(int(xx*hlxi)))
         ix=int(glxui*(xx-xmin))
         qjx(ix)=qjx(ix)-dq*sign(one,dy(i))
         !Note: qjx gives the jump going from ix to ix+1
      endif
   else
      !Here, there is no segment (i,next(i)) to consider:
      crossx(i)=.false.
   endif
enddo
 !Above, ybeg is very slightly greater than ymin to detect boundary crossings

 !Sum q jumps to obtain the gridded q along lower boundary:
qbot(0)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=0,nxum2
   qbot(ix+1)=qbot(ix)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
do ix=0,nxum1
   do iy=0,nyup1
      qa(iy,ix)=zero
   enddo
enddo

 !Determine x grid line crossings and accumulate q jumps:
do i=il1q(iz),il2q(iz)
   if (crossx(i)) then
      jump=sign(1,nxc(i))
      ixbeg=ixc(i)+(jump-1)/2+nxu
      sdq=dq*sign(one,dx(i))
      ncr=0
      do while (ncr /= nxc(i))
         ix=mod(ixbeg+ncr,nxu)
         xx=xgu(ix)-xq(i)
         px0=(xx-ellx*dble(int(xx*hlxi)))/dx(i)
         !The contour crossed the fine grid line ix at the point
         !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
         iy=int(one+dyyui*(yq(i)+px0*dy(i)-ybeg))
         !Increment q jump between the grid lines iy-1 & iy:
         qa(iy,ix)=qa(iy,ix)+sdq
         !Go on to consider next x grid line (if there is one):
         ncr=ncr+jump
      enddo
   endif
enddo

 !Get q values by sweeping through y:
do ix=0,nxum1
   qa(0,ix)=qbot(ix)
   do iy=1,nyu
      qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
   enddo
enddo

 !Restore average (use qjx as temp array):
do ix=0,nxum1
   qjx(ix)=f12*(qa(0,ix)+qa(nyu,ix))
   do iy=1,nyu-1
      qjx(ix)=qjx(ix)+qa(iy,ix)
   enddo
enddo

qavg0=zero
do ix=0,nxum1
   qavg0=qavg0+qjx(ix)
enddo
qavg0=qavg0/dble(nxu*nyu)

qadd=qavg(iz)-qavg0
qa=qa+qadd

return
end subroutine con2ugrid

!==========================================================================

 !Main end module
end module congen
