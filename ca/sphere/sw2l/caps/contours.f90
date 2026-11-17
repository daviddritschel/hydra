module contours

! Module contains all subroutines related to contour advection, for the 
! spe suite of f90 codes.

use constants

implicit none

 !Maximum number of PV levels:
integer,parameter:: nlevm=10000

 !Contour arrays:
double precision:: x(npm),y(npm),z(npm)
integer:: i1(nm),i2(nm),np(nm),ind(nm),lay(nm)
integer:: il1(nz),il2(nz),jl1(nz),jl2(nz)
integer:: next(0:npm),n,npt

 !Contour to grid conversion quantities (fine grid used in this module):
double precision:: clonf(ntf),slonf(ntf)
double precision,parameter:: dlf=twopi/dble(ntf),dlfi=dble(ntf)/(twopi+small)

 !Contour to grid conversion quantities (ultra-fine grid used in congen):
double precision:: clonu(ntu),slonu(ntu)
double precision,parameter:: dlu=twopi/dble(ntu),dlui=dble(ntu)/(twopi+small)

 !Weight arrays for grid to ultra-fine grid interpolation:
double precision:: w00(mgu,mgu),w10(mgu,mgu)
double precision:: w01(mgu,mgu),w11(mgu,mgu)
integer:: ixfw(ntu),iyfw(ngu)
integer:: ix0w(ntu),iy0w(ngu)
integer:: ix1w(ntu),iy1w(ngu)

 !Grid -> Contour arrays:
double precision:: xgu(ntu+1),ygu(0:ngu)
double precision:: qlev(2*nlevm),qoff
integer(kind=dbleint):: ibx(ntu,0:1)

 !Half-grid -> Full-grid tri-diagonal arrays:
double precision:: etd(nt),htd(nt),ptd(nt),xndeno

 !Coriolis frequency:
double precision:: fcor(ng)

!Constants used in contour node placement (Fontane & Dritschel JCP 2009):
double precision,parameter:: amu=0.2d0, ell=6.25d0*dl
double precision,parameter:: elf=one/ell**2, densf=one/(amu*sqrt(ell))
double precision,parameter:: dm=amu**2*ell/four, dm2=dm**2, dmsq=four*dm2
double precision,parameter:: dmi=two/dm, d4small=small*dl**4


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains 

!=======================================================================

subroutine init_contours
! Initialises all quantities needed for contour advection.

implicit none

double precision:: fac,pxu,pyu,pxc,pyc,rlonu,rlonf
integer(kind=dbleint):: ngulong
integer:: ixu,iyu,ix,iy,ixx,iyy,lev,i,j

!--------------------------------------------------------------------
 !Initialise area weights for interpolation of residual PV (qr) 
 !onto the ultra-fine horizontal grid:
fac=one/dble(mgu)
do ixu=1,mgu
  pxu=fac*dble(ixu-1)
  pxc=one-pxu
  do iyu=1,mgu
    pyu=fac*dble(iyu-1)
    pyc=one-pyu
    w00(iyu,ixu)=pyc*pxc
    w10(iyu,ixu)=pyu*pxc
    w01(iyu,ixu)=pyc*pxu
    w11(iyu,ixu)=pyu*pxu
  enddo
enddo

 !modulo values to access above weights:
do ix=1,ntu
  ixx=(ix-1)/mgu
  ix0w(ix)=1+ixx
  ix1w(ix)=2+ixx-nt*(ix0w(ix)/nt)
  ixfw(ix)=ix-mgu*ixx
enddo

do iy=1,ngu-1
  iyy=iy/mgu
  iy0w(iy)=iyy
  iy1w(iy)=1+iyy
  iyfw(iy)=1+iy-mgu*iyy
enddo

!--------------------------------------------------------------
 !Quantities needed in generating new contours (ugrid2con):
qoff=dq*dble(nlevm)
! qoff: should be a large integer multiple of the 
!       contour interval, dq.  The multiple should exceed 
!       the maximum expected number of contour levels.

do lev=1,2*nlevm
  qlev(lev)=(dble(lev)-f12)*dq-qoff
enddo

 !Fixed arrays:
do i=1,ntu
  rlonu=dlu*dble(i-1)-pi
  clonu(i)=cos(rlonu)
  slonu(i)=sin(rlonu)
enddo

 !Coordinates of grid lines (longitudes and latitudes):
do ix=1,ntu+1
  xgu(ix)=glxu*dble(ix-1)-pi
enddo
do iy=0,ngu
  ygu(iy)=glyu*dble(iy)-hpi
enddo

 !Grid box reference indices:
ngulong=ngu
do ix=1,ntu
  ibx(ix,1)=ngulong*(ix-1)
enddo
do ix=2,ntu
  ibx(ix,0)=ibx(ix-1,1)
enddo
ibx(1,0)=ibx(ntu,1)

!-----------------------------------------------
 !Fixed arrays:
do i=1,ntf
  rlonf=dlf*dble(i-1)-pi
  clonf(i)=cos(rlonf)
  slonf(i)=sin(rlonf)
enddo

 !Coriolis frequency (f):
do j=1,ng
  fcor(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

!----------------------------------------------------------------------
 !Tri-diagonal arrays for half grid -> full grid interpolation of qr
 !(used in congen.f90)

 !Initialise periodic tridiagonal problem:
htd(1)=one
ptd(1)=-f16*htd(1)
etd(1)=ptd(1)

do j=2,nt
  htd(j)=one/(one+f16*etd(j-1))
  ptd(j)=-f16*ptd(j-1)*htd(j)
  etd(j)=-f16*htd(j)
enddo

ptd(ntm1)=etd(ntm1)+ptd(ntm1)
do j=ntm2,1,-1
  ptd(j)=etd(j)*ptd(j+1)+ptd(j)
enddo

xndeno=one/(one-etd(nt)*ptd(1)-ptd(nt))

return
end subroutine init_contours

!==========================================================================

subroutine velint(u1,u2,v1,v2,u,v,w)
! Interpolates the gridded velocity field (uu,vv) in each layer to 
! the PV contours as (u,v,w) in Cartesian components.

implicit none

 !Passed variables:
double precision:: u1(ng,nt),u2(ng,nt)
double precision:: v1(ng,nt),v2(ng,nt)
double precision:: u(npt),v(npt),w(npt)

 !Local variables:
double precision:: uu(0:ng+1,nt,nz),vv(0:ng+1,nt,nz)
double precision:: ri,rj,abl,bbl,cbl,dbl
double precision:: ulon,ulat,rr,rri,zula
integer:: i,ic,ip1,j,jp1,k,iz

 !Extend the velocity field adjacent to poles (j = 1 and ng)
 !with a pi shift in longitude to simplify interpolation below:
uu(1:ng,:,1)=u1
uu(1:ng,:,2)=u2
vv(1:ng,:,1)=v1
vv(1:ng,:,2)=v2
do i=1,ng
  ic=i+ng
  uu(0,i,1)=-u1(1,ic)
  uu(0,i,2)=-u2(1,ic)
  vv(0,i,1)=-v1(1,ic)
  vv(0,i,2)=-v2(1,ic)

  uu(0,ic,1)=-u1(1,i)
  uu(0,ic,2)=-u2(1,i)
  vv(0,ic,1)=-v1(1,i)
  vv(0,ic,2)=-v2(1,i)

  uu(ngp1,i,1)=-u1(ng,ic)
  uu(ngp1,i,2)=-u2(ng,ic)
  vv(ngp1,i,1)=-v1(ng,ic)
  vv(ngp1,i,2)=-v2(ng,ic)

  uu(ngp1,ic,1)=-u1(ng,i)
  uu(ngp1,ic,2)=-u2(ng,i)
  vv(ngp1,ic,1)=-v1(ng,i)
  vv(ngp1,ic,2)=-v2(ng,i)
enddo

 !Loop over layers:
do iz=1,nz
   !Only process if the layer contains contours:
  if (jl2(iz) > 0) then
     !Bi-linearly zonal and meridional velocity at the contour nodes (x,y,z):
    do k=il1(iz),il2(iz)
      ri=dli*(pi+atan2(y(k),x(k)))
      i=1+int(ri)
      ip1=1+mod(i,nt)
      bbl=dble(i)-ri
      abl=one-bbl

      rj=dli*(hpidl+asin(z(k)))
      j=int(rj)
      jp1=j+1
      cbl=rj-dble(j)
      dbl=one-cbl

      ulon=bbl*(dbl*uu(j,i,iz)  +cbl*uu(jp1,i,iz))+  &
           abl*(dbl*uu(j,ip1,iz)+cbl*uu(jp1,ip1,iz))
      ulat=bbl*(dbl*vv(j,i,iz)  +cbl*vv(jp1,i,iz))+  &
           abl*(dbl*vv(j,ip1,iz)+cbl*vv(jp1,ip1,iz))

      rr=sqrt(x(k)**2+y(k)**2+small)
      rri=one/rr
      zula=-z(k)*ulat
       !Convert to Cartesians at the contour nodes:
      u(k)=(x(k)*zula-y(k)*ulon)*rri
      v(k)=(y(k)*zula+x(k)*ulon)*rri
      w(k)=rr*ulat
    enddo
  endif
enddo

return
end subroutine velint

!==========================================================================

subroutine con2grid(qc)
! Calculates the PV anomaly field (qc) from the PV contours
! (x,y,z) in each layer.  Takes away Coriolis frequency (fcor).

implicit none

 !Passed array:
double precision:: qc(ng,nt,nz)

 !Local variables:
double precision:: qa(0:ngf+1,ntf),qaend(ngf/2)
double precision:: cx(npt),cy(npt),cz(npt),sq(npt)
double precision:: sig,rlatc,p
integer:: ilm1(npt),ntc(npt)
integer:: iz,k,ka,jump,ncr,i,j,ioff,ip1
integer:: ic,io,ie,je,ngh,nth,nghh,nthh,nghp1

!----------------------------------------------------------------
 !Begin a major loop over layers
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,nz
  if (jl2(iz) > 0) then
     !There are contours in this layer.  Initialise crossing information:
    do k=il1(iz),il2(iz)
      ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
    enddo

    do k=il1(iz),il2(iz)
      ka=next(k)
      cx(k)=z(k)*y(ka)-y(k)*z(ka)
      cy(k)=x(k)*z(ka)-z(k)*x(ka)
      cz(k)=x(k)*y(ka)-y(k)*x(ka)
      ntc(k)=ilm1(ka)-ilm1(k)
    enddo

    do k=il1(iz),il2(iz)
      sig=sign(one,cz(k))
      sq(k)=dq*sig
      ntc(k)=ntc(k)-ntf*((2*ntc(k))/ntf)
      if (sig*dble(ntc(k)) < zero) ntc(k)=-ntc(k)
      if (abs(cz(k)) > zero) then
        cx(k)=cx(k)/cz(k)
        cy(k)=cy(k)/cz(k)
      endif
    enddo

     !Initialise PV jump array:
    qa=zero

     !Determine crossing indices:
    do k=il1(iz),il2(iz)
      if (ntc(k) .ne. 0) then
        jump=sign(1,ntc(k))
        ioff=ntf+ilm1(k)+(1+jump)/2
        ncr=0
        do while (ncr .ne. ntc(k))
          i=1+mod(ioff+ncr,ntf)
          rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
          j=int(rlatc)+1
          p=rlatc-dble(j-1)
          qa(j,i)=  qa(j,i)+(one-p)*sq(k)
          qa(j+1,i)=qa(j+1,i)+    p*sq(k)
          ncr=ncr+jump
        enddo
      endif
    enddo

     !Get PV values, at half latitudes, by sweeping through latitudes:
    do i=1,ntf
      do j=2,ngf
        qa(j,i)=qa(j,i)+qa(j-1,i)
      enddo
    enddo
     !Here, qa(j,i) stands for the PV at latitude j-1/2,
     !from j = 1, ..., ngf.

     !----------------------------------------------------------------------
     !Average PV values on the fine grid to get corresponding 
     !values on the inversion grid (ng,nt):
    ngh=ngf
    nth=ntf

    do while (ngh > ng)
       !Pre-store PV adjacent to poles at complementary longitudes (+pi):
      nthh=nth/2
      nghp1=ngh+1
      do i=1,nthh
        ic=i+nthh
        qa(0,i)=qa(1,ic)
        qa(0,ic)=qa(1,i)
        qa(nghp1,i)=qa(ngh,ic)
        qa(nghp1,ic)=qa(ngh,i)
      enddo

       !Work from SP to NP to define PV at full latitudes from averages
       !at adjacent half latitudes:
      do i=1,nth
        do j=0,ngh
          qa(j,i)=f12*(qa(j+1,i)+qa(j,i))
        enddo
      enddo

       !Now qa(j,i) is the PV at latitude j*(pi/ngh)-pi/2

       !Next 1-2-1 average these values to define PV at half latitudes
       !on a grid twice as coarse:
      nghh=ngh/2
      do i=1,nth
        do j=1,nghh
          je=2*j
          qa(j,i)=f12*qa(je-1,i)+f14*(qa(je-2,i)+qa(je,i))
        enddo
      enddo

       !Now perform analogous longitudinal 1-2-1 average:
      do j=1,nghh
        qaend(j)=f12*(qa(j,nth)+qa(j,1))
      enddo
      do i=1,nth-1
        ip1=i+1
        do j=1,nghh
          qa(j,i)=f12*(qa(j,i)+qa(j,ip1))
        enddo
      enddo
      do j=1,nghh
        qa(j,nth)=qaend(j)
      enddo
       !Now qa(j,i) gives the PV at the half-longitudes i + 1/2.

       !Average these on the twice coarser grid:
      do j=1,nghh
        qa(j,1)=f12*(qa(j,nth)+qa(j,1))
      enddo
      do i=2,nthh
        io=2*i-1
        ie=io-1
        do j=1,nghh
          qa(j,i)=f12*(qa(j,ie)+qa(j,io))
        enddo
      enddo

      ngh=nghh
      nth=nthh

    enddo

     !Finalise and take away f to define PV anomaly:
    do i=1,nt
      do j=1,ng
        qc(j,i,iz)=qa(j,i)-fcor(j)
      enddo
    enddo

  else
     !There are no contours in this layer; assign qc = 0:
    qc(:,:,iz)=zero
  endif
enddo
 !Ends loop over layers
    
return
end subroutine con2grid

!========================================================================
subroutine surgery
! Performs surgery and, afterwards, redistributes nodes

! Major revision 01/01/2001 by D. G. Dritschel to accelerate surgery
! using local boxes to minimise search costs.  All surgery is now
! done in a single way.

! Here we use rectangular, equal area z-longitude boxes, except at
! each pole where we use a circular patch of radius bwz, where bwz 
! is the z spacing of the grid lines.  This special treatment of
! the poles was introduced 07/07/2012 by dgd @ Wildwood Crest, NJ.

! Extended 7 November 2026 to multiple layers by DGD @ Kyoto
  
implicit none

 !Local arrays and parameters:
integer,parameter:: ngbs=nt*ng+1
! nlevm: max number of distinct vorticity levels
! ngbs:  (+1) max number of boxes used in surgery below
integer,parameter:: nsegm=npm
! nsegm: maximum number of segments belonging to contours of a 
!        given q level crossing through all the boxes.
!        The factor multiplying npm should be proportional to 
!        C*((dbar+amu*ell/2)/bw)^2, where C is equal to the 
!        ratio of the maximum number of nodes of a given q
!        level to npm, dbar is the average segment length,
!        and bw is the box width.

double precision:: dx(npt),dy(npt),dz(npt)
double precision:: xd(nprm),yd(nprm),zd(nprm)
double precision:: lonc(npt),dlonc(npt),dsq(npt)
double precision:: xa(npm),ya(npm),za(npm)
double precision:: xx,fnq,fnbl,fnbz,bw,bwl,bwli,fac,zbm,bwz,bwzi
double precision:: c0,adl,adz,delx,dely,delz,aa,cc
double precision:: dxa,dya,dza
integer:: jq1(nlevm),jq2(nlevm),iq1(nlevm),iq2(nlevm),levq(nlevm)
integer:: nspb(0:ngbs),kb1(0:ngbs),kb2(0:ngbs)
integer:: loc(nsegm),list(nsegm),node(nsegm)
integer:: i1a(nm),i2a(nm),nexta(npm)
integer:: i,ia,iz,j1,j2,levp,nlev,j,lev,levt,ibeg,iend
integer:: nptq,nbl,nbz,nbzp1,nblm1,nblp1,nbox,nseg,ib
integer:: mbl,mbl1,mbl2,mbz1,mbz2,mbz,mmbz,mb,k,ks,isb,is
integer:: jbeg,npd
 !Logicals:
logical:: avail(npt)

!------------------------------------------------------------------------
 !Get work arrays for efficient surgery:
do i=1,npt
  xx=atan2(y(i),x(i))
  lonc(i)=xx-twopi*dble(int(xx*hlxi))
enddo

do i=1,npt
  ia=next(i)
  xx=lonc(ia)-lonc(i)
  dlonc(i)=xx-twopi*dble(int(xx*hlxi))
   !dlonc: the change in longitude from node i to node ia
  dx(i)=x(ia)-x(i)
  dy(i)=y(ia)-y(i)
  dz(i)=z(ia)-z(i)
enddo

do i=1,npt
  dsq(i)=dx(i)**2+dy(i)**2+dz(i)**2
enddo

 !Zero out counters for contours and nodes for building new contours:
n=0
npt=0

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !              Begin a major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,nz
  if (jl2(iz) .gt. 0) then
     !There are contours in this layer.  Calculate beginning and 
     !ending contours (jq1,jq2) for each distinct value of ind:
    j1=jl1(iz)
    j2=jl2(iz)

    levp=ind(j1)
    nlev=1
    jq1(1)=j1
    levq(1)=levp
    do j=j1+1,j2
      lev=ind(j)
      if (lev .gt. levp) then
        jq2(nlev)=j-1
        nlev=nlev+1
        jq1(nlev)=j
        levq(nlev)=lev
        levp=lev
      endif
    enddo
    jq2(nlev)=j2
     !Note: levq(lev) gives the q level (ind) of contours.

    do lev=1,nlev
      iq1(lev)=i1(jq1(lev))
      iq2(lev)=i2(jq2(lev))
    enddo

     !Reset initial layer counters for contours and nodes:
    jl1(iz)=n+1
    il1(iz)=npt+1

     !-----------------------------------------------------------------
     !              Begin a major loop over q levels
     !-----------------------------------------------------------------

    do lev=1,nlev
       !Current contour level:       
      levt=levq(lev)
      ibeg=iq1(lev)
      iend=iq2(lev)

       !Work out optimal box number, nbox, for fast surgery:
      nptq=iend-ibeg+1
       !nptq: number of nodes having this q level.
      fnq=dble(nptq)

       !Balance boxing costs, (3+15*nseg/nbox)*nbox, with surgery 
       !search costs, 24*nptq*nseg/nbox, to work out optimal box number;
       !here we estimate nseg, the total number of segments counted
       !in all boxes, as (1+3*eps*nb)*nptq, where nb=sqrt(nbox);
       !this has been verified in realistically complex tests).
       !This balance results in a quartic equation for a = nb/sqrt(nptq)
       !whose solution only depends on r = 3*eps*sqrt(nptq), where
       !eps=(average distance between adjacent nodes)/(domain area).
       !Fortunately, a only varies from 1.129 to 1.265 over the entire
       !range of r (from 0 to infinity).  Here, therefore, we simply
       !take a = 1.2, i.e. nb = 1.2*sqrt(nptq), as an approximation:
      fnbl=1.2d0*sqrt(fnq*pi)
      fnbz=fnbl/pi
      nbl=max(min(nint(fnbl),nt),1)
       !nbl: number of boxes in longitude (lambda); nbl <= nt.
      nbz=max(min(nint(fnbz),ng),1)
       !nbz: number of boxes in sine(latitude) (z); nbz <= ng.
      nblm1=nbl-1
      nbzp1=nbz+1
      nbox=nbl*nbz+1
       !nbox: (+1) total number of boxes (including circular patches
       !      at each pole (box numbers 0 and nbox)
      bwl=twopi/dble(nbl)
      bwli=one/bwl
       !Ensure radius of polar "boxes" equals the z box division:
      fac=two/dble(nbz)
      zbm=one/sqrt(one+fac**2)
      bwz=zbm*fac
      bwzi=one/bwz

       !Box all segments [i,next(i)] to reduce search costs in surgery:
      do mb=0,nbox
        nspb(mb)=0
      enddo
      nseg=0

      c0=pi+bwl*dble(nbl+1)
      !above, pi is the domain half-width (in longitude).
      do ib=ibeg,iend
        i=next(ib)
         !find range of boxes spanned by segment (i,next(i)):
        adl=c0+lonc(i)
        if (dlonc(i) > zero) then
          mbl1=int(bwli*adl)-nbl
          mbl2=int(bwli*(adl+dlonc(i)))-nbl
        else
          mbl1=int(bwli*(adl+dlonc(i)))-nbl
          mbl2=int(bwli*adl)-nbl
        endif
         !mbl1,mbl2 is the longitude range spanned by the segment

        adz=zbm+z(i)
        if (dz(i) > zero) then
          mbz1=int(bwzi*adz+one)
          mbz2=int(bwzi*(adz+dz(i))+one)
        else
          mbz1=int(bwzi*(adz+dz(i))+one)
          mbz2=int(bwzi*adz+one)
        endif
         !mbl1,mbl2 is the z range spanned by the segment

        do mbz=mbz1,mbz2
          if (mbz == 0) then
             !This is the "box" at the south pole:
            nspb(0)=nspb(0)+1
            nseg=nseg+1
            loc(nseg)=0
            list(nseg)=ib
          else if (mbz == nbzp1) then
             !This is the "box" at the north pole:
            nspb(nbox)=nspb(nbox)+1
            nseg=nseg+1
            loc(nseg)=nbox
            list(nseg)=ib
          else
             !All other boxes:
            mmbz=mbz-1
            do mbl=mbl1,mbl2
              mb=nbl*mmbz+mod(mbl+nblm1,nbl)+1
              !nblm1=nbl-1 above
              nspb(mb)=nspb(mb)+1
              nseg=nseg+1
              loc(nseg)=mb
              list(nseg)=ib
            enddo
          endif
        enddo
      enddo
       !nspb counts number of segments that cross through box mb.
       !nseg counts the total number of segments crossing all boxes.
       !loc gives the box location of this segment
       !list gives the node before the one at the segment origin (i).

      kb1(0)=1
      do mb=0,nbox-1
        kb1(mb+1)=kb1(mb)+nspb(mb)
      enddo
      do mb=0,nbox
        kb2(mb)=kb1(mb)-1
      enddo

      do k=1,nseg
        mb=loc(k)
        ks=kb2(mb)+1
        node(ks)=list(k)
        kb2(mb)=ks
      enddo

       !With the above, segments [node(kb1(mb)),next(node(kb1(mb))],
       ![node(kb1(mb)+1),next(node(kb1(mb)+1)], ..., 
       ![node(kb2(mb)),next(node(kb2(mb))] cross through box mb.

       !Note: roughly (3+15*nseg/nbox)*nbox operations are required
       !      up to this point (in the loop over q levels) to do
       !      all the segment boxing.

       !Now search for surgery with segments crossing through the box 
       !containing i:
      do ib=ibeg,iend
        i=next(ib)

        mbz=int(bwzi*(zbm+z(i))+one)
        if (mbz == 0) then
           !This is the "box" at the south pole:
          mb=0
        else if (mbz == nbzp1) then
           !This is the "box" at the north pole:
          mb=nbox
        else
           !All other boxes:
          mb=nbl*(mbz-1)+int(bwli*(pi+lonc(i)))+1
           !mb: box containing the node i.
        endif

        do k=kb1(mb),kb2(mb)
          isb=node(k)
          is=next(isb)
           !Segment [is,next(is)] lies in box mb.  
           !Exclude the segments having node i at either endpoint:
          if ((is-ib)*(is-i) == 0) cycle
           !We will see if the node i gets within a distance dm 
           !to the line segment between is and next(is):
          delx=x(is)-x(i)
          dely=y(is)-y(i)
          delz=z(is)-z(i)
          aa=delx*dx(is)+dely*dy(is)+delz*dz(is)

           !Note: roughly 24*nptq*(nseg/nbox) operations are required 
           !up to the following statement (which is rarely satisfied); 
           !this is assumed to be the dominant cost of surgery.  
           !Even if each node i surgically reconnects on average 
           !once, the above estimate holds.

          if (aa*(aa+dsq(is))+d4small < zero) then
             !Passing this condition, a perpendicular can be dropped 
             !onto the line segment from is to next(is).  
             ![Tests show that this is satisfied only 6.4% of the time.]
             !Check distance to line segment:
            cc=(delx*dy(is)-dely*dx(is))**2+ &
               (dely*dz(is)-delz*dy(is))**2+ &
               (delz*dx(is)-delx*dz(is))**2
              !NB: sqrt(cc/dsq(is)) is the distance to the segment and 
              !    dm2=dm**2 below.  This distance must be strictly 
              !    positive, yet less than dm, to permit surgery:
            if (cc*(cc-dm2*dsq(is)) < zero) then
               !Surgery is now possible between node i and the segment 
               ![is,next(is)].
               ![Tests show that this is satisfied only 0.45% of the time.]
               !Move node i to a point midway between i and either is or 
               !next(is), whichever is closer:
              if (dsq(is)+two*aa < zero) then
                dxa=f12*(delx+dx(is))
                dya=f12*(dely+dy(is))
                dza=f12*(delz+dz(is))
                isb=is
                is=next(is)
              else
                dxa=f12*delx
                dya=f12*dely
                dza=f12*delz
              endif
               !Move node i & is to a common node; first deal with node i:
              x(i)=x(i)+dxa
              y(i)=y(i)+dya
              z(i)=z(i)+dza
              dx(i)=dx(i)-dxa
              dy(i)=dy(i)-dya
              dz(i)=dz(i)-dza
              dsq(i)=dx(i)**2+dy(i)**2+dz(i)**2
               !isb becomes the node before i:
              next(isb)=i
              dx(isb)=dx(isb)-dxa
              dy(isb)=dy(isb)-dya
              dz(isb)=dz(isb)-dza
              dsq(isb)=dx(isb)**2+dy(isb)**2+dz(isb)**2

               !Now deal with node is:
              x(is)=x(i)
              y(is)=y(i)
              z(is)=z(i)
              dx(is)=dx(is)+dxa
              dy(is)=dy(is)+dya
              dz(is)=dz(is)+dza
              dsq(is)=dx(is)**2+dy(is)**2+dz(is)**2
               !ib becomes the node before is:
              next(ib)=is
              dx(ib)=dx(ib)+dxa
              dy(ib)=dy(ib)+dya
              dz(ib)=dz(ib)+dza
              dsq(ib)=dx(ib)**2+dy(ib)**2+dz(ib)**2
                !Now update i and look for further possible surgery:
              i=next(ib)
              if (i == ib) exit
               !i = ib when a contour consists of a single node
            endif
          endif
        enddo
      enddo 

       !---------------------------------------------------------------
       !It remains to rebuild the contours using the next() information
       !and to eliminate small contours (having np < 5):
      do i=ibeg,iend
        avail(i)=.true.
      enddo

      jbeg=n+1
      do i=ibeg,iend
        if (avail(i)) then
           !avail(i) is true if node i has not yet been associated with
           !a contour. Start a new contour:
          n=n+1

           !First point on the contour:
          npd=1
          fac=one/sqrt(x(i)**2+y(i)**2+z(i)**2)
          xd(1)=fac*x(i)
          yd(1)=fac*y(i)
          zd(1)=fac*z(i)

          avail(i)=.false.
          is=next(i)
           !Generate the contour using next() and finish when the
           !original node (i) is reached:
          do while (is .ne. i) 
            npd=npd+1
            fac=one/sqrt(x(is)**2+y(is)**2+z(is)**2)
            xd(npd)=fac*x(is)
            yd(npd)=fac*y(is)
            zd(npd)=fac*z(is)
            avail(is)=.false.
            is=next(is)
          enddo

           !Re-distribute nodes on this contour:
          np(n)=npd
          if (npd > 3) &
            call renode(xd,yd,zd,npd,xa(npt+1),ya(npt+1),za(npt+1),np(n))
          if (np(n) > 3) then
            i1a(n)=npt+1
            npt=npt+np(n)
            i2a(n)=npt
            ind(n)=levt
            lay(n)=iz
            do is=i1a(n),npt-1
              nexta(is)=is+1
            enddo
            nexta(npt)=i1a(n)
          else
             !Delete contour if deemed too small (see renode):
            n=n-1
          endif
        endif
      enddo
       !Now go back and consider the next q level:
    enddo
     !--------------------------------------------------------------------
     !This ends the loop over q levels; surgery is complete for this layer
     !--------------------------------------------------------------------

     !Reset final layer counters for contours and nodes:
    jl2(iz)=n
    il2(iz)=npt
  endif
enddo
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !          The above ends the major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !Copy nodes back to those in the argument of the subroutine:
do i=1,npt
  x(i)=xa(i)
  y(i)=ya(i)
  z(i)=za(i)
  next(i)=nexta(i)
enddo
 !next(i) is not altered until surgery.

do j=1,n
  i1(j)=i1a(j)
  i2(j)=i2a(j)
enddo

return
end subroutine surgery

!=======================================================================

subroutine renode(xd,yd,zd,npd,xr,yr,zr,npr)
! Re-nodes a single closed contour (xd(i),yd(i),zd(i)), i = 1,...,npd 
!   and returns the new contour in (xr(i),yr(i),zr(i)), i = 1,...,npr

! Note: a closed contour closes on itself

! If npr = 0, the contour is too small and should be removed

implicit none

 !Passed variables:
double precision:: xd(nprm),yd(nprm),zd(nprm),xr(nprm),yr(nprm),zr(nprm)
integer:: npd,npr

 !Local variables:
double precision:: dx(npd),dy(npd),dz(npd),a(npd),b(npd),c(npd),d(npd),e(npd)
double precision:: dsa(npd),dsb(npd)
double precision:: bdif,ww,sum,fac,acc,p,eta,del,xtmp,ytmp,ztmp
double precision:: afac,ax,ay,az,sx,sy,sz
integer:: node(npd+1),nextr(0:npd)
integer:: i,ia,ib,ncorn,im,k,ibeg,iend,npseg
logical:: corner(npd)

!------------------------------------------------------------------
 !Define the node following a node:
do i=0,npd-1
  nextr(i)=i+1
enddo
nextr(npd)=1

 !Compute the contour increments:
do i=1,npd
  ia=nextr(i)
  dx(i)=xd(ia)-xd(i)
  dy(i)=yd(ia)-yd(i)
  dz(i)=zd(ia)-zd(i)
enddo

do i=1,npd
  dsa(i)=dx(i)**2+dy(i)**2+dz(i)**2
  e(i)=sqrt(dsa(i))
enddo

do ib=1,npd
  i=nextr(ib)
  dsb(i)=dsa(ib)
  a(i)=-dx(ib)
  b(i)=-dy(ib)
  c(i)=-dz(ib)
enddo

ncorn=0
do i=1,npd
  corner(i)=dx(i)*a(i)+dy(i)*b(i)+dz(i)*c(i) > zero
  if (corner(i)) then 
     !Keep track of corner locations for use in renoding below:
    ncorn=ncorn+1
    node(ncorn)=i
     !Set curvature to zero at corners:
    d(i)=zero
  else
    d(i)=(xd(i)*(dy(i)*c(i)-b(i)*dz(i))+ &
          yd(i)*(dz(i)*a(i)-c(i)*dx(i))+ &
          zd(i)*(dx(i)*b(i)-a(i)*dy(i)))/ &
          sqrt((a(i)*dsa(i)-dx(i)*dsb(i))**2+ &
               (b(i)*dsa(i)-dy(i)*dsb(i))**2+ &
               (c(i)*dsa(i)-dz(i)*dsb(i))**2+small3)
  endif
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npd
  ia=nextr(i)
  dsb(i)=e(i)*(d(ia)+d(i))
  bdif=e(i)*(d(ia)-d(i))
  a(i)=f16*bdif-f12*dsb(i)
  b(i)=f12*(dsb(i)-bdif)
  c(i)=f13*bdif
enddo

!------------------------------------------------------------------------
 !Use the spherical curvature expression (radius of the sphere = ell)
 !to ensure an adequate node density in low curvature regions.
do i=1,npd
  ww=one/(dsa(i)+dmsq)
  dsb(i)=ww*sqrt(elf*dsa(i)+dsb(i)**2)
  dsa(i)=ww*e(i)
enddo
 !NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; e(i)=sqrt{v(i)};
 !    dsb(i)/e(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

 !Re-assign curvature at a node from weighted average on either side
 !(v above is the weight):
do ib=1,npd
  i=nextr(ib)
  d(i)=(dsb(ib)+dsb(i))/(dsa(ib)+dsa(i))
enddo

 !Re-average to get interval value (effectively, four curvature
 !values go into getting the final interval value, dsb(i)):
do i=1,npd
  ia=nextr(i)
  dsb(i)=f12*(d(i)+d(ia))
enddo

 !Compute fractional number of nodes to be placed between old
 !nodes i and i+1:
do i=1,npd
  e(i)=e(i)*min(dmi,densf*sqrt(dsb(i))+dsb(i))
enddo
 !NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

sum=zero
do i=1,npd
  sum=sum+e(i)
enddo

 !Number of points on renoded contour:
npr=nint(sum)+1
if (npr < 3) then
   !Contour too small - remove it:
  npr=0
  return
endif

!-------------------------------------------------------------
 !Redistribute nodes making sure to preserve corner locations:
if (ncorn == 0) then
   !No corners - simplest case

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sum
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node is fixed; find the remaining node positions:
  xr(1)=xd(1)
  yr(1)=yd(1)
  zr(1)=zd(1)
  acc=zero
  i=0
  do im=2,npr
    do while (acc < one)
      i=i+1
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    del=f12*p*(one-p)
    xtmp=yd(i)*dz(i)-zd(i)*dy(i)
    ytmp=zd(i)*dx(i)-xd(i)*dz(i)
    ztmp=xd(i)*dy(i)-yd(i)*dx(i)
    afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
    ax=afac*xtmp
    ay=afac*ytmp
    az=afac*ztmp
    sx=dy(i)*az-dz(i)*ay
    sy=dz(i)*ax-dx(i)*az
    sz=dx(i)*ay-dy(i)*ax
    xr(im)=xd(i)+p*dx(i)+eta*ax+del*sx
    yr(im)=yd(i)+p*dy(i)+eta*ay+del*sy
    zr(im)=zd(i)+p*dz(i)+eta*az+del*sz
  enddo

else if (ncorn == 1) then
   !A single corner - start new contour at the corner:

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sum
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node (the corner) is fixed; find the remaining nodes:
  i=node(1)
  xr(1)=xd(i)
  yr(1)=yd(i)
  zr(1)=zd(i)
  acc=zero
  i=i-1
  do im=2,npr
    do while (acc < one)
      i=nextr(i)
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    del=f12*p*(one-p)
    xtmp=yd(i)*dz(i)-zd(i)*dy(i)
    ytmp=zd(i)*dx(i)-xd(i)*dz(i)
    ztmp=xd(i)*dy(i)-yd(i)*dx(i)
    afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
    ax=afac*xtmp
    ay=afac*ytmp
    az=afac*ztmp
    sx=dy(i)*az-dz(i)*ay
    sy=dz(i)*ax-dx(i)*az
    sz=dx(i)*ay-dy(i)*ax
    xr(im)=xd(i)+p*dx(i)+eta*ax+del*sx
    yr(im)=yd(i)+p*dy(i)+eta*ay+del*sy
    zr(im)=zd(i)+p*dz(i)+eta*az+del*sz
  enddo

else
   !Multiple corners - start new contour at the first corner and
   !preserve locations of all corners:
  npr=0
  node(ncorn+1)=node(1)
  ibeg=node(1)

  do k=1,ncorn
    iend=node(k+1)
    i=ibeg
    sum=zero
    do while (i .ne. iend)
      sum=sum+e(i)
      i=nextr(i)
    enddo

    if (sum > small) then
       !Adequate spacing exists between corners to renode this segment

       !Number of points on renoded contour segment:
      npseg=nint(sum)+1

       !Make the sum of e(i) equal to npseg
      fac=dble(npseg)/sum
      i=ibeg
      do while (i .ne. iend)
        e(i)=fac*e(i)
        i=nextr(i)
      enddo

       !Fix the first node along the segment as the corner:
      xr(npr+1)=xd(ibeg)
      yr(npr+1)=yd(ibeg)
      zr(npr+1)=zd(ibeg)

      if (npseg > 1) then
         !Find the remaining points along this segment:
        acc=zero
        i=ibeg-1
        do im=npr+2,npr+npseg
          do while (acc < one)
            i=nextr(i)
            acc=acc+e(i)
          enddo
          acc=acc-one
          p=one-acc/e(i)
          eta=p*(a(i)+p*(b(i)+p*c(i)))
          del=f12*p*(one-p)
          xtmp=yd(i)*dz(i)-zd(i)*dy(i)
          ytmp=zd(i)*dx(i)-xd(i)*dz(i)
          ztmp=xd(i)*dy(i)-yd(i)*dx(i)
          afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
          ax=afac*xtmp
          ay=afac*ytmp
          az=afac*ztmp
          sx=dy(i)*az-dz(i)*ay
          sy=dz(i)*ax-dx(i)*az
          sz=dx(i)*ay-dy(i)*ax
          xr(im)=xd(i)+p*dx(i)+eta*ax+del*sx
          yr(im)=yd(i)+p*dy(i)+eta*ay+del*sy
          zr(im)=zd(i)+p*dz(i)+eta*az+del*sz
        enddo
      endif

      npr=npr+npseg

    endif
    ibeg=iend

     !go on and consider the next segment (if any) between corners:
  enddo
endif

 !Normalise coordinates to unit sphere:
do i=1,npr
  fac=one/sqrt(xr(i)**2+yr(i)**2+zr(i)**2)
  xr(i)=fac*xr(i)
  yr(i)=fac*yr(i)
  zr(i)=fac*zr(i)
enddo

return
end subroutine renode

!==========================================================================

 !Main end module
end module contours
