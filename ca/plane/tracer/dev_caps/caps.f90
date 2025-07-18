!#########################################################################
!                    The Doubly-Periodic Single-Layer 
!        Quasi-Geostrophic Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!        Code adapted from the code in imhd on 17 Jan 2014 @ St Andrews

!          This code simulates the following system of equations:

!             Dq/Dt = S 

!          where q is the QG PV. The PV source term is

!             S = -r*zeta + (psi-psi_eq)/(tau*L_D^2) + F

!          The velocity field (u,v) is found by inverting 

!             Lap{psi} - psi/L_D^2 = q - beta*y   
!             u = -dpsi/dy ; v = dpsi/dx                   

!          where
!             r      is the Ekman damping rate
!             tau    is the thermal damping rate
!             L_D    is the Rossby deformation length
!             psi_eq is the thermal equilibrium streamfunction
!             F      is stochastic forcing (by vortex injection)
!             beta   is the planetary vorticity gradient

!          We split the PV evolution equation into *three* equations,

!             Dq_c/Dt = 0  :  contour advection
!             Dq_s/Dt = 0  :  pseudo-spectral
!             Dq_d/Dt = S  :  pseudo-spectral

!          such that q is a weighted sum of these,

!             q = F(q_s) + (1-F)*(q_c) + q_d

!          where F is a low-pass filter defined in spectral.f90 and 1-F is
!          a complementary high pass filter (see Dritschel & Fontane, JCP,
!          2010).

!          Hence advection at large to intermediate scales is controlled 
!          by the pseudo-spectral method, whereas advection at intermediate
!          to small scales is controlled by contour advection (where it 
!          is most accurate).  The source term S is handled entirely by
!          the pseudo-spectral method.

!          At the beginning of each time step, q_s is replaced by q,
!          while q_d is replaced by (1-F)(q-q_c).  q_c remains as
!          contours for a period of time determined by twistmax below.
!          After this period, q is obtained on an ultra-fine grid and
!          re-contoured so that the accumulated forcing in q_d is
!          given to the contours in q_c (to the extent possible).

!     The full algorithm consists of the following modules:
!        casl.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        generic.f90   : Generic service routines for CASL;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded 
!                        fields using a PS method along with contours.
!----------------------------------------------------------------------------
program casl

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Obtain new PV contours:
  call recont

   !Advect PV and other fields until next recontouring or end:
  call evolve

enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local variables:
double precision:: ff(ny,nx)
integer:: ix,iy,kx,ky

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!-----------------------------------------------------------------
 !Read in full PV, possibly subtract beta*y, and convert to spectral space:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)

if (beffect) then
   !Subtract beta*y to define the PV anomaly (use ff below for FFT):
  do ix=1,nx
    do iy=1,ny
      qr(iy,ix)=qr(iy,ix)-bety(iy)
    enddo
  enddo
   !Choose an integral number of PV jumps (required in recontour):
  qjump=beta*elly/dble(ncontq)
  write(*,*)
  write(*,'(a,f14.8,a,f12.8)') 'beta = ',beta,'   qjump = ',qjump
else
   !Choose contour interval based on range of PV values:
  call contint(qr,ncontq,qjump)
  write(*,*)
  write(*,'(a,1x,f13.8)') ' qjump = ',qjump
endif

 !Copy qr into ff before taking Fourier transform (ff is overwritten):
do ix=1,nx
  do iy=1,ny
    ff(iy,ix)=qr(iy,ix)
  enddo
enddo
 !Convert PV anomaly (in ff) to spectral space as qs:
call ptospc(nx,ny,ff,qs,xfactors,yfactors,xtrig,ytrig)

 !Read in thermal equilibrium streamfunction (if present) and convert:
if (heating) then
  open(12,file='psieq.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) dum,ff
  close(12)
  call ptospc(nx,ny,ff,ppeq,xfactors,yfactors,xtrig,ytrig)
endif

 !Initialise random number generator:
do i=1,iseed
  dum=rand(0)
enddo

!------------------------------------------------------------
 !Initially there are no contours:
nq=0
nptq=0

 !Initialise stochastic vortex forcing variables if used:
if (stoch) then
  if (ivor .eq. 1) then
     !Point vortices of mean (grid) vorticity of +/-vorvor are
     !added at an average rate of dnvor vortices per unit time:
    dnvor=two*esr*(three*pi/vorvor)**2/garea
  else
     !Point vortex dipoles (concentrated to points) are added
     !at an average rate of dnvor dipoles per unit time, with 
     !a maximum absolute grid vorticity of vorvor:
    dnvor=six*esr*(pi/vorvor)**2/garea
  endif
!Above, esr is the enstrophy input rate.

   !Initialize random # generator on first call:
  do i=1,iseed
    uni=rand(0)
  enddo

   !Initialise total added vortices:
  totnvor=zero
endif

 !Initialise time step so that subroutine adapt chooses a suitable one:
dt=zero

 !Set final time for simulation end:
itime=int((t+small)/tgsave)
tgrid=tgsave*dble(itime)
tfin=tgrid+tsim

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='complexity.asc',status='unknown')
open(15,file='ene.asc',status='unknown')
open(16,file='norms.asc',status='unknown')
open(17,file='monitor.asc',status='unknown')

 !Open file for 1d PV & current density spectra:
open(51,file='spectra.asc',status='unknown')

 !Open files for coarse grid saves:
open(31,file='qq.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

 !Open files for contour writes:
open(80,file='cont/qqsynopsis.asc',status='unknown')
open(83,file='cont/qqresi.r4',form='unformatted',access='direct', &
                          & status='replace',recl=nbytes)

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect PV until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour(qr)
write(*,'(a,i8,a,i9,a,f9.5)') '   nq = ',nq,'   nptq = ',nptq,'   dq = ',qjump

return 
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(16)
close(17)
close(31)
close(51)
close(80)
close(83)

return
end subroutine

 !End main program
end program
!=======================================================================
