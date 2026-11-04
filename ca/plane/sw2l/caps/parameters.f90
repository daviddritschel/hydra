module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!----------------------------------------------------------------------------
! ==> Numerical parameters <==

integer,parameter:: nz=2, ng=N_G
double precision,parameter:: ncont=N_CONTQ
 !Simulation time length etc..
double precision,parameter:: dt=T_STEP, tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE, tcsave=T_CSAVE
! nz     : number of layers (here 2, but the contour modules are general)
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box.)
! ncont  : Number of PV jumps (contours) spanning q_min to q_max
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

!----------------------------------------------------------------------------
! ==> Physical parameters <==

! *** Note, the Coriolis frequency f (cof in the code) is 4*pi by default ***
double precision,parameter:: vm11=V_11, vm12=V_12, vm21=V_21, vm22=V_22
double precision,parameter:: hbar1=H_BAR1, hbar2=H_BAR2, alpha=DEN_RAT
double precision,parameter:: gravity=ACC_GRAV
double precision,parameter:: cdamp=C_DAMP, nnu=POW_HYPER
! vmij   : Coefficients to convert layer quantities to mode quantities;
!          for a field F_j in layers j = 1 & 2, its projection onto
!          vertical mode 1 is F_1*vm11 + F_2*vm12, while for mode 2,
!          its projection is  F_1*vm21 + F_2*vm22
! hbar1  : Mean fluid depth of layer 1, H_1
! hbar2  : Mean fluid depth of layer 2, H_2
!          *** Note, wlog we take hbar1 + hbar2 = 1 throughout ***
! alpha  : density ratio, rho_2/rho_1
!          *** Note, alpha must be strictly less than 1 ***
! gravity: The acceleration due to gravity (see hinvert script)
! cdamp  : This times f is the damping rate on wavenumber ng/2
! nnu    : The power of hyperviscosity (applied to residual PV,
!          divergence and acceleration divergence).

!----------------------------------------------------------------------------
! For implementing topography if desired:
logical,parameter:: topogr=TOP_FLAG
 !Random number seed:
integer,parameter:: iseed=I_SEED

end module parameters
