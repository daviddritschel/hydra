module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: nz=2, ng=N_G, nt=2*ng
integer,parameter:: ncont=N_CONTQ
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! nz     : number of layers (here 2, but the contour modules are general)
! ng     : number of  latitudinal grid points (using a half grid)
! nt     : number of longitudinal grid points (always equal to 2*ng)
! ncont  : number of PV jumps used to represent f = 2*Omega*sin(lat)
! dt     : time step (*fixed*)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: vm11=V_11, vm12=V_12, vm21=V_21, vm22=V_22
double precision,parameter:: hbar1=H_BAR1, hbar2=H_BAR2, denrat=DEN_RAT
double precision,parameter:: gravity=ACC_GRAV, ctsq1=C_MODE1, ctsq2=C_MODE2
double precision,parameter:: rth=R_THERM, drate=D_RATE, nnu=POW_HYPER
double precision,parameter:: fqrms=FQ_RMS, tfq=T_FQ
integer,parameter:: nbeg=N_BEG, nend=N_END
! vmij       : Coefficients to convert layer quantities to mode quantities;
!              for a field F_j in layers j = 1 & 2, its projection onto
!              vertical mode 1 is F_1*vm11 + F_2*vm12, while for mode 2,
!              its projection is  F_1*vm21 + F_2*vm22
! hbar1      : Mean fluid depth of layer 1, H_1
! hbar2      : Mean fluid depth of layer 2, H_2
!              *** Note, wlog we take hbar1 + hbar2 = 1 throughout ***
! denrat     : density ratio, rho_2/rho_1
!              *** Note, denrat must be strictly less than 1 ***
! gravity    : The acceleration due to gravity (see hinvert script)
! ctsq1      : Short-scale gravity-wave speed associated with 1st vertical mode
! ctsq2      : Short-scale gravity-wave speed associated with 2nd vertical mode
! rth        : thermal damping rate
! drate      : drate*fpole is the hyperviscous damping rate on wavenumber ng
! nnu        : power of hyperviscosity
! fqrms      : r.m.s. PV forcing amplitude
! tfq        : Markovian de-correlation time of the forcing
! nbeg, nend : range in order n of spherical harmonics Y_n^m used

!----------------------------------------------------------------------------
! For implementing topography if desired:
logical,parameter:: topogr=TOP_FLAG
 !Random number seed:
integer,parameter:: iseed=I_SEED

end module parameters
