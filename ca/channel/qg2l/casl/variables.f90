module variables

 !Include all variables which may change during the course of a simulation.

 !Time, time step, half time step and time when gridded data are saved:
double precision:: t,dt,hfdt,tgrid

 !Twist parameter for timing of surgery:
double precision:: twist

 !Counters for writing direct-access data (gridded and contours):
integer:: igrec,icrec

end module variables
