#!/usr/bin/env python3

# This script plots the fields for a selected field, and either the
# full field, the balanced field, or the imbalanced field.  Compares
# results in 4 directories, indicated below (2 in SW and 2 in GN),
# at a selected resolution.

#========== Perform the generic imports =========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# set tick label size:
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 2

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    mpow=0
    rmult=fmax-fmin
    while rmult < 10.0:
       mpow+=1
       rmult=rmult*10.0

    while rmult >= 100.0:
       mpow-=1
       rmult=rmult/10.0

    emag=10.0**(float(-mpow))

    kmult=int(rmult/10.0)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2.0*emag
    elif kmult < 4:
       ci=4.0*emag
    elif kmult < 8:
       ci=10.0*emag
    else:
       ci=20.0*emag

    return ci

#=================================================================
# Select data to compare:
field_list=['h','zeta','delta','gamma','gamma-tilde']
field_acro=['hh','zz','dd','gg','gt']

print
print ' Select the field type from the following options:'
print
print ' (1) full;'
print ' (2) balanced;'
print ' (3) imbalanced.'
print
option=int(raw_input(' Option (default 1)? ') or 1)

print
t=float(raw_input(' Time to show (default 25)? ') or 25.0)

print
ng=int(raw_input(' Resolution (default 256)? ') or 256)
print

# Define grid:
#xg=np.linspace(-np.pi,np.pi,ng+1)
#yg=xg
#X,Y=np.meshgrid(xg,yg)
N=ng*ng

#=================================================================
# List of directories to compare (need the final /):
dirend=str(ng)+'/'
direc=['sw/ng'+dirend,'sw/bal_ng'+dirend,'gn/ng'+dirend,'gn/bal_ng'+dirend]
# Corresponding labels on the subplots:
label_list=['SW','SW-bal','GN','GN-bal']

# Open ene.asc file in one directory to get time between frames:
in_file=open(direc[0]+'ene.asc','r')
time, etot = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]
# Frame corresponding to time chosen:
frame=int((t+0.0001)/dt)

#=================================================================
# Loop over fields (h, zeta, delta and gamma):
for k in range(5):
   field=field_list[k]
   acron=field_acro[k]

   print ' ================================='
   print '  *** Processing data for ',field
   print ' ================================='

   # Read data into arrays for plotting:
   Z1=np.empty([ng+1,ng+1])
   Z2=np.empty([ng+1,ng+1])
   Z3=np.empty([ng+1,ng+1])
   Z4=np.empty([ng+1,ng+1])
   if option==1:
      # Read full fields in all directories:
      in_file=open(direc[0]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[1]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[2]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[3]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z4[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      # Create unique plot output file name:
      outfile=field+'_n'+str(ng)+'_t'+str(int(t+0.01))+'.eps'

      # Create title of entire plot:
      field_labl=['$\\tilde{h}$','$\zeta$','$\delta$','$\gamma$']
      title_label=field_labl[k]

   elif option==2:
      # Read balanced fields in all directories and add periodic edges:
      in_file=open(direc[0]+'b'+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[1]+'b'+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[2]+'b'+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[3]+'b'+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z4[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      # Create unique plot output file name:
      outfile=field+'_bal_n'+str(ng)+'_t'+str(int(t+0.01))+'.eps'

      # Create title of entire plot:
      field_labl=['$\\tilde{h}_b$','$\zeta_b$','$\delta_b$','$\gamma_b$']
      title_label=field_labl[k]

   else:
      # Read full and balanced fields in all directories to create imbalanced
      # fields and add periodic edges:
      in_file=open(direc[0]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      in_file=open(direc[0]+'b'+acron+'.r4','r')
      raw_arrayb=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T-raw_arrayb[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[1]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      in_file=open(direc[1]+'b'+acron+'.r4','r')
      raw_arrayb=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T-raw_arrayb[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[2]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      in_file=open(direc[2]+'b'+acron+'.r4','r')
      raw_arrayb=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T-raw_arrayb[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      in_file=open(direc[3]+acron+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      in_file=open(direc[3]+'b'+acron+'.r4','r')
      raw_arrayb=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      Z4[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T-raw_arrayb[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

      # Create unique plot output file name:
      outfile=field+'_imb_n'+str(ng)+'_t'+str(int(t+0.01))+'.eps'

      # Create title of entire plot:
      field_labl=['$\\tilde{h}_i$','$\zeta_i$','$\delta_i$','$\gamma_i$']
      title_label=field_labl[k]

#====================================================================
   # Add periodic edges:
   Z1[ng,0:ng]=Z1[0,0:ng]
   Z1[0:ng+1,ng]=Z1[0:ng+1,0]
   Z2[ng,0:ng]=Z2[0,0:ng]
   Z2[0:ng+1,ng]=Z2[0:ng+1,0]
   Z3[ng,0:ng]=Z3[0,0:ng]
   Z3[0:ng+1,ng]=Z3[0:ng+1,0]
   Z4[ng,0:ng]=Z4[0,0:ng]
   Z4[0:ng+1,ng]=Z4[0:ng+1,0]

#====================================================================
   # Work out the overall min/max values:
   zmin1=np.amin(Z1)
   zmin2=np.amin(Z2)
   zmin3=np.amin(Z3)
   zmin4=np.amin(Z4)
   zmax1=np.amax(Z1)
   zmax2=np.amax(Z2)
   zmax3=np.amax(Z3)
   zmax4=np.amax(Z4)

   print
   print ' Minimum and maximum field values for each simulation:'
   print
   print ' Simulation    Min field value      Max field value'
   print ' ----------    ---------------      ---------------'
   print '   SW         ',"{:14.10f}".format(zmin1),'     ',"{:14.10f}".format(zmax1)
   print '   SW-bal     ',"{:14.10f}".format(zmin2),'     ',"{:14.10f}".format(zmax2)
   print '   GN         ',"{:14.10f}".format(zmin3),'     ',"{:14.10f}".format(zmax3)
   print '   GN-bal     ',"{:14.10f}".format(zmin4),'     ',"{:14.10f}".format(zmax4)

   if option < 3:
      # For full and balanced fields, use same limits for SW & SW-bal, and
      # for GN and GN-bal:
      zminSW=min(zmin1,zmin2)
      zminGN=min(zmin3,zmin4)
      zmaxSW=max(zmax1,zmax2)
      zmaxGN=max(zmax3,zmax4)

      zmag=max(abs(zminSW),zmaxSW)
      zmin1=-zmag
      zmax1= zmag
      zmin2=-zmag
      zmax2= zmag

      zmag=max(abs(zminGN),zmaxGN)
      zmin3=-zmag
      zmax3= zmag
      zmin4=-zmag
      zmax4= zmag
   else:
      # For the imbalanced fields, use different limits for each simulation:
      zmag=max(abs(zmin1),zmax1)
      zmin1=-zmag
      zmax1= zmag

      zmag=max(abs(zmin2),zmax2)
      zmin2=-zmag
      zmax2= zmag

      zmag=max(abs(zmin3),zmax3)
      zmin3=-zmag
      zmax3= zmag

      zmag=max(abs(zmin4),zmax4)
      zmin4=-zmag
      zmax4= zmag

   print
   zmin1=float(raw_input('Minimum value to show in   SW   simulation? (default '+str(zmin1)+') ') or zmin1)
   zmax1=float(raw_input('Maximum value to show in   SW   simulation? (default '+str(zmax1)+') ') or zmax1)

   print
   zmin2=float(raw_input('Minimum value to show in SW-bal simulation? (default '+str(zmin2)+') ') or zmin2)
   zmax2=float(raw_input('Maximum value to show in SW-bal simulation? (default '+str(zmax2)+') ') or zmax2)

   print
   zmin3=float(raw_input('Minimum value to show in   GN   simulation? (default '+str(zmin3)+') ') or zmin3)
   zmax3=float(raw_input('Maximum value to show in   GN   simulation? (default '+str(zmax3)+') ') or zmax3)

   print
   zmin4=float(raw_input('Minimum value to show in GN-bal simulation? (default '+str(zmin4)+') ') or zmin4)
   zmax4=float(raw_input('Maximum value to show in GN-bal simulation? (default '+str(zmax4)+') ') or zmax4)

   # Obtain contour levels for plotting the colorbars:
   dz=contint(zmin1,zmax1)
   jmin=-int(-zmin1/dz)
   jmax= int( zmax1/dz)
   clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   dz=contint(zmin2,zmax2)
   jmin=-int(-zmin2/dz)
   jmax= int( zmax2/dz)
   clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   dz=contint(zmin3,zmax3)
   jmin=-int(-zmin3/dz)
   jmax= int( zmax3/dz)
   clevels3=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   dz=contint(zmin4,zmax4)
   jmin=-int(-zmin4/dz)
   jmax= int( zmax4/dz)
   clevels4=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

#==============================================================================
   # Set up figure:
   fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=[10,10], nrows=2, ncols=2)

   # Customise tick values:
   ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax4.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax4.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

   ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax4.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax4.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

   # Plot the images in an array with individual colourbars:
   im1=ax1.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax1)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
   setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=16)

   im2=ax2.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax2)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
   setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=16)

   im3=ax3.imshow(Z3,cmap=cm.seismic,vmin=zmin3,vmax=zmax3,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax3)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im3, cax=cax, ticks=clevels3)
   setp(cbar.ax.yaxis.set_ticklabels(clevels3), fontsize=16)

   im4=ax4.imshow(Z4,cmap=cm.seismic,vmin=zmin4,vmax=zmax4,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax4)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im4, cax=cax, ticks=clevels4)
   setp(cbar.ax.yaxis.set_ticklabels(clevels4), fontsize=16)

   # Add a title:
   fig.suptitle(title_label, fontsize=32)

   # Add labels for each simulation:
   ax1.set_title('SW', fontsize=20, fontname='Times New Roman')
   ax2.set_title('SW-bal', fontsize=20, fontname='Times New Roman')
   ax3.set_title('GN', fontsize=20, fontname='Times New Roman')
   ax4.set_title('GN-bal', fontsize=20, fontname='Times New Roman')

   ax1.label_outer()
   ax2.label_outer()
   ax3.label_outer()
   ax4.label_outer()

   # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
   plt.setp(ax1.get_xticklabels(), visible=False)
   plt.setp(ax2.get_xticklabels(), visible=False)

   #fig.subplots_adjust(wspace=0.8, hspace=0.5, top=0.6, bottom=0.05)
   fig.subplots_adjust(wspace=0.8, hspace=-0.1)

#=========================================================================
   # Save image:
   fig.savefig(outfile, format='eps', dpi=300)

   print
   print ' To view the image, type'
   print
   print ' gv ',outfile
   print
