#!/usr/bin/env python3

# This script plots the results in davg.asc and uavg.asc, created by
# profile.f90 (which reads data created by pam.f90).

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# set tick label size:
label_size = 24
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 3

#=================================================================
# Set up figure:
fig, ((ax1, ax2, ax3)) = plt.subplots(figsize=[24.0,7.6], nrows=1, ncols=3)

ax1.set_xlabel('$\\hat{\\delta}_1,~\\hat{u}_1/(2\\Omega)$', fontsize=30)
ax1.set_ylabel('$a$', fontsize=30)
ax1.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax1.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax2.set_xlabel('$\\hat{\\delta}_2,~\\hat{u}_2/(2\\Omega)$', fontsize=30)
ax2.set_ylabel('$a$', fontsize=30)
ax2.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax2.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax3.set_xlabel('$\\hat{\\delta}_3,~\\hat{u}_3/(2\\Omega)$', fontsize=30)
ax3.set_ylabel('$a$', fontsize=30)
ax3.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax3.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

#=================================================================
# Open input files and read data:
in_file=open('d1.asc','r')
x1, y1 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
ax1.plot(x1,y1,c='k',lw=3)

in_file=open('u1.asc','r')
x1, y1 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
ax1.plot(x1,y1,c='b',lw=3)

in_file=open('d2.asc','r')
x1, y1 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
ax2.plot(x1,y1,c='k',lw=3)

in_file=open('u2.asc','r')
x1, y1 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
ax2.plot(x1,y1,c='b',lw=3)

in_file=open('d3.asc','r')
x1, y1 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
ax3.plot(x1,y1,c='k',lw=3)

in_file=open('u3.asc','r')
x1, y1 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
ax3.plot(x1,y1,c='b',lw=3)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0.1)

fig.savefig('evec.eps',  format='eps', dpi=300)

print (' Type:   gv evec.eps')
