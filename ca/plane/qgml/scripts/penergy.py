#!/usr/bin/python3

#=====================================================================
# Plots the kinetic, potential and total energy as a function of time.
#=====================================================================

#=====perform various generic imports=====
import warnings
import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})

## global settings

# set tick label size:
label_size = 25
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

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)
#=========================================

# Read data and plot:
in_file=open('evolution/energy.asc','r')
t, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

#-------------------------------------------------------------------------
# Set up figure:
fig = plt.figure(1,figsize=[12,5])
ax1 = fig.add_subplot(111)

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('Energy', fontsize=30)

ax1.plot(t,ekin,color='b',lw=2,label='$\\mathcal{K}$')
ax1.plot(t,epot,color='r',lw=2,label='$\\mathcal{P}$')
ax1.plot(t,etot,color='k',lw=2,label='$\\mathcal{K}+\\mathcal{P}$')

ax1.legend(loc='best',prop={'size':25})

plt.savefig('energy.png', bbox_inches='tight', pad_inches = 0.025, dpi=300)
