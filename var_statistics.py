#!/usr/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

data = np.loadtxt('output')
plt.figure(figsize=(8,8),dpi=80)

ax = []
gs = gridspec.GridSpec(3, 1)
gs.update(left=0.08, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.0)
ax.append(plt.subplot(gs[0:2, 0]))
ax.append(plt.subplot(gs[2, 0]))

ax[0].set_ylim(0,110)
ax[0].set_yticks(np.arange(0,110,20))
ax[0].set_yticklabels(np.arange(0,110,20),fontsize=14)
ax[0].set_xlim(0,20)
ax[0].set_xticks(np.arange(0,21,2.5))
ax[0].set_xticklabels([])
ax[0].set_ylabel('S/N',fontsize=16)
ax[0].plot(data[:,0], data[:,1], 'ro', markersize=8, label=r'$(S/N)_{\rm{I}}$')
ax[0].plot(data[:,0], data[:,2], 'b*', markersize=8, label=r'$(S/N)_{\rm{var}}$')
#ax[0].plot(data[:,0], data[:,1], ls='-', lw=2, color='red', label=r'$\delta\nu_{\rm{DISS}}$ = 10')
#ax[0].plot(data[:,0], data[:,2], ls='--', lw=2, color='blue', label=r'$\delta\nu_{\rm{DISS}}$ = 10')

legend = ax[0].legend(loc='upper right',numpoints=1,fontsize=16)

#ax.text(2.8, 80, r'$\delta\nu=100$')

ax[1].set_ylim(2,8)
ax[1].set_yticks(np.arange(2,8,2))
ax[1].set_yticklabels(np.arange(2,8,2),fontsize=14)
ax[1].set_xlim(0,20)
ax[1].set_xticks(np.arange(0,21,2))
ax[1].set_xticklabels(np.arange(0,21,2),fontsize=14)
ax[1].set_ylabel(r'$(S/N)_{\rm{I}}/(S/N)_{\rm{var}}$',fontsize=16)
ax[1].set_xlabel('Noise level (arbitary units)',fontsize=16)
ax[1].plot(data[:,0], data[:,3], 'k+', markersize=8)
#ax[1].plot(data[:,0], data[:,3], ls='-', lw=2, color='black', label=r'$\delta\nu_{\rm{DISS}}$ = 10')
plt.savefig('varStatistics.ps',dpi=80)
plt.show()
