#!/usr/bin/python
import math
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

x0 = [120.0]
y0 = [2.0]

data0 = np.loadtxt('out1')
m,n = data0.shape
#print m,n
#z = np.arange(0,1000).reshape((10,100))
#x = np.arange(0,100)
#y = np.arange(0,100)

data = np.sort(data0.view('f8,f8,f8,f8,f8'), order=['f0','f1'], axis=0)

#for i in np.arange(0,m):
#	print data['f0'][i],data['f1'][i],data['f2'][i]
z = data['f2'].reshape((54,55))
x = data['f0'].reshape((54,55))
y = data['f1'].reshape((54,55))

#print x[:,0],y[0,:]
ax1 = plt.subplot(111)
ax1.set_xlim(10,5000)
ax1.set_ylim(10,5000)
ax1.set_xlabel(r'$\delta t_{\rm{ISS}}\,(\rm{arbitary unit})$',fontsize=14)
ax1.set_ylabel(r'$\delta\nu_{\rm{ISS}}\,(\rm{arbitary unit})$',fontsize=14)
#plt.hist(valHis, 100, normed=1, facecolor='blue', alpha=0.5)
#plt.hist(valNHis, 40, normed=1, facecolor='red', alpha=0.75)
#plt.plot(xp,yp,'-',color='blue',label='PSR')

#legend = ax.legend(loc='upper right', numpoints=1)
#plt.imshow(z)
#plt.pcolor(x,y,z,norm=LogNorm(vmin=z.min(), vmax=z.max()))
print z.min(), z.max()
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xticks([10,100,1000])
ax1.set_xticklabels([10,100,1000],fontsize=12)
ax1.set_yticks([10,100,1000])
ax1.set_yticklabels([10,100,1000],fontsize=12)
#plt.imshow(z)
plt.pcolor(x,y,z,norm=LogNorm(vmin=0.6, vmax=z.max()))
#plt.pcolor(x,y,z,norm=LogNorm(vmin=z.min(), vmax=z.max()))
#plt.pcolormesh(x,y,z,norm=LogNorm(vmin=z.min(), vmax=z.max()))
plt.title('Detection Sensitivity Map')
cbar = plt.colorbar()

ticks = []
ticks0 = [0.6,1,2,3,4]
for i in np.arange(0,5):
	ticks.append((math.log10(ticks0[i])-math.log10(0.6))/(math.log10(z.max())-math.log10(0.6)))
#ticks = np.log10(np.arange(1,13))/math.log10(z.max())
cbar.ax.get_yaxis().set_ticks(ticks)
#cbar.ax.get_yaxis().set_ticks([0,math.log10(5.0)/math.log10(z.max()),1.0/math.log10(z.max())])
cbar.ax.set_yticklabels([0.6,'1 (noise)',2,3,4], rotation=270)
cbar.set_label('Flux density', rotation=270, fontsize=14)

plt.scatter(x0,y0,color='red',marker='o',s=20,edgecolor='none')
#plt.savefig('histogram.ps',dpi=80)
plt.show()
