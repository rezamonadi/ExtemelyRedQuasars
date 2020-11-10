from astropy.table import Table, Column
import time
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
from scipy.stats import kde
from tqdm import tqdm          
sys.path.insert(0, '/home/reza/erq/')
from erqml import *
import os
from  sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde


# Reading data
# smp=Table.read('/home/reza/erq/sampling/org_sample2.fits')
smp=Table.read('/home/reza/erq/sampling/sample_Lum_flat.fits')
iW3_0 = smp['i-w3']
kt80_0= smp['kurt80_gf']
rew_0 = smp['rew_gf']
rew_0  = np.log10(rew_0)
frat_0 = smp['frat_nv/civ']
fwhm_0 = smp['fwhm_gf']
iz_0 = smp['i-z']
data = np.array(list(zip(iW3_0, rew_0, kt80_0)))
x0, y0, z0= np.median(data, axis=0)
x,y,z= data.T
start=0
ERQ = data[((iW3_0>4.6) & (rew_0>2))]
x_erq, y_erq, z_erq= np.median(ERQ, axis=0)
r_erq = np.sqrt((x_erq-x0)**2 + (y_erq-y0)**2 + (z_erq-z0)**2)
ngrid=100; ex=1.2
x_line =np.linspace(x0,x_erq*ex, ngrid)
y_line = np.linspace(y0,y_erq*ex, ngrid)
z_line= np.linspace(z0,z_erq*ex, ngrid)
r = np.sqrt((x_line-x0)**2+(y_line-y0)**2 + (z_line-z0)**2)
r/=r_erq
# # ---------------3D---------------------------------
u_erq = np.array([x_erq-x0, y_erq-y0, z_erq-z0])
u_erq = u_erq/np.linalg.norm(u_erq)


#  building random cone
nTrial=10
angle=np.empty(nTrial)
xf_cone=np.empty(nTrial)
yf_cone=np.empty(nTrial)
zf_cone=np.empty(nTrial)
Ax=0.01; Ay=0.05; j=0; s=25

for i in tqdm(range(nTrial)):


    dx, dy = (np.random.random_sample() -0.5)*Ax, (np.random.random_sample() -0.5)*Ay
    xf_cone[i] = x_erq+dx
    yf_cone[i] = y_erq+dy
    d = r_erq**2 -(xf_cone[i]-x0)**2 -(yf_cone[i]-y0)**2
    if (d>0):
        zf_cone[i] = np.sqrt(d) + z0
        u_cone = np.array([xf_cone-x0, yf_cone-y0, zf_cone-z0])
        u_cone=u_cone/np.linalg.norm(u_cone)
        # angle[i] = np.arccos(np.dot(u_erq, u_cone))
        print(i, np.arccos(np.dot(u_erq, u_cone)))


#   density in 3d  for different bandwidths 
for bw in [0.75, 1, 1.1, 1.25]:
    
    kde= gaussian_kde(data.T)
    kde.set_bandwidth(bw_method=kde.factor*bw)
    den_con=[]
    print(bw)
    for i in range(len(xf_cone)):
        x=np.linspace(x0,xf_cone[i]*ex, ngrid)
        y=np.linspace(y0,yf_cone[i]*ex, ngrid)
        z=np.linspace(z0,zf_cone[i]*ex, ngrid)
        line = np.array(list(zip(x,y,z)))
        density= np.log10(kde(line.T))
        den_con.append(density)


    den_med = np.median(den_con, axis=0)
    print(den_med.shape)
    plt.plot(r[s:-1], den_med[s:-1], alpha=0.7, label='bandwidth=%.2f' %(bw))
plt.legend()
plt.xlabel(r'r/r_{ERQ}')
plt.ylabel()
plt.savefig('cone-density.pdf')
plt.show()