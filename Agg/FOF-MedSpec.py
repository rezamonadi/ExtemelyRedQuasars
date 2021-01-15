#!/usr/bin/env python
# coding: utf-8

# In[10]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
from sklearn.decomposition import PCA
from astropy.table import Table, Column
import sys
sys.path.insert(0, '/home/reza/erq/')
from erqml  import *
from sklearn.cluster import AgglomerativeClustering
from sklearn import svm
import sys
sys.path.insert(0, '/home/reza/erq')

from scipy import interpolate
import pyfits
import math
import scipy.ndimage
from scipy import interpolate
from numpy import nanmean
from numpy import nanmedian
from readSDSSspectrafast import *
from erqml import *
from tqdm import notebook
import os
import sys
sys.path.insert(0, '/home/reza/erq')
from scipy import interpolate
import pyfits
from numpy import *
import math
import scipy.ndimage
import numpy as np
from erqml import *
from scipy import ndimage
from line_db import line_db
import os


# In[ ]:


def scale(x):

        y= x - np.mean(x)
        y = y/np.std(x)
        return y, np.mean(x), np.std(x)


# In[3]:


# reading ...
smp=Table.read('/home/reza/erq/sampling/sample_Lum_flat.fits') 

iW3 = smp['i-w3']
# iW3 = 10.0**(iW3/2.5)
kt80 = smp['kurt80_gf']
rew = smp['rew_gf']
rew=np.log10(rew)
fwhm = smp['fwhm_gf'] 
fwhm=np.log10(fwhm)
N5_C4 = smp['frat_nv/civ']
iz=smp['i-z']
sdss_name=smp['sdss_name']
z_dr12=smp['z_dr12']
plate = smp['Plate']
mjd = smp['MJD']
fiberid= smp['FiberID']

iW3_sc,m1,s1 = scale(iW3)
rew_sc,m2,s2 = scale(rew)
kt80_sc,a,b = scale(kt80)
fwhm_sc,a,b = scale(fwhm)
N5_C4_sc,a,b = scale(N5_C4)
iz_sc,a,b = scale(iz)


X = np.array(list(zip(iW3, rew)))
# X = np.array(list(zip(iW3, rew, kt80, kt80 , fwhm, N5_C4, iz)))
X_sc = np.array(list(zip(iW3_sc, rew_sc, kt80_sc , fwhm_sc, N5_C4_sc, iz_sc)))
# X_sc = np.array(list(zip(iW3_sc, rew_sc  )))

# X_tr = np.delete(X_sc, j, axis=1)
X_tr = X_sc[:,0:2]
ERQ_sc = X_tr[((iW3>4.6) & (rew>2))]
ERQ = X[((iW3>4.6) & (rew>2))]
link='complete'


# In[4]:


ncl =8 # number of clusters 
print(ncl, link)
model = AgglomerativeClustering(n_clusters=ncl, linkage=link, 
                                    affinity='euclidean').fit(X_tr)
labels = model.labels_


# In[6]:


colors = ['C0', 'C1', 'C2','C3', 'C4', 'C5', 'C6', 'C7', 'C8']
color_code=[]
for i in range(len(X_sc)):
    color_code.append(colors[labels[i]])


# # calculting median spectra

# # in[7]:


cl_pop=np.zeros([ncl])
cl_pop=np.array(cl_pop, dtype=np.int64)
for i in range(len(labels)):
    for j in range(ncl):
        if labels[i]==j:
            cl_pop[j]+=1
print(cl_pop)
# plt.clf()
# plt.cla()
# plt.scatter(x[:,0], x[:,1], c=color_code, alpha=0.7, s=4)
# plt.xlabel('i-w3')
# plt.ylabel('log10(rew)')
# plt.savefig('agg-8-cl.pdf')


# In[ ]:



# ## Define a log wavelength grid for the composite spectrum
# step = 1.00015
# bb = arange(0,8813,1)
# wgrid = 800.0 * step**bb
# nw = len(bb)

# # Insert code here to read data for all possible quasars (redshift, colors, line data, etc.)

# from astropy.table import Table, Column

# # selection parameters
# tab=Table.read('/home/reza/erq/sampling/sample_Lum_flat.fits')



# # nqsos is the total number of quasars in the catalog.
# # zdr12 is the redshift from the DR12 catalog. You might use something else.
# # zp1 is z+1

# # This loop goes through the quasar selection only to determine the total number that will go into the median
# # nqsos = len(z_dr12)
# nqsos = len(z_dr12)
# # for i in range(nqsos):
# # # Example criteria for EXCLUDING quasars:
# #     if (fwhm[i] <= 2000) or ((imw3[i]) < - 0.4) or (z_dr12[i]<2.0) or (z_dr12[i]>3.4) or \
# #             (ABw3_snr[i]<=3) or (ratio_rew[i]<=4.0) or (ratio_fwhm[i]<=5) or (qflag[i]!=0) or (nvflag[i]!=0) \
# #             or (bal_flag_vi[i] != 0) or (bi_civ[i] > 500 and bi_civ[i]/err_bi_civ[i] > 3): continue
# #     nspec = nspec + 1
# #     if nspec >5 : break  # uncomment this for testing to limit number contributing to median
# # print('nspec =',nspec)

# # This 2D array will hold all the spectra used for the median below

# for cl_label in range(8):
#     print('cl_label:', cl_label)
            
#     sp = zeros([int(cl_pop[cl_label]), nw])
#     k=-1
#     for i in notebook.tqdm(range(nqsos)):
#         if(labels[i]==cl_label):
#             k+=1
#         # Retrieve the spectra:
#             file = '/home/reza/erq/fred/sdss/%d/spec-%d-%d-%04d.fits'             % (plate[i], plate[i],mjd[i],fiberid[i])
#             spec = readDR10spec(file)
#             wave = spec['wl']
#             wz = wave/(z_dr12[i]+1)
#             flux = spec['flux']
#             mask = (wz > 1680.0) & (wz < 1730.0)
#             fnorm = median(flux[mask])
#             fluxn = flux/fnorm
#         # interpolate the rest-frame spectrum onto the standard grid
#             f = interpolate.interp1d(wz,fluxn,bounds_error=False,fill_value=float('nan'))
#             sp[k] = f(wgrid)
#         # calculate the median spectrum

#     med1 = nanmedian(sp,axis=0)
#     savetxt( 'agg-medspec-nCl-8-cl_label-%d.txt' %(cl_label), med1)



# Plotting the median Spectra for each cluster

# In[22]:




## Define a log wavelength grid for the composite spectrum
step = 1.00015
bb = arange(0,8813,1)
wgrid = 800.0 * step**bb

# Insert code here to read data for all possible quasars (redshift, colors, line data, etc.)

from astropy.table import Table, Column
# SMALL_SIZE = 20
# MEDIUM_SIZE = 20
# BIGGER_SIZE = 20

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# # 



# selection parameters
tab=Table.read('/home/reza/erq/sampling/sample_Lum_flat.fits')
# This 2D array will hold all the spectra used for the median below
# selection parameters
z_dr12=tab['z_dr12']
plate = tab['Plate']
mjd = tab['MJD']
fiberid= tab['FiberID']


# nqsos is the total number of quasars in the catalog.
# zdr12 is the redshift from the DR12 catalog. You might use something else.
# zp1 is z+1

# This loop goes through the quasar selection only to determine the total number that will go into the median
# nqsos = len(z_dr12)
nqsos = len(z_dr12)

# This 2D array will hold all the spectra used for the median below
# while(wedge<360):


#     bin_label=np.loadtxt('bw/enc-r-%.2f/%.1f/bin_label-bw-%.1f.txt' %(enclosinf_ratio, bw, bw))
#     tip_label=np.loadtxt('bw/enc-r-%.2f/%.1f/tip_label-bw-%.1f.txt' %(enclosinf_ratio, bw, bw) )
#     bin_pop=np.loadtxt('bw/enc-r-%.2f/%.1f/bin_pop-bw-%.1f.txt' %(enclosinf_ratio, bw, bw))
#     tip_pop = int(sum(tip_label))
c=['C0', 'C1', 'C2','C3', 'C4', 'C5', 'C6', 'C7',  'C8']
j=0
if(j==0): ymax=19.5; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==1 and enclosing_ratio==0.85): ymax=18.0; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==1 and enclosing_ratio==0.9): ymax=19; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==2): ymax=14; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==2 and enclosing_ratio==0.95): ymax=14; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==1 and enclosing_ratio==0.95): ymax=18; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==3): ymax=11; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.7
# if(j==4): ymax=6; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.3
# if(j==5): ymax=5; ymin=0.01*ymax; fs=10; fs1=9; x_median=1280; y_median = 0.83*ymax; med_step=0.3

for ii in range(3,4):

#     if(ii==0): lambda_min = 1000; lambda_max = 1980
#     if(ii==1): lambda_min = 1150; lambda_max = 1980
#     if(ii==2): lambda_min = 1150; lambda_max = 2850
    if(ii==3): lambda_min = 1150; lambda_max = 1700
    fig = plt.figure(figsize=(11.5,5.6))
#     ax=fig.add_subplot(111)
    # plt.axes().set_aspect('equal')
    # plt.xlim(lambda_min, lambda_max)
    plt.ylim(ymin,ymax)
    line_db(ymax, fs, fs1)
   
    
    
    for cl_label in range(7):
            
        print('Cl:', cl_label, lambda_min, 'to', lambda_max)

        med1=loadtxt('agg-medspec-nCl-8-cl_label-%d.txt' %(cl_label))
        plt.ylabel('Normalized Flux')
        plt.xlabel(r'$\lambda  (\AA)$')

        if(cl_pop[cl_label]<100):
            sm_med1 = ndimage.filters.gaussian_filter1d(med1,2.0)
        else:
            sm_med1=med1

        # masking for plot ranges
        ind = (wgrid>lambda_min) & (wgrid<lambda_max)
        plt.plot(wgrid[ind],  sm_med1[ind], lw = 1,  c=colors[cl_label],
                 label='cl-%d: #%d' %(cl_label+1, cl_pop[cl_label]))
    plt.legend(loc=7)
    
    plt.savefig('agg-nCl-8-med-spec-%d-%d.pdf' %(lambda_min, lambda_max), format='pdf',  bbox_inches='tight')
    plt.savefig('/home/reza/erq/paper/draft/fig/agg-nCl-8-med-spec-%d-%d.pdf' %(lambda_min, lambda_max), format='pdf',
                bbox_inches='tight')
   



