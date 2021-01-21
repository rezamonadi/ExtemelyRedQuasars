import sys
sys.path.insert(0, '/home/reza/erq')

from scipy import interpolate
import pyfits
from numpy import *
import math
import scipy.ndimage
from scipy import interpolate
from numpy import nanmean
from numpy import nanmedian
from readSDSSspectrafast import *
import matplotlib.pylab as plt
import numpy as np
from erqml import *
from scipy import ndimage
from line_db import line_db
import os

## Define a log wavelength grid for the composite spectrum
step = 1.00015
bb = arange(0,8813,1)
wgrid = 800.0 * step**bb

# Insert code here to read data for all possible quasars (redshift, colors, line data, etc.)

from astropy.table import Table, Column

# selection parameters
                # /home/reza/erq/sampling/org_sample2.fits
tab=Table.read('/home/reza/erq/sampling/sample_Lum_flat.fits')
z_dr12=tab['z_dr12']
nqsos = len(z_dr12)
plate = tab['Plate']
mjd = tab['MJD']
fiberid= tab['FiberID']
W3 = tab['ABw3']
imw3 =tab['i-w3']
REW = tab['rew_gf']
FWHM = tab['fwhm_gf']
kt80 = tab['kurt80_gf']
frat_civ_nv = tab['frat_nv/civ']
# This 2D array will hold all the spectra used for the median below

# tip_id= np.loadtxt('/home/reza/erq/kde/kde-tip-label.txt')
# tip_pop = int(sum(tip_id))
# ex1=1.2
# ex2=1.9
# for nearest_neighbor_ratio in [0.1, 0.25, 0.5, 0.75, 0.9]:
enclosing_ratio =0.9; nn=0.25; nbins=7; nearest_neighbor_ratio = 0.25
bin_id = np.loadtxt( '/home/reza/erq/lof/MedianSpectra/' +str(enclosing_ratio)+'/lof_bin_id_nn'+str(nn) + '-nb-'+str(nbins)+'.txt')
bin_pop= np.loadtxt('/home/reza/erq/lof/MedianSpectra/' +str(enclosing_ratio)+ '/lof_bin_pop_nn' + str(nn) + '-nb-' + 
str(nbins) +'.txt')
c=[ 'black', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'r']
print(nqsos, len(bin_id))
ymax=9; fs=10; fs1=9
# ymax=8.0; fs=10; fs1=9
for ii in range(0,3):

    if(ii==0): l = 1000; u = 1980
    if(ii==1): l = 1150; u = 1980
    if(ii==2): l = 1150; u = 2850
    fig = plt.figure(figsize=(11.5,5.6))
    plt.xlim(l,u)
    # plt.ylim(0.1*ymax, ymax)
    plt.ylim(-0.05*ymax, ymax)
    # plt.axes().set_aspect('equal')
    w3_med=[];imw3_med=[]; REW_med=[]; fwhm_med=[]; kt80_med=[]; frat_med=[]
    
    
    line_db(ymax, fs, fs1)

    print(len(bin_pop))
    for b in range(0,nbins):
        if(bin_pop[b]>1):   
            print('bin:', b, l, 'to', u)
            w3_med=[];imw3_med=[]; REW_med=[]; fwhm_med=[]; kt80_med=[]; frat_med=[]

            for i in range(nqsos):
                if(bin_id[i]==b+1):
                    w3_med.append(W3[i]) 
                    imw3_med.append(imw3[i])
                    REW_med.append(REW[i])
                    fwhm_med.append(FWHM[i])
                    kt80_med.append(kt80[i])
                    frat_med.append(frat_civ_nv[i])
                            
            med1=loadtxt('/home/reza/erq/lof/MedianSpectra/'+str(enclosing_ratio) + '/MedianSpectrumlof-nn-'+str(nn)+'-nb-'+str(b)+'.txt')
            plt.ylabel('Normalized Flux')
            plt.xlabel(r'$\lambda$')
            
            if(bin_pop[b]<100):
                sm_med1 = ndimage.filters.gaussian_filter1d(med1,5.0)
            else:
                sm_med1=med1


            ind = (wgrid>l) & (wgrid<u)

            plt.plot(wgrid[ind],  sm_med1[ind], lw = 1, c=c[b])
            plt.text(1280, 0.83*ymax -(b+1)*0.7, 'bin%d  #=%7d    W3=%.2f  i-w3=%.1f  EW=%3.0f  FW=%.0f  kt80=%.2f  nv/civ=%.2f' % (b+1, bin_pop[b], np.median(w3_med),  np.median(imw3_med), np.median(REW_med), \
            np.median(fwhm_med),   np.median(kt80_med), np.median(frat_med) ), ha='left', va='center', color=c[b], fontsize=10)

            # plt.text(1280, 0.83*ymax -(b+1)*0.2, 'bin%d  #=%7d    W3=%.2f  i-w3=%.1f  EW=%3.0f  FW=%.0f  kt80=%.2f nv/civ=%.2f' % (b+1, bin_pop[b], \
            # plt.text(1280, 0.83*ymax -(b+1)*0.35, 'bin%d  #=%7d    W3=%.2f  i-w3=%.1f  EW=%3.0f  FW=%.0f  kt80=%.2f nv/civ=%.2f' % (b+1, bin_pop[b], \
            # np.median(w3_med),  np.median(imw3_med), np.median(REW_med), np.median(fwhm_med),   np.median(kt80_med), np.median(frat_med) ), ha='left', va='center', color=c[b], fontsize=9)
    
    # tip_spec = loadtxt('/home/reza/erq/med-spec/kde/medspec-tip-kde.txt')
    
    # plt.plot(wgrid[ind],  tip_spec[ind], lw = 1, label = 'typical, #15929', c = 'gray', alpha=0.5)

    
    # plt.text(2800,  0.9*ymax,'MGII',ha='center',va='center',fontsize=fs)
                # /home/reza/erq/kde/1.31/1.31kde-bin-wedge-192-1.31.pdf
    plt.savefig('LofMedianSpectra-nn-'+str(nearest_neighbor_ratio) + '-'+str(l)+'-to-'+str(u)+'.pdf', format='pdf',  bbox_inches='tight')
    # plt.savefig('kde-med-spec-'+ 'wedge-'+str(wedge) + '-'+str(l)+'-to-'+str(u), format='pdf')
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()
