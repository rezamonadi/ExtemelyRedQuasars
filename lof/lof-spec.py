from astropy.table import Table, Column
import numpy as np
from random import shuffle
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm          
sys.path.insert(0, '/home/reza/erq/')
from erqml import *
import os

# Reading data
# smp=Table.read('/home/reza/erq/sampling/org_sample2.fits')
tab=Table.read('/home/reza/erq/sampling/sample_Lum_flat.fits')
iW3 = tab['i-w3']

kt80 = tab['kurt80_gf']
rew = tab['rew_gf']
rew  = np.log10(rew)
kt80 = tab['kurt80_gf']


z_dr12=tab['z_dr12']
plate = tab['Plate']
mjd = tab['MJD']
fiberid= tab['FiberID']
nqso=len(z_dr12)


# p=plt.scatter(iW3, rew, c=kt80, cmap='plasma', s=.1, alpha=0.3)
# cbar = plt.colorbar(p)
# cbar.set_label('kt80')

# plt.xlabel('color')
# plt.ylabel('REW')

# plt.savefig('kurt3d.pdf', format='pdf', bbox_inches='tight', dpi=1800)

# sys.exit()

X=np.array(list(zip(iW3, rew)))

center = np.empty([2])
center[0]=np.median(X[:,0])
center[1]=np.median(X[:,1])


# ##########################
ERQ=[]
for i in range(len(X)):
    if(X[i,0]>4.6) and (X[i,1]>2):
        ERQ.append(X[i,:])


# print(len(ERQ))
ERQ= np.array(ERQ)
# ###########################


# wedge_vec_angle=np.pi*2 -0.1
#  LOF scores 

NN = np.sqrt(len(X))
i=-1
# r= np.linspace(0.2, 0.3, 11)
# print(r)
# n = len(r)
# min_s=np.empty([n])
# mean_s=np.empty([n])
# max_s=np.empty([n])
# var_s=np.empty([n])
enclosing_ratio =0.9
nn=0.25
nbins=7

def doBin():
    for enclosing_ratio in [ 0.9]:
        opening_angle, wedge_vec_angle, ERQ_center = wedge(ERQ, center, enclosing_ratio)
        # print(wedge_vec_angle*180/np.pi)

        for nbins in [7]:
            for nn in [0.25]:

                print('lof.....')
                prediction, LofScores = LOF(X, 0.02, int(nn*len(X)))

                #  Triangle for wedge needed for binning 
                A, B = edger(X, center, wedge_vec_angle, opening_angle)



                #  LOF binning 
                print('binning...')
                binning_method= 'linear'
                bin_id, bin_population, LOF_score_edge = LofBinner(X, LofScores, center,A, B, nbins, binning_method)
                np.savetxt('MedianSpectra/' +str(enclosing_ratio)+'/lof_bin_id_nn'+str(nn) + '-nb-'+str(nbins)+'.txt', bin_id)
                np.savetxt('MedianSpectra/' +str(enclosing_ratio)+ '/lof_bin_pop_nn' + str(nn) + '-nb-' + str(nbins) +'.txt', bin_population)

                print(bin_population)
                print(max(bin_id))
                #  lof bin plotting
                print('plotting...')
                path = 'Bins/'+str(enclosing_ratio)+ '/lof-bin-' + str(nn)+   '-nb-' +str(nbins)+'.pdf'
                
                Lof_bin_plotter(X, bin_id, bin_population, LOF_score_edge, center, A, B, nbins, nn , path)


def doMed(nqso, enclosing_ratio, nbins, nn):
    #  median spectra 
    print('median spectra...')
    bin_id=np.loadtxt('MedianSpectra/' +str(enclosing_ratio)+'/lof_bin_id_nn'+str(nn) + '-nb-'+str(nbins)+'.txt')
    bin_population=np.loadtxt('MedianSpectra/' +str(enclosing_ratio)+ '/lof_bin_pop_nn' + str(nn) + '-nb-' + str(nbins) +'.txt')
    print(max(bin_id))
    print(bin_population)
    for i_bin in range(0,nbins):
        index_bin=[] # pool of indeces of qsos for each bin 
        nqso_bin = int(bin_population[i_bin])
        print('bin-pop', nqso_bin)

        for j in (range(nqso)):
            if(int(bin_id[j])==i_bin):
                index_bin.append(i_bin)
        z_dr12_bin=z_dr12[index_bin]
        plate_bin=plate[index_bin]
        mjd_bin=mjd[index_bin]
        fiberid_bin = fiberid[index_bin]

        MedianSpectrum=medspec(z_dr12_bin, plate_bin,mjd_bin, fiberid_bin)
        # plt.plot(MedianSpectrum)
        np.savetxt('MedianSpectra/'+str(enclosing_ratio) + '/MedianSpectrumlof-nn-'+str(nn)+'-nb-'+str(i_bin)+'.txt', MedianSpectrum)
        # plt.show()


# --------------------------------------------
#  Script 
# doBin()
doMed(nqso, enclosing_ratio, nbins, nn)
            
r=np.loadtxt('range.txt')
plt.plot()



