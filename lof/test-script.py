from  erqml import *
from astropy.table import Table, Column
import numpy as np
import matplotlib.pyplot as plt
# Reading data
smp=Table.read('/home/reza/erq/sampling/org_sample2.fits')
iW3 = smp['i-w3']
kt80 = smp['kurt80_gf']
rew = smp['rew_gf']
rew  = np.log10(rew)
X=np.array(list(zip(iW3, rew)))


NN= int(0.5*len(X))
alpha=np.empty([9]); i=-1
y0=[]
for cont  in [0.001, 0.002, 0.005, 0.007, 0.01, 0.015]:
    i+=1
    print('-----',cont, '----')
    # ERQ, MCL, out, score, score_ERQ, score_in, score_out= LOF_Separater(X, cont, NN)
    # np.savetxt('bins/ERQ'+str(cont) + '.txt', ERQ) 
    # np.savetxt('bins/MCL'+str(cont) + '.txt', MCL)
    # np.savetxt('bins/OUT'+str(cont) + '.txt', out)
    # np.savetxt('bins/score'+str(cont) + '.txt', score)
    # np.savetxt('bins/score_out'+str(cont) + '.txt', score_out)
    # np.savetxt('bins/score_in'+str(cont) + '.txt', score_in)
    # np.savetxt('bins/score_ERQ'+str(cont) + '.txt', score_ERQ)

    ERQ = np.loadtxt('bins/ERQ'+str(cont) + '.txt') 
    MCL = np.loadtxt('bins/MCL'+str(cont) + '.txt')
    out= np.loadtxt('bins/OUT'+str(cont) + '.txt')
    score = np.loadtxt('bins/score'+str(cont) + '.txt')
    score_out= np.loadtxt('bins/score_out'+str(cont) + '.txt')
    # score_in = np.loadtxt('bins/score_in'+str(cont) + '.txt')
    score_ERQ = np.loadtxt('bins/score_ERQ'+str(cont) + '.txt')
    MCL_center=np.empty([2])
    MCL_center[0] =np.mean(MCL[:,0])
    MCL_center[1] =np.mean(MCL[:,1])
    ERQ_center=np.empty([2])
    ERQ_center[0] =np.mean(ERQ[:,0])
    ERQ_center[1] =np.mean(ERQ[:,1])
    # alpha[i], ERQ_center=wedge(ERQ, MCL_center, 0.90)
    path = 'cont-'+str(cont)+ '.pdf'
    # def boundary(sample, ERQ, MCL, out, y_range, score, score_ERQ, score_out, path):

    # y0=boundary(X, ERQ, MCL, out, [10,50], score, score_ERQ, score_out, path)
    
    # def binning(MCL_center, ERQ_center, n_bin, theta, sample, ERQ, MCL, path):
    theta, xx = wedge(ERQ, MCL_center, 0.90)
    binning(MCL_center, ERQ_center,  4, theta, X, ERQ, MCL, path)