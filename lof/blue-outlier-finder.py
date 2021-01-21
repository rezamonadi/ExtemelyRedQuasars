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


NN= int(0.5*len(X)); cont=0.01

ERQ, MCL, out, blue, blue_ind, score, score_ERQ, score_in, score_out, score_blue= LOF_Separater(X, cont, NN)
np.savetxt('blue-cont-0.01.txt', blue_ind)
