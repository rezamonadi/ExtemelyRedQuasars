from erqml import *
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

fig, ax=plt.subplots()
p=plt.scatter(iW3, rew, c=kt80, cmap='coolwarm', s=3, alpha=0.7)
cbar=plt.colorbar(p)
cbar.set_label('kt80')
plt.xlabel('i-W3')
plt.ylabel('log10(REW)')
plt.savefig('kt80-cmap.pdf', dpi=1200, bbox_inches='tight')
plt.show()
#--------------------------------------

# LOF nearest neighbors input
NN= int(0.5*len(X))
# wedge of ERQ angle initialization 
alpha=np.empty([9])
i=-1

for cont in [0.001, 0.002, 0.005, 0.007, 0.01, 0.015]:
    print('----',cont, '----')
    # Seprating ERQs and getting scores of LOF
    ERQ, MCL, outlier, sample_score, score_ERQ, score_in, score_out = LOF_Separater(X, cont, NN)

    # plotting the boundary between ERQs and inliers
    boundary(X, ERQ, MCL, outlier, [10,50], sample_score, score_ERQ, score_out, 'boundary-cont' + str(cont) + '.pdf')

    # Binning
    bin_label, bin_pop = binning(5, X, ERQ, MCL, 0.9, 'bin-cont'+str(cont) + '.pdf')
    np.savetxt('bin-label-'+str(cont) + '.txt', bin_label)
    np.savetxt('bin-pop-'+str(cont) + '.txt', bin_pop)