import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
from sklearn.decomposition import PCA
from astropy.table import Table, Column


def scale(x):

        y= x - np.mean(x)
        y = y/np.std(x)
        return y





# reading ...
smp=Table.read('LumMatch/LumMatch.fits') 

iW3 = smp['i-w3']
# iW3 = 10.0**(iW3/2.5)
kt80 = smp['kurt80_gf']
rew = smp['rew_gf']
rew=np.log10(rew)

iW3_sc = scale(iW3)
rew_sc = scale(rew)
kt80_sc = scale(kt80)



X = np.array(list(zip(iW3, rew))) #, kt80, fwhm, N5_C4, iz)))
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca.fit(X)


def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)

# plot data
plt.scatter(X[:, 0], X[:, 1], alpha=0.2, s=1)
for length, vector in zip(pca.explained_variance_, pca.components_):
    v = vector *4* np.sqrt(length)
    draw_vector(pca.mean_, pca.mean_ + v)
# plt.axis('equal');

MainCenter = np.mean(X, axis=0)
T1CERQ = X[(iW3>=4.6) & (rew>=2)]
T1CERQ_Center = np.median(T1CERQ, axis=0)
# C_cerq = 
print(T1CERQ_Center)
plt.arrow(MainCenter[0],MainCenter[1], T1CERQ_Center[0]-MainCenter[0],
          T1CERQ_Center[1] -MainCenter[1] ,  head_width=0.05, head_length=0.1, 
          fc='lightblue', ec='black')
plt.annotate('C$_{T1LM}$', (MainCenter[0]-0.1,MainCenter[1]+0.1),fontsize=10)
plt.annotate('C$_{T1CERQ}$', (T1CERQ_Center[0]+0.1,T1CERQ_Center[1]),fontsize=10)
plt.scatter(T1CERQ[:, 0], T1CERQ[:, 1], alpha=1, s=1)

plt.ylabel('log10(REW)')
plt.xlabel('i-w3')
plt.axis('equal')
plt.savefig('pca.pdf')

plt.show()