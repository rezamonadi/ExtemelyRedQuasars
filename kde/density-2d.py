# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from astropy.table import Table, Column
import time
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import sys
from scipy.stats import gaussian_kde
from tqdm import tqdm          
sys.path.insert(0, '/home/reza/erq/')
from  erqml import *
import os
from matplotlib import cm
from  sklearn.neighbors import KernelDensity
# from scipy.stats import gaussian_kde
import matplotlib.ticker as ticker
from matplotlib import cm

# %matplotlib


# %%
## # Reading data
# smp=Table.read('/home/reza/erq/sampling/org_sample2.fits')
smp=Table.read('/home/reza/erq/sampling/LM_sample.fits')
iW3_0 = smp['i-w3']
kt80_0= smp['kurt80_gf']
rew_0 = smp['rew_gf']
rew_0  = np.log10(rew_0)
frat_0 = smp['frat_nv/civ']
fwhm_0 = smp['fwhm_gf']
# iz_0 = smp['i-z']
data_0=np.array(list(zip(iW3_0, rew_0)))
data = np.array(list(zip(iW3_0, rew_0)))
(data, minData, rangeData) =  MinMaxScaler(data)
# data, m, s= scale(data)
x0, y0= np.median(data, axis=0)
x,y= data.T
ERQ = data[((iW3_0>=4.6) & (rew_0>=2))]
x_erq, y_erq= np.median(ERQ, axis=0)  #center of ERQ population 
r_erq = np.sqrt((x_erq-x0)**2 + (y_erq-y0)**2)  # ERQ radius 
r_start = np.sqrt(((4.6-minData[0])/rangeData[0]-x0)**2 + ((2-minData[1])/rangeData[1]-y0)**2)  # ERQ radius 
print(r_start)


# %%

# #  Checking the extension 
# plt.scatter(data[:,0], data[:,1], alpha=0.1)
# plt.scatter(x0, y0, s=100)
# plt.scatter(x_erq, y_erq, s=100)
# plt.scatter(x_line, y_line, s=0.1, alpha=0.3)
# plt.axis('equal')

ti = time.time()
ngrid=200
Bw=[.18, 0.18*1.25, 0.18*1.5]
n=[1, 1.25, 1.5]
nn=-1
for bw in Bw:
    nn+=1
    enclosing_ratio=0.75; resolution=100; ex=1.4

    opening_angle, Wedge_direction = opening_angle_finder(ERQ, [x0,y0],     enclosing_ratio, resolution)
    d_theta= opening_angle
    nCone=100
    theta_cone = np.linspace(Wedge_direction-d_theta, Wedge_direction+d_theta, nCone)

    rf= r_erq*ex
    xf_cone = rf*np.cos(theta_cone) + x0
    yf_cone= rf*np.sin(theta_cone) + y0

    line_segments=100; 
    x_line = np.zeros([nCone, line_segments])
    y_line = np.zeros([nCone, line_segments])
    for i in range(nCone):
        x_line[i,:] =np.linspace(x0,xf_cone[i], line_segments) # x positions along a line as long as ERQ radius times the expansion which we calculate density 
        y_line[i,:] = np.linspace(y0,yf_cone[i], line_segments) # y positions ...
    r = np.sqrt((x_line-x0)**2+(y_line-y0)**2) # radius of each point needed for density 

    print(bw)
    fig = plt.figure()
    plt.clf()
    plt.cla()

    kde =gaussian_kde(data.T)
    kde.set_bandwidth(bw)
    # # Contour plot
    x, y = data.T
    xi, yi = np.mgrid[x.min():x.max():ngrid*1j, y.min():y.max():ngrid*1j]
    # xii, yii = np.mgrid[x.min()*s[0]+m[0]:x.max()*s[0]+m[0]:ngrid*1j,
    #                     y.min()*s[1]+m[1]:y.max()*s[1]+m[1]:ngrid*1j]


    zi = kde(np.vstack([xi.flatten(), yi.flatten()]))
    # zi = np.log(zi)
    zi/=np.max(zi)

    mask = zi >0.015
    zi[mask]=0.015

    # peak_ind = np.where(zi.reshape(xi.shape)==1)
    # x_peak = xi[peak_ind[0], peak_ind[1]]
    # y_peak = yi[peak_ind[0], peak_ind[1]]
    # c = axes.contour(xi, yi, zi.reshape(xi.shape), levels=[0.010])
    # plt.clabel(c, fontsize=15)

    # # plt.scatter(x_peak, y_peak, marker='X', s=100)


    # %%
    plt.clf()
    plt.cla()
    den_plt=plt.pcolormesh(xi, yi, 1000*zi.reshape(xi.shape), cmap=cm.coolwarm)
    cbar= plt.colorbar(den_plt)
    cbar.set_label(r'$\rho\times 10^{-3}$')
    # # den_plt=plt.c
    ax = fig.add_subplot(111)
    # ax.scatter(x,y, alpha=0.4, s=2, c='C0')
    # # # ax.scatter(ERQ[:,0], ERQ[:,1], alpha=0.4, c='r', s=2)
    # #     ax.plot([ex*np.array(xf_cone),x00], 
    # #                  [ex*np.array(yf_cone),y00], lw=.5, c='grey', alpha=0.5)
    plt.plot([xf_cone[0], x0], [yf_cone[0], y0], ls='--', lw=1, c='black' ,alpha=0.6)
    plt.plot([xf_cone[-1], x0], [yf_cone[-1], y0], ls='--', 
            lw=1, c='black', alpha=0.6 )
    plt.arrow(x0, y0, xf_cone[int(nCone*0.5)]-x0, yf_cone[int(0.5*nCone)]-y0,
            head_width=0.01, head_length=0.01, fc='black', ec='black')
    # # # plb.axis('equal')
    ticks_x = ticker.FuncFormatter(lambda x, 
                                    pos: '{0:g}'.format(round(x*rangeData[0]+ minData[0], 2)))
    ax.xaxis.set_major_formatter(ticks_x)

    ticks_y = ticker.FuncFormatter(lambda x, 
                                    pos: '{0:g}'.format(round(x*rangeData[1]+minData[1],2)))
    ax.yaxis.set_major_formatter(ticks_y)
    plt.xlabel('MMN(i-W3)')
    plt.ylabel('MMN(log10(REW))')
    # plb.title('Gaussian Kernel Bandwith: %.3f' %( Bw[nn]))
    # #     plb.savefig('/home/reza/erq/paper/draft/fig/cone_density-bw%2d.png' %(bw*100),
    # #                 dpi=1200, bbox_inches='tight', format='png')

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title(r'$\Omega=%.2f^{\circ}, h=%.2fh_{Scott}=%.2f$' %(2*np.rad2deg(opening_angle), n[nn], bw))
    # plt.scatter(iW3_0, rew_0, s=1, alpha=0.2)
    plt.savefig('cone_density-bw%2d.png' %(bw*100), dpi=1200, 
                bbox_inches='tight', format='png')
    # print(x0, (min(xf_cone)-max(xf_cone))/x0)
    plt.clf()
    plt.cla()

#######################################################################################################################
# enclosing_ratio=0.75; resolution=100; ex=1.3
# opening_angle, Wedge_direction = opening_angle_finder(ERQ, [x0,y0],     enclosing_ratio, resolution)
# print(Wedge_direction)
# d_theta=opening_angle
# nCone=100
# theta_cone = np.linspace(Wedge_direction-d_theta, Wedge_direction+d_theta, nCone)

# rf= r_erq*ex
# xf_cone = rf*np.cos(theta_cone) + x0
# yf_cone= rf*np.sin(theta_cone) + y0

# line_segments=100
# x_line = np.zeros([nCone, line_segments])
# y_line = np.zeros([nCone, line_segments])
# for i in range(nCone):
#     x_line[i,:] =np.linspace(x0,xf_cone[i], line_segments) # x positions along a line as long as ERQ radius times the expansion which we calculate density 
#     y_line[i,:] = np.linspace(y0,yf_cone[i], line_segments) # y positions ...
# r = np.sqrt((x_line-x0)**2+(y_line-y0)**2) # radius of each point needed for density 


# for bw in Bw:
        
#     kde =gaussian_kde(data.T)
#     kde.set_bandwidth(bw)
#     den_med=np.zeros([line_segments])
#     den_con=[]
#     den_max = max(kde(data.T))
#     for i in tqdm(range(nCone)):
#         line = np.array(list(zip(x_line[i,:], y_line[i,:])))
#         density= kde(line.T)
#         den_con.append(density) #density for each point along a line stored in a row of density_cone


#     den_med = np.median(den_con, axis=0) # median of each column gives the median density for a secific radius from center 

#     # for i in range(3): 
#     #     print(i)
#     plt.semilogy(r[0,:]/r_erq, den_med/den_max, alpha=0.7, 
#                 label=r'bandwidth: %.2f' %(bw))
# plt.xlabel(r'r/r$_{CERQ}$')
# plt.ylabel(r'$\rho$')
# plt.title(r'2D density profile, $\Omega$=%.2f$^{\circ}$'  %(2*np.rad2deg(opening_angle)))
# # plt.axvline(x=r_start/r_erq, ls='--', alpha=0.5, label='r(4.6,2)', c='black', lw=1)
# plt.legend()
# # plt.savefig('/home/reza/erq/paper/draft/fig/opp-cone-2d-density-d%2d.pdf' %(d_theta))
# plt.savefig('cone-2d-density-bw%2d.pdf' %(bw*100), bbox_inches='tight', format='pdf')



# # %%
