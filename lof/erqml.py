"""
This module contains the functions needed for ERQ clustering
"""






def scale(x):
        """
        Scaling to mean and variance
        output has zero mean and variance of 1
        x is an N-dim array and output will be overwritten on x
        """
        import numpy as np


        dim = x.shape[1]
        y = np.empty([len(x), dim])
        mean= np.empty([dim])
        var= np.empty([dim])
        for i in range(0,dim):
                mean[i] = np.mean(x[:,i])
                var[i] = np.std(x[:,i])
                y[:,i]=(x[:,i] - mean[i])/var[i]

        return [y, mean, var]

def scale_inv(x, mean, var):
        """
        Scaling back to mean and variance
        output has zero mean and variance of 1
        x is an N-dim array and output will be overwritten on x
        """
        import numpy as np
        y=np.empty([len(x)])
        for i in range(len(y)):
                y[i]=x[i] *var + mean

        return y




def dbf(X, e, minPts, plotting, path,dim, verbose):
        """
        DBSCAN clustering
        inputs: -search radius -> e
                -minimum number of points in a cluster
                -an nd-array -> X
        """
        from sklearn.cluster import DBSCAN
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        db = DBSCAN(eps=e, min_samples=minPts, algorithm= 'auto', metric='euclidean' , leaf_size= 30, n_jobs= -1).fit(X)

        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        if(plotting==True):

                # dim = input('2 for 2d or 3 for 3d')
                colors = set_colors(labels)
                if(dim==3):
                        fg =plt.figure()
                        ax= Axes3D(fg)
                        for i in range(len(labels)):
                                ax.scatter(X[i,0], X[i,1], X[i,2],  s=5, color=colors[i] )
                else:
                        for i in range(len(labels)):
                                plt.scatter(X[i,0], X[i,1],   s=5, color=colors[i] )

                plt.savefig(path)
                plt.close()
                print('r->', e, 'mpt->', minPts)
        if(verbose==True):
                if(n_clusters_>0):
                        count=np.zeros([n_clusters_])
                        for i in range(len(labels)):
                                for j in range(0,n_clusters_):
                                        if(labels[i]==j):
                                                count[j]+=1
                        print('nCL: ', n_clusters_)
                        print('populations:', count)
                else:
                        print('No cluster found!')
        return [n_clusters_, labels]


def db_separation(X, e, minPts, verbose, plotting):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import numpy as np

        dim = X.shape[1]
        X_sc = scale(X)
        n_cl, labels= dbf(X_sc, e, minPts)

        count=np.zeros([n_cl])
        for i in range(len(labels)):
                for j in range(0,n_cl):
                        if(labels[i]==j):
                                count[j]+=1
        ind_count_sorted = np.argsort(count)
        ind_count_sorted = ind_count_sorted[::-1]
        n_main_cluster =0; n_noise=0; n_noise_db=0
        for i in range(len(labels)):
                if(labels[i]==ind_count_sorted[0]):
                        n_main_cluster+=1
                else:
                        n_noise+=1
                if(labels[i]==-1):
                        n_noise_db+=1
        noise_scaled = np.empty([n_noise, dim])
        noise = np.empty([n_noise, dim])
        main_cluster_scaled = np.empty([n_main_cluster,dim])
        main_cluster = np.empty([n_main_cluster,dim])

        i_noise=0; i_main_cluster=0
        for i in range(len(labels)):
                if(labels[i]==ind_count_sorted[0]):
                        main_cluster[i_main_cluster,:] = X[i,:]
                        main_cluster_scaled[i_main_cluster,:] = X_sc[i,:]
                        i_main_cluster+=1
                else:
                        noise[i_noise,:]=X[i,:]
                        noise_scaled[i_noise,:] = X_sc[i,:]
                        i_noise+=1

        if (verbose==True):
                        print('number of clusters: ', n_cl)
                        print('pop. of each cluster', count)
                        print('pop. of noise', n_noise_db )
        if(plotting == True):
                # path = input('fig name?' )
                path ='2d-rew-kt80'
                dim = 2
                # dim = input('2 for 2D, 3 for 3D')
                lx = 'iw3'
                # lx = input('xlabel?' )
                ly ='kt80'
                # ly = input('y label?' )
                # d1 = int(input('x axis?'))
                d1 = 0
                # d1 = int(input('x axis?'))
                d2 = 1
                # dim = X.shape[1]
                # d1=0; d2=1; d3=2
                if(dim==3):
                        d3 = int(input('z axis?' ))
                        lz= input('z label?' )
                        fg = plt.figure()
                        ax = Axes3D(fg)
                        ax.scatter(main_cluster[:,d1], main_cluster[:,d2], main_cluster[:,d3],  s=2, color='b', marker = 'x', alpha=1 )
                        ax.scatter(noise[:,d1], noise[:,d2], noise[:,d3],  s=2, color='C6' )
                        ax.set_xlabel(lx)
                        ax.set_ylabel(ly)
                        ax.set_zlabel(lz)
                        fg.savefig(path + '.pdf')

                else:
                        plt.scatter(noise[:,d1], noise[:,d2],  s=2, color='C6')
                        plt.scatter(main_cluster[:,d1], main_cluster[:,d2],  s=2, color='b', alpha=1, marker='x')
                        plt.xlabel(lx)
                        plt.ylabel(ly)
                        plt.savefig(path+'.pdf')

                plt.show()
                plt.close()



        return [noise, noise_scaled, main_cluster, main_cluster_scaled]

def set_colors(labels):
        colors=[ 'b', 'g','C5', 'C8', 'C3', 'C2', 'C1', 'C4', 'C7',
                'C5', 'C6', 'b', 'g','C5', 'C8', 'C3', 'C2', 'C1', 'C4', 'C7', 'C5', 'C6']

        colored_labels = []
        for label in labels:
                colored_labels.append(colors[label])

        return colored_labels



def plot_db(labels, n_cl, x,ERQ_cluster_number,  path):

        """
        plotting clustering results based on DBSCAN algorithm
        """
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        import numpy as np
        colors = set_colors(labels)
        fg = plt.figure()

        dim=np.shape(x)[1]
        if (dim==3):
                d1 = 0
                d2 = 1
                lx ='color'
                ly = 'rew'
                d3 = 2
                lz= 'kt80'

                ax = Axes3D(fg)
                for i in range(len(labels)):
                        if (labels[i]==ERQ_cluster_number):
                                ax.scatter(x[i,d1], x[i,d2], x[i,d3],  s=10, color=colors[i], marker = 'x', alpha=1 )

                        else:
                                ax.scatter(x[i,d1], x[i,d2], x[i,d3],  s=10, color=colors[i],  alpha=0.6 )
                        if(x[i,d1]>=4.6 and x[i,d2]>=100):
                                ax.scatter(x[i,d1], x[i,d2], x[i,d3], marker='o', s=70, facecolors='none', edgecolors='#000000')
                ax.set_xlabel(lx)
                ax.set_ylabel(ly)
                ax.set_zlabel(lz)
        else:

                for i in range(len(labels)):
                        if (labels[i]==ERQ_cluster_number):

                                plt.scatter(x[i,d1], x[i,d2],  s=12, color=colors[i], alpha=1, marker='x')
                        else:
                                plt.scatter(x[i,d1], x[i,d2],  s=4, color=colors[i], alpha=0.6)
                        plt.xlabel(lx)
                        plt.ylabel(ly)


        fg.savefig(path + '.pdf')
        plt.show()
        plt.close()

        return 0

def ERQ_finder( sample, sample_scaled, e, minPts,  plotting, verbose ):
        """
        Finds ERQ cluster in noise:
        the biggest cluster that is on average redder than i - w3 = 4.6
        is recognized as ERQ cluster
        """

        import numpy as np
        import sys
        # dbscan
        n_cl, labels= dbf(sample_scaled, e, minPts)
        dim = sample.shape[1]
        # Deciding which cluster is ERQ
        iW3_mean=np.zeros([n_cl])
        count=np.zeros([n_cl])

        for j in range(n_cl):
                for i in range(len(labels)):
                        if(labels[i]==j):
                                iW3_mean[j] += sample[i,0]
                                count[j]+=1

        iW3_mean=iW3_mean/count
        sample_size=len(sample)

        index_sorted=np.argsort(count)
        index_sorted=index_sorted[::-1]



        # i_w3_max=0
        # for i in index_sorted[0:2]:
        #         if(iW3_mean[i]>i_w3_max):
        #                 i_w3_max=iW3_mean[i]
        #                 ERQ_cluster_number = i
        n_sample = len(sample)
        ERQ_cluster_number=-1
        for i in index_sorted:
                if((iW3_mean[i]>4.0) and (float(count[ERQ_cluster_number])/float(n_sample)> 0.01 )):
                        ERQ_cluster_number = i

        if (ERQ_cluster_number>-1):

                ERQ=np.empty([int(count[ERQ_cluster_number]),dim])
                ERQ_scaled=np.empty([int(count[ERQ_cluster_number]),dim])
                i_ERQ=0
                ind_ERQ=[]
                for i in range(len(labels)):
                        if(labels[i]==ERQ_cluster_number):
                                ERQ[i_ERQ,:] = sample[i,:]
                                ERQ_scaled[i_ERQ,:] = sample_scaled[i,:]
                                i_ERQ+=1
                                ind_ERQ.append(i)


        if(verbose==True):
                print('-----ERQ cluster------')
                print('Search rad=', e, 'MinPts=', minPts)
                print('DBSCAN on the noise points')
                print('---------------------------------------------')
                print('number of clusters in noise:', n_cl)
                print('ERQ population:', count[ERQ_cluster_number])
                print('pops:', count)
                print('color:', iW3_mean)


        if(plotting==True):
                plot_db( labels, n_cl, sample, ERQ_cluster_number, str(e) + '-' + str(minPts) )

        if(ERQ_cluster_number==-1):
                print('No significant ERQ cluster found')
                return [0,0]
        else:
                return [ERQ, ERQ_scaled]



def ERQ_Shell(ERQ_scaled, dr, max_probe_distance, main_cluster_scaled, plotting, verbose):
        """
        ERQ thickness
        """
        import numpy as np
        import matplotlib.pyplot as plt



        dim=ERQ_scaled.shape[1]
        main_cluster_center=np.zeros([dim])

        for i in range(dim):
                main_cluster_center[i] = np.mean(main_cluster_scaled[:,i])

        r_low=0; r_high=0; dr=0.1
        nERQ=[]; rr=[]
        while(r_high<max_probe_distance):
                r_high = r_low + dr
                cc=0
                for i in range(len(ERQ_scaled)):
                        d=0.0
                        for k in range(0,dim):
                                d+=(main_cluster_center[k]-ERQ_scaled[i,k])**2
                        d=np.sqrt(d)

                        if ((r_low<d) and (d<r_high)):
                                cc+=1
                nERQ.append(cc)
                rr.append(r_low + dr*0.5)
                r_low+=dr


        for i in range(1,len(rr)-1):

                if((nERQ[i]>0) and (nERQ[i-1] == 0)):
                        r1=rr[i]

                if((nERQ[i]>0) and (nERQ[i+1] == 0)):
                        r2=rr[i]

        ERQ_thickness=r2-r1

        if(plotting==True):
                plt.plot(rr, nERQ , 'o', linewidth=3)
                plt.show()

        if (verbose==True):
                print('-----ERQ Shell info. -----')
                print('')
                print('r< ',r1)
                print('r> ', r2)
                print('ERQ shell thickness: ', ERQ_thickness)
        return [r1, r2, ERQ_thickness]




def Rad_P(pts, P, R_max, dr, verbose):

        """
        Probablity radius finder
        """

        import numpy as np


        dim = pts.shape[1]
        center=np.empty([dim])
        for i in range(0,dim):
                center[i] = np.mean(pts[:,i])
        r=0
        while(r<R_max):
                r+=dr
                count=0
                for i in range(0,len(pts)):
                        d=0
                        for k in range(0,dim):
                                d+=(center[k] - pts[i,k])**2
                        d=np.sqrt(d)
                        if(d<r):
                                count+=1
                if(float(count)/float(len(pts))>=P):
                        R=r
                        break

        if(verbose==True):
                print(" ------Probability Radius-----")
                print("")
                print('Number of points: ', len(pts))

                print('Number of points within R_ ', P,'=', count)
                print('R_', P, '= ', R)
        return R




def MeanShift(sample, plotting, dim, path, quantile, verbose, sampling_ratio):

        """
        ===================================
        mean-shift clustering
        ===================================
        """
        import numpy as np
        from sklearn.cluster import MeanShift, estimate_bandwidth
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from itertools import cycle


        # Scaling sample
        # sample_scaled = scale(sample)

        # bandwidth = estimate_bandwidth(sample_scaled, quantile=quantile, n_samples=int(sampling_ratio*len(sample)))
        bandwidth = estimate_bandwidth(sample, quantile=quantile, n_samples=int(sampling_ratio*len(sample)))
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(sample)
        # ms.fit(sample_scaled)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)

        count=np.zeros([n_clusters_])
        for i in range(0,n_clusters_):
                count[i]=0
                for j in range(0,len(labels)):
                        if(labels[j]==i):
                                count[i]+=1



        if(verbose==True):

                print('-----mean-shift------')
                print('pops.: ', count)
                print('centers: ', cluster_centers)
        if(plotting==True):

                colors = set_colors(labels)
                if(dim==3):
                        fg =plt.figure()
                        ax= Axes3D(fg)
                        for i in range(len(labels)):
                                ax.scatter(sample[i,0], sample[i,1], sample[i,2],  s=5, color=colors[i] )
                else:
                        for i in range(len(labels)):
                                plt.scatter(sample[i,0], sample[i,1],   s=5, color=colors[i] )

                plt.savefig(path)
                plt.show()
                plt.close()





def MeanShiftSeparation(sample, plotting, quantile, verbose, sampling_ratio):

        """
        ===================================
        mean-shift main cluster removal
        ===================================
        """
        import numpy as np
        from sklearn.cluster import MeanShift, estimate_bandwidth
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from itertools import cycle




        # Scaling sample
        # sample_scaled = scale(sample)

        # bandwidth = estimate_bandwidth(sample_scaled, quantile=quantile, n_samples=int(sampling_ratio*len(sample)))
        bandwidth = estimate_bandwidth(sample, quantile=quantile, n_samples=int(sampling_ratio*len(sample)))
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(sample)
        # ms.fit(sample_scaled)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)

        count=np.zeros([n_clusters_])
        for i in range(0,n_clusters_):
                count[i]=0
                for j in range(0,len(labels)):
                        if(labels[j]==i):
                                count[i]+=1


        # scaled main cluster and noise points
        dim = sample.shape[1]
        noise_pop = int(sum(count[1:]))
        # noise_scaled = np.zeros([noise_pop, dim])
        noise = np.zeros([noise_pop, dim])

        main_cluster=np.zeros([int(count[0]), dim])
        # main_cluster_scaled=np.zeros([int(count[0]), dim])
        i_noise=0
        i_main_cluster=0
        for i in range(len(labels)):
                if (labels[i]==0):
                        # main_cluster_scaled[i_main_cluster,:] = sample_scaled[i,:]
                        main_cluster[i_main_cluster,:] = sample[i,:]
                        i_main_cluster+=1
                else:
                        noise[i_noise,:] = sample[i, :]
                        # noise_scaled[i_noise,:] = sample_scaled[i,:]
                        i_noise+=1


        if(verbose==True):
                print("-----------Mean shift clustering----------")
                print("==========================================")
                print("Quantail= ", quantile)
                print('Estimated bandwith for sample:', bandwidth)
                print("number of estimated clusters : %d" % n_clusters_)
                print('-----------------clusters----------------------')
                for j in range(0, n_clusters_):
                        c=0
                        for i in range(0, len(labels)):
                                if(labels[i]==j):
                                        if(sample[i,0]>=4.6):
                                                c+=1
                        print(str(j+1) , ') ' , str(int(count[j])), '\t' , 'nERQ=', str(int(c)) , '\n')

        if(plotting==True):
                # path = input('fig name?' )
                path ='2d-rew-kt80'
                dim = 2
                # dim = input('2 for 2D, 3 for 3D')
                lx = 'iw3'
                # lx = input('xlabel?' )
                ly ='kt80'
                # ly = input('y label?' )
                # d1 = int(input('x axis?'))
                d1 = 0
                # d1 = int(input('x axis?'))
                d2 = 1
                # dim = X.shape[1]
                # d1=0; d2=1; d3=2

                if(dim==3):
                        d3 = int(input('z axis?' ))
                        lz= input('z label?' )
                        fg = plt.figure()
                        ax = Axes3D(fg)
                        ax.scatter(main_cluster[:,d1], main_cluster[:,d2], main_cluster[:,d3],  s=4, color='r', marker = 'x', alpha=1 )
                        ax.scatter(noise[:,d1], noise[:,d2], noise[:,d3],  s=4, color='b' )
                        ax.set_xlabel(lx)
                        ax.set_ylabel(ly)
                        ax.set_zlabel(lz)
                        fg.savefig(path + '.pdf')

                else:
                        plt.scatter(noise[:,d1], noise[:,d2],  s=4, color='r')
                        plt.scatter(main_cluster[:,d1], main_cluster[:,d2],  s=4, color='b', alpha=1, marker='x')
                        plt.xlabel(lx)
                        plt.ylabel(ly)
                        plt.savefig(path+'.pdf')

                plt.show()
                plt.close()

        # return [noise, noise_scaled, main_cluster, main_cluster_scaled]
        return 0



def qhull1(pts, verbose, plotting):

        """
        qhull1 gives Convex hull
        plots convex hull in 2D/3D
        """
        from scipy.spatial import ConvexHull
        import matplotlib.pyplot as plt


        dim = pts.shape[1]
        hull = ConvexHull(pts)
        V =hull.volume
        den = float(len(pts))/V
        if (plotting==True):
                if(dim==2):
                        plt.plot(pts[:,0], pts[:,1], 'o')
                        for simplex in hull.simplices:
                                plt.plot(pts[simplex, 0], pts[simplex, 1], 'k-')
                else:
                        plt.plot(pts[:,0], pts[:,1], pts[:,2], 'o')
                        for simplex in hull.simplices:
                                plt.plot(pts[simplex, 0], pts[simplex, 1], pts[simplex, 2],  'k-')


                plt.axis('equal')
        if(verbose==True):
                print('-------Convex hull info.--------')
                print('Vol. =', V)
                print('Density = ', den)
        return [V, den]



def random_rotation(p, center):
        """
        Rotates a set of points randomly
        """


        from scipy.stats import special_ortho_group
        import numpy as np
        dim=p.shape[1]
        T = special_ortho_group.rvs(dim)
        p_rot=np.empty([p.shape[0],p.shape[1]])

        p = np.array(p)
        center  = np.array(center)
        p = p - center
        for i in range(0, len(p)):

                p_rot[i,:]=np.matmul(T, p[i,:])
        return p_rot + center


def manual_rotation_2D(p, center, angle):
        """
        p , center, angle -> p-rot
        """
        import numpy as np

        T = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

        p = np.array(p)


        center = np.array(center)
        p = p -center
        p_rot=np.empty([p.shape[0],p.shape[1]])
        for i in range(0, len(p)):

                p_rot[i,:]=np.matmul(T, p[i,:])
        return p_rot + center


def manual_rot2D_ellipse1(pts, axis_ratio, rot_angel, center):
        """
        axis ratio is b/a
        """
        import numpy as np

        pts = np.array(pts) - np.array(center)
        n_pts = pts.shape[0]
        pts_rot = np.empty([n_pts,2])

        for i in range(n_pts):
                x = pts[i,0]; y = pts[i,1]

                # print(x,y)
                if ((x>=0.0) and (y>=0.0)):
                        theta = np.arctan(y/x)
                        # print('1st Quadrant ')

                # 2nd Quadrant
                if((x<=0.0) and (y>=0.0)):
                        theta = np.pi - np.arctan(-y/x)
                        # print('2nd Quadrant')

                # 3rd Quadrant
                if((x<=0.0) and (y<=0.0)):
                        theta = np.pi + np.arctan(y/x)
                        # print('3rd Quadrant')
                # 4th Quadrant
                if((x>=0.0) and (y<=0.0)):
                        theta = 2*np.pi - np.arctan(-y/x)
                        # print('4th Quadrant')
                # print(theta)
                theta = theta + rot_angel # rotated angel
                # print(theta)
                a = np.sqrt(x**2 + (y/axis_ratio)**2)
                b = a*axis_ratio
                f = np.tan(theta) # ratio of y_rotated to x_rotated
                # print(f)
                # if(f>=0):
                        # pts_rot[i,0] = 1.0/np.sqrt((f/b)**2 + 1.0/a**2)
                # else:
                pts_rot[i,0] = 1.0/np.sqrt((f/b)**2 + 1.0/a**2)

                pts_rot[i,1] = f*pts_rot[i,0]

        return pts_rot

def manual_rot2D_ellipse2(pts, axis_ratio, x_transfer, center):
        """
        axis ratio is b/a
        """
        import numpy as np
        import sys

        pts = np.array(pts) - np.array(center)
        n_pts = pts.shape[0]
        pts_rot_p = np.empty([n_pts,2])
        pts_rot_n = np.empty([n_pts,2])
        for i in range(n_pts):
                x = pts[i,0]; y = pts[i,1]
                a = np.sqrt(x**2 + (y/axis_ratio)**2)
                b = a*axis_ratio
                x=x+x_transfer
                if(abs(x)>a):
                        print('invalid transformation, decrees the x_transfer value ')
                        sys.exit()
                y=b*np.sqrt(1-(x/a)**2)
                pts_rot_p[i,0]=x
                pts_rot_n[i,0]=x
                pts_rot_p[i,1]=y
                pts_rot_n[i,1]=-y

        return [pts_rot_p, pts_rot_n]





def in_hull(hull, x):
        #  Checking which points reside inside the contour
        from scipy.optimize import linprog
        import numpy as np
        n_points = len(hull)
        n_dim = len(x)
        c = np.zeros(n_points)
        A = np.r_[hull.T,np.ones((1,n_points))]
        b = np.r_[x, np.ones(1)]
        lp = linprog(c, A_eq=A, b_eq=b)
        return lp.success




def cen(pts):
        import numpy as np
        dim = pts.shape[1]
        center=np.empty([dim])
        for i in range(0,dim):
                center[i] = np.mean(pts[:,i])
        return center



def poisson(density, num_pts):
        """
        Poisson probabolity:
        `density` is the mean density of
        the space we are gonna check overdensity
        `num_pts` is the number of points in the subarea
        residing in the region of interst

        """
        import numpy as np

        if(num_pts!=0):
                lnp = num_pts*np.log(density) - num_pts*np.log(num_pts) + num_pts -0.5*np.log(2*np.pi*num_pts) -density
                y = np.exp(lnp)
        else:
                y=np.exp(-density)
        return y




def optics_ini(path_sample, eps, minPts):
        import random
        from pyclustering.cluster import cluster_visualizer
        from pyclustering.cluster.optics import optics, ordering_analyser, ordering_visualizer
        from pyclustering.utils import read_sample, timedcall
        from pyclustering.samples.definitions import SIMPLE_SAMPLES, FCPS_SAMPLES
        from astropy.table import Table, Column
        import numpy as np

        def template_clustering(path_sample, eps, minpts, amount_clusters = None, visualize = True, ccore = True):
                sample = read_sample(path_sample)

                optics_instance = optics(sample, eps, minpts, amount_clusters, ccore)
                (ticks, _) = timedcall(optics_instance.process)

                print("Sample: ", path_sample, "\t\tExecution time: ", ticks, "\n")

                if (visualize is True):
                        clusters = optics_instance.get_clusters()
                        noise = optics_instance.get_noise()
                        visualizer = cluster_visualizer()
                        visualizer.append_clusters(clusters, sample)
                        visualizer.append_cluster(noise, sample, marker = 'x')
                        visualizer.show()
                        ordering = optics_instance.get_ordering()
                        analyser = ordering_analyser(ordering)
                        ordering_visualizer.show_ordering_diagram(analyser, amount_clusters)

        def cluster_reza():
                template_clustering( path_sample, eps, minPts)


        cluster_reza()


def LOF(sample,  contamination, n_neighbors, ):
        """
        Local Outlier detection, uising Mahalanobis metric which works
        based on covariance matrix. This metric  uses variance and correlation 
        based distances for data set so there is no need for normalization. 
        """

        from sklearn.neighbors import LocalOutlierFactor
        import numpy as np
        V = np.cov(sample.T)
        clf = LocalOutlierFactor(n_neighbors=n_neighbors, contamination=contamination,
        novelty= False, metric_params={'V': V}, metric='mahalanobis')
        y_pred = clf.fit_predict(sample)
        X_scores = clf.negative_outlier_factor_
        return  [y_pred, X_scores]


def modified_arctan(x,y):
    import numpy as np

    if(x>0 and y>=0):
        z=np.arctan(y/x)
    if(x>0 and y<=0):
        z = 2*np.pi - np.arctan(-y/x)
    if(x<0 and y<=0):
        z = np.pi + np.arctan(y/x)
    if(x<0 and y>=0):
        z = np.pi - np.arctan(-y/x)
    return z

def Line(P,m, x):
        """
        Line defined by a Point (P) and a slope (m)
        function returns y value at x 

        """
        return m*(x-P[0]) + P[1]

def LOF_Separater(sample, contamination, NN):
        """
        Seprating sample to ERQs, non ERQ outliers, and inliers
        """

        import numpy as np

        iW3 = sample[:,0]
        rew = sample[:,1]
        iw3_mean = np.mean(iW3)
        rew_mean = np.mean(rew)
        
        
        # Numpy.cov: Each row of m represents a variable, and each column a single observation of all those variables
        [y_pred, sample_score] = LOF(sample, contamination, NN)
        
        rew_inlier=[]
        iW3_inlier=[]
        iW3_ERQ=[]
        rew_ERQ=[]
        iW3_blue=[]
        rew_blue=[]
        iw3_outlier=[]; rew_outlier=[]
        score_out=[]; score_ERQ=[]; score_in=[]; score_blue=[]; blue_ind=[]
        outlier=0; inlier=0
        for i in range(len(sample)):
                if (y_pred[i]==1):
                        inlier+=1
                        rew_inlier.append(sample[i,1])
                        iW3_inlier.append(sample[i,0])
                        score_in.append(sample_score[i])
                        
                else:
                        outlier+=1
                        iw3_outlier.append(sample[i,0])
                        rew_outlier.append(sample[i,1])
                        score_out.append(sample_score[i])
                        if rew[i]>rew_mean and iW3[i]>iw3_mean:
                                rew_ERQ.append(sample[i,1])
                                iW3_ERQ.append(sample[i,0])
                                score_ERQ.append(sample_score[i])
                        if(rew[i]<rew_mean and iW3[i]<iw3_mean):
                                blue_ind.append(i)
                                rew_blue.append(sample[i,1])
                                iW3_blue.append(sample[i,0])
                                score_blue.append(sample_score[i])

        ERQ=np.array(list(zip(iW3_ERQ, rew_ERQ)))
        MCL = np.array(list(zip(iW3_inlier, rew_inlier)))
        outlier = np.array(list(zip(iw3_outlier, rew_outlier)))
        blue = np.array(list(zip(iW3_blue, rew_blue)))
        return ERQ, MCL, outlier, blue, blue_ind, sample_score, score_ERQ, score_in, score_out, score_blue


def wedge(ERQ,  MCL_center, enclosing_ratio):
        import numpy as np
        ERQ_center=np.empty([2])
        ERQ_center[0]= np.mean(ERQ[:,0])
        ERQ_center[1]= np.mean(ERQ[:,1])
        # --------- ERQ Vector --------------
        ERQ_vector = ERQ_center - MCL_center
        # calculating bin points
        theta = np.arctan(ERQ_vector[1]/ERQ_vector[0]) # ERQ vector angle
        A = np.linspace(0, np.pi/3.0, 500)
        for alpha in A:
                m_up = np.tan(theta + alpha)
                m_down = np.tan( theta - alpha)
                pop=0
                for i in range(len(ERQ)):
                        # upper line of the weddge
                        x,y = ERQ[i,:]
                        if (y< Line(MCL_center, m_up, x) and y > Line(MCL_center, m_down, x)):
                                pop+=1
                
                if(float(pop)/float(len(ERQ))>=enclosing_ratio):
                        return alpha, ERQ_center
                                

        



def binning( n_bin, sample, ERQ, MCL, enclosing_ratio, path):
        import numpy as np
        import matplotlib.pyplot as plt

        MCL_center = np.zeros([2])
        MCL_center[0] =np.mean(MCL[:,0])
        MCL_center[1] =np.mean(MCL[:,1])
        theta, ERQ_center=wedge(ERQ, MCL_center, enclosing_ratio)
        ERQ_vector = ERQ_center - MCL_center
        # calculating bin points
        m = ERQ_vector[1]/ERQ_vector[0]
        m_prep = -1.0/m

        
        # Wedge angle
        # Creating points which are on the bin line point
        P_bin=np.empty([n_bin+1,2])
        P_bin[0,0]=MCL_center[0]
        P_bin[0,1]=MCL_center[1]
        alpha = np.arctan(ERQ_vector[1]/ERQ_vector[0])
        End_x = np.max(ERQ[:,0])
        End_y = End_x * np.tan(alpha)
        d = np.sqrt((MCL_center[0] -End_x)**2 + (MCL_center[1]- End_y)**2)/float(n_bin)
        for i in range(1,n_bin+1):
                P_bin[i,0]=P_bin[i-1,0] + d*np.cos(alpha)
                P_bin[i,1]=P_bin[i-1,1] + d*np.sin(alpha)


        P_rot_low = manual_rotation_2D(P_bin, MCL_center, -theta )
        P_rot_up = manual_rotation_2D(P_bin, MCL_center, theta)

        m_up = np.tan(alpha + theta)
        m_down = np.tan(alpha -theta)
        for i in range(n_bin+1):
                x_up = (MCL_center[1] - P_bin[i, 1] + m_prep*P_bin[i,0] - m_up*MCL_center[0])/(m_prep - m_up)
                x_down = (MCL_center[1] - P_bin[i, 1] + m_prep*P_bin[i,0] - m_down*MCL_center[0])/(m_prep - m_down)
                y_up = Line([P_bin[i,0], P_bin[i,1]], m_prep, x_up)
                y_down = Line([P_bin[i,0], P_bin[i,1]], m_prep, x_down)
                plt.plot([x_up, x_down], [y_up, y_down], linewidth=1)
        plt.plot(P_rot_low[:,0], P_rot_low[:,1],  c='r', linewidth=1)
        plt.plot(P_rot_up[:,0], P_rot_up[:,1],  c='r', linewidth=1)
        # Plotting ERQ_vector o top of data  + bin points
        imw3=sample[:,0]; rew_gf=sample[:,1]
        plt.scatter(imw3, rew_gf, s=2, color='C3', alpha=0.5, label='Other ouliers')
        plt.scatter(MCL[:,0], MCL[:,1], s=5, color='C4', label='Main Cluster')
        plt.scatter(ERQ[:,0], ERQ[:,1] , s=7, color='C9', label='ERQs')
        plt.scatter(MCL_center[0], MCL_center[1], marker='x', s=50, color='g', label='MCL center')
        plt.scatter(ERQ_center[0], ERQ_center[1], marker='x', s=50, color='r', label='ERQ center')
        plt.legend()
        plt.savefig(path)
        plt.clf()
        # plt.show()


        # Determining points' bin
        bin_label=np.zeros([len(sample)])
        bin_pop=np.zeros([n_bin-1])
        for i in range(len(sample)):
                # upper line of the wedge
                x,y = sample[i,:]
                if (y< Line(MCL_center, m_up, x) and y > Line(MCL_center, m_down, x)):
                        for j in range(1,n_bin):
                                if(y> Line(P_bin[j-1,:], m_prep, x) and y < Line(P_bin[j,:], m_prep, x)):
                                        bin_label[i]=j
                                        bin_pop[j-1]+=1
        return bin_label, bin_pop


def boundary(sample, ERQ, MCL, out, y_range, score, score_ERQ, score_out, path):
        """
        Fits the the y intercept of the line
        which separates MCL from ERQ. No point
        remains to the left of the line from MCL. 
        """
        import numpy as np
        import matplotlib.pyplot as plt


        ERQ_center=np.zeros([2])
        ERQ_center[0] = np.median(ERQ[:,0])
        ERQ_center[1] = np.median(ERQ[:,1])

        MCL_center = np.zeros([2])
        MCL_center[0] =np.mean(MCL[:,0])
        MCL_center[1] =np.mean(MCL[:,1])


        ERQ_vec = ERQ_center - MCL_center
        slope = -ERQ_vec[0]/ERQ_vec[1]
        yi= y_range[0]; yf=y_range[1]
        y0 = np.linspace(yi,yf, 100)
        up = np.zeros([len(y0)])
        for j in range(len(y0)):

                for i in range(len(MCL)):
                        if Line([0,y0[j]],slope, MCL[i,0])< MCL[i,1]:
                                up[j]+=1

                for  i in range(0,len(y0)):
                        if(up[i]<1):
                                j1=i
                                break
        y0_best=y0[j1]
        y_line = np.linspace(MCL_center[1] - 0.5, ERQ_center[1] + 0.5, 10)
        x_line = (y_line - y0_best)/slope
        plt.clf()
        plt.plot(x_line, y_line)
                
        p=plt.scatter(sample[:,0], sample[:,1], c=score, s=0.5, cmap='plasma', alpha=0.7)
        plt.scatter(MCL_center[0], MCL_center[1], s=70, marker='x', c='black')
        plt.scatter(ERQ_center[0], ERQ_center[1], s=70, marker='x', c= 'black')
        plt.scatter(out[:,0], out[:,1], c=score_out,  marker ='v' ,cmap = 'plasma', label= 'none ERQ outlier', s=4)
        plt.scatter(ERQ[:,0], ERQ[:,1], marker='s', label = 'ERQ', c= score_ERQ, cmap='plasma', s=4)
        plt.xlabel(r'i - $W_3$')
        plt.ylabel(r'$\log(REW)$')
        plt.title('y = '+str(round(slope,3)) + 'x + ' + str(round(y0_best,3))   )
        # plt.axis('equal')
        plt.legend()
        plt.colorbar(p)
        plt.savefig(path)
        plt.close()
        return y0_best