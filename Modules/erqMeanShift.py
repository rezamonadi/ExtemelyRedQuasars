"""
Mean-Shitf algorith tools for ERQ clustering 
"""
def ERQ_Shell(ERQ_scaled, dr, max_probe_distance, main_cluster_scaled, plotting, verbose):
        """
        ERQ thickness finder.
	Given ERQ population and Main CLuster
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
        we didn't use it because it depends on quantile which is a distance
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
        rotates the hull in 2d around the center point
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
        """
        Checking which points reside inside the convex hull
        """
        from scipy.optimize import linprog
        import numpy as np
        n_points = len(hull)
        n_dim = len(x)
        c = np.zeros(n_points)
        A = np.r_[hull.T,np.ones((1,n_points))]
        b = np.r_[x, np.ones(1)]
        lp = linprog(c, A_eq=A, b_eq=b, method='revised_simplex')
        return lp.success
 = num_pts*np.log(density) - num_pts*np.log(num_pts) + num_pts -0.5*np.log(2*np.pi*num_pts) -density
                y = np.exp(lnp)
        else:
                y=np.exp(-density)
        return y

