""" DBSCAN tools for ERQ clustering"""

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
        """
        Finds the main population of inliers
        in a data set X and plots the results
        """

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
	""" assigne color to labeled data"""

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



