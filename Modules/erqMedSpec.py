""" Tools for median Spectra """

def readDR10spec(input):
        """
        Load a 1D SDSS DR10 spectrum and extract wavelength solution.
        @param input: input file name
        @type input: string
        """
        import astropy.io.fits as pyfits
        import numpy as np
        dat = pyfits.open(input)
        wl = 10.0**(dat[1].data['loglam'])
        temperr = dat[1].data['ivar']
        air = np.mean(dat[2].data['airmass'])
        mask = temperr <= 0.0
        temperr[mask] = 1.0e-5
        err = 1./np.sqrt(temperr)
        return {'wl':wl,'flux':dat[1].data['flux'],'error':err,'z':dat[2].data['Z'],'air':air}

# def readDR10specnoerr(input):
#             """
#             Load a 1D SDSS DR10 spectrum and extract wavelength solution.
#             @param input: input file name
#             @type input: string
#             """
#             dat = pyfits.open(input)
#             wl = 10.0**(dat[1].data['loglam'])
#         #    temperr = dat[1].data['ivar']
#             air = mean(dat[2].data['airmass'])
#         #    mask = temperr <= 0.0
#         #    temperr[mask] = 1.0e-5
#         #    err = 1./sqrt(temperr)
#             return {'wl':wl,'flux':dat[1].data['flux'],'air':air}

# def readDR7spec(input):
#         """
#         Load a 1D SDSS DR7 spectrum and extract wavelength solution.
#         @param input: input file name
#         @type input: string
#         """
#         dat = pyfits.open(input)
#         w1 = dat[0].header['CRVAL1']
#         wp = w1
#         dw = dat[0].header['CD1_1']
#         z = dat[0].header['z']
#         wl = ones(n)
#         for k in range(1,n):
#             wp = wp + dw
#             wl[k] = 10.0**wp
#         return {'wl':wl,'wlz':wlz,'flux':dat[0].data[0],'error':dat[0].data[1],'z':z}


def cn_PnPoly(V, P):
        """
        check if a point is inside the polygon
        made by set of points or not"""

        cn = 0    # the crossing number counter

        # repeat the first vertex at end
        V = tuple(V[:])+(V[0],)

        # loop through all edges of the polygon
        for i in range(len(V)-1):   # edge from V[i] to V[i+1]
            if ((V[i][1] <= P[1] and V[i+1][1] > P[1])   # an upward crossing
                or (V[i][1] > P[1] and V[i+1][1] <= P[1])):  # a downward crossing
                # compute the actual edge-ray intersect x-coordinate
                vt = (P[1] - V[i][1]) / float(V[i+1][1] - V[i][1])
                if P[0] < V[i][0] + vt * (V[i+1][0] - V[i][0]): # P[0] < intersect
                    cn += 1  # a valid crossing of y=P[1] right of P[0]

        return cn % 2   # 0 if even (out), and 1 if odd (in)



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


def in_contour(contour, x):
        """
        checks a point is in the given contour or not
        """
        import numpy as np

        center_x = np.mean(contour[:,0])
        center_y = np.mean(contour[:,1])
        center = np.array([center_x,center_y])
        # alpah = modified_arctan(x[0],x[1])
        contour-= center
        d1 = np.sqrt(contour[:,0]**2 + contour[:,1]**2)
        x-=center
        d = np.sqrt(x[0]**2 + x[1]**2)
        for i in range(len(d1)):
                if d1[i]<d:
                        return False
                        break


        return True



def modified_arctan(x,y):
        """
        artctan just works for the first quadrant
        so I fixed the problem if the angle is more than 90 degrees
        """

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


def Line_y(P,m, y):

        """
        Line defined by a Point (P) and a slope (m)
        function returns x value at y
        """
        return  P[0] + (y - P[1])/m


def wedge(ERQ,  MCL_center, enclosing_ratio, resolution):

        """
        Finds the right opening angle that encloses the enclosing_ratio*len(ERQ)
        centered at MCL_center (center of inliers)

        theta -> the angle between the central vector of the wedge and x-axis
        """
        from scipy.stats import cumfreq
        import numpy as np


        ERQ_center=np.median(ERQ, axis=1)
        center_to_ERQ = ERQ_center - MCL_center
        # --------- ERQ Vector --------------
        center_to_ERQ_unit = center_to_ERQ/ np.linalg.norm(center_to_ERQ)
        ERQ_vectors = sample_uniter(ERQ - MCL_center)
        ERQ_angles = np.dot(ERQ_vectors, ERQ_center)
        ERQ_angles = np.rad2deg(np.arccos(ERQ_angles))



        res = cumfreq(ERQ_angles, numbins=resolution)
        x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size,
                                 res.cumcount.size)
        count = res.cumcount
        # return np.rad2deg(np.arccos(x[np.where(count==int(enclosing_ratio*len(ERQ_dir)))]))
        angle = x[np.where(count>=int(enclosing_ratio*len(ERQ_angles)))]
        return angle[0]


def boundary(sample, ERQ, MCL, out, y_range, score, score_ERQ, score_out, path):
        """
        Fits the the y intercept of the line
        which separates MCL from ERQ. No point
        remains to the left of the line from MCL.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        SMALL_SIZE = 15
        MEDIUM_SIZE = 18
        BIGGER_SIZE = 18

        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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
        plt.plot(x_line, y_line, label='Boundary')

        p=plt.scatter(sample[:,0], sample[:,1], c=score, s=0.5, cmap='plasma', marker='.', alpha=0.5, label='Inlier')
        # plt.scatter(MCL_center[0], MCL_center[1], s=70, marker='x', c='black')
        # plt.scatter(ERQ_center[0], ERQ_center[1], s=70, marker='x', c= 'black')
        plt.scatter(out[:,0], out[:,1], c=score_out,  marker ='s' ,cmap = 'plasma', label= 'Outlier', s=4, alpha=1)
        # plt.scatter(ERQ[:,0], ERQ[:,1], marker='s', label = 'ERQ', c= score_ERQ, cmap='plasma', s=4, alpha=)
        plt.axvline(x=4.6,ymin=0.65, ls='--')
        plt.axhline(y=2, xmin=0.5, ls='--')
        plt.xlabel(r'i - $W_3$')
        plt.ylabel(r'$\log(REW)$')
        # plt.title('y = '+str(round(slope,3)) + 'x + ' + str(round(y0_best,3))   )
        plt.title('LOF scoring and ERQ boundary')
        # plt.axis('equal')
        plt.legend()
        cbar=plt.colorbar(p)
        cbar.set_label(r'${{\leftarrow}}  \ Outliers \ \ \ \ \ \ \  \ \ \ \  Inliers {{\rightarrow}}$')

        plt.savefig(path, bbox_inches='tight', quality=100, dpi=1200, format='jpg')
        plt.close()
        return y0_best


def ContourExpantion(Contour, expansion):

        import numpy as np

        center = np.mean(Contour, axis=0)
        Contour_0 = Contour - center
        exp_Contour = np.empty([len(Contour), 2])
        for i in range(len(Contour)):
                r=np.sqrt(Contour_0[i,0]**2+Contour_0[i,1]**2)
                theta=modified_arctan(Contour_0[i,0], Contour_0[i,1])
                exp_Contour[i,0] = np.cos(theta)*r*expansion
                exp_Contour[i,1] = np.sin(theta)*r*expansion

        return exp_Contour + center


def in_triangle(A, B, C, P, eps):
        """
        This checks whether a points
        is located inside a triangle or not

        Inputs:
        tri -> coordinates of the edges
        P    -> test point coordinate
        """

        import numpy as np

        A = np.array(A)
        B = np.array(B)
        C = np.array(C)
        P = np.array(P)

        # --------------------

        S = np.linalg.norm(np.cross(A-B, A-C))
        S1 = np.linalg.norm(np.cross(P -A, P-B))
        S2 = np.linalg.norm(np.cross(P -A, P-C))
        S3 = np.linalg.norm(np.cross(P -C, P-B))


        if (abs(S1+S2+S3-S)<=eps):
                return True
        else:
                return False


def edger(sample, center, theta, opening_angle):
        """
        This function finds the edges of a triangle
        needed for the wedge centered at main population center.
        This edges should make a triangle which encompasses all of the
        points in the sample and leave no point outside of the triangle.

        We find a maximum distance "r" from the center
        then rotate this rad according to the central vec and opening angle...

        """

        import numpy as np

        rmax=-100
        sample0=sample- center

        for i in range(len(sample)):
                r = np.linalg.norm(sample0[i,:])
                if(r>=rmax): rmax=r

        alpha_p = theta + opening_angle*0.5
        alpha_n = theta - opening_angle*0.5

        rmax*=1.02
        A = [rmax*np.cos(alpha_p) + center[0], rmax*np.sin(alpha_p) + center[1]]
        B = [rmax*np.cos(alpha_n) + center[0], rmax*np.sin(alpha_n) + center[1]]
        return A, B


def stacker(z_dr12, plate, mjd, fiberid):
        """
        Median spectrum for a given set of spectra.
        nanmedian considers nan points
        flux is normalized to somewhere in continuim
        """
        
        from tqdm import tqdm
        import numpy as np
        import scipy
        ## Define a log wavelength grid for the composite spectrum
        step = 1.00015
        bb = np.arange(0,8813,1)
        wgrid = 800.0 * step**bb
        nw = len(bb)

        nqsos = len(z_dr12)
        #  initialize the spectrum with zeros
        sp = np.zeros([int(nqsos), nw])
        
        for i in tqdm(range(nqsos)):
        # Retrieve the spectra:
                file = '/media/reza/My Passport/erq/fred/sdss/%d/spec-%d-%d-%04d.fits' % (plate[i], plate[i],mjd[i],fiberid[i])
                spec = readDR10spec(file)
                wave = spec['wl']
                wz = wave/(z_dr12[i]+1)
                flux = spec['flux']
                mask = (wz > 1680.0) & (wz < 1730.0)
                fnorm = np.median(flux[mask])
                fluxn = flux/fnorm
        # interpolate the rest-frame spectrum onto the standard grid
                f = scipy.interpolate.interp1d(wz,fluxn,bounds_error=False,fill_value=float('nan'))
                sp[i] = f(wgrid)
        # calculate the median spectrum
        med1 = np.nanmedian(sp,axis=0)
        return med1

def KDE_Bin2D(sample, rangeData, minData, x_erq, y_erq,ngrid, bw,levels,A, B, \
 ERQ_color, ERQ_rew, expansion_handle, expansion, path, tip, x_label, y_label, tit):
        """
        This function gives the labels for each bin
        given the opening angle and direction of the wedge centeral vector"""
        
        import numpy as np
        from scipy.stats import kde
        from tqdm import tqdm
        import sys
        import matplotlib.pyplot as plt
        from  sklearn.neighbors import KernelDensity
        import matplotlib.ticker as ticker
        print('KDE density estimation in 2D...')
        data =sample
        x, y = data.T
        # xi, yi, zi =kde2D(x,y, bw)
        k = kde.gaussian_kde(data.T)
        k.set_bandwidth(bw_method=k.factor*bw)
        print('Contour plot...')
        xi, yi = np.mgrid[x.min():x.max():ngrid*1j, y.min():y.max():ngrid*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        zi /=max(zi)
        density = np.array(k(sample.T))
        density/=max(density)
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        fig=plt.figure()
        plt.cla()
        ax = fig.add_subplot(111)
        # plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        # plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        # plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        # plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        # plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        # plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        # plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # fg, axes = plt.subplots()
        c = ax.contour(xi, yi, zi.reshape(xi.shape), levels=levels, alpha=0.6, c='black')
        ax.clabel(c, fontsize=5)
        segments = c.allsegs

        # finding the contours
        lines = []
        nContours = len(levels)
        for i in range(nContours):

                l = np.size(np.array(segments[i]))
                line = np.reshape(segments[i], (int(l/2),2))
                lines.append(line) # contours are added to lines regardless of their length

        if expansion_handle==True:
                for i in range(len(expansion)):

                        line_exp = ContourExpantion(lines[i], expansion[i])
                        lines.insert(0, line_exp)


        lines=np.array(lines)
        nContours = len(lines)
        # for i in range(nContours):

        center=np.empty([2])
        center[0] = np.median(sample[:,0])
        center[1] = np.median(sample[:,1])

        bin_label=np.zeros([len(sample)])
        nBins = nContours +1
        bin_pop=np.zeros([nBins])
        bin_med_position=np.zeros([nBins,2])
        tip_label=np.zeros([len(sample)])
        tip_pop=0



        #A, B = edger(sample, center, wedge_direction, opening_angle)
        color=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']
        print('hull...')
        sample0=sample-center
        r_max_plotting=-100
        ind_in_triangle=[]
        for i in tqdm(range(len(sample))):

                if(tip==True):
                        # if(in_hull(lines[-1], sample[i,:])==True):
                        if(cn_PnPoly(lines[nContours-1], sample[i,:])==1):
                                tip_pop+=1
                                tip_label[i]=1

                if (in_triangle(A,B, center, sample[i,:], 1e-8)==True):
                        ind_in_triangle.append(i)
                        # r=np.sqrt(sample0[i,0]**2+sample0[i,1]**2)
                        # if(r>r_max_plotting): r_max_plotting=r
                        for j in range(0, nBins):
                                Contour_ind = (nContours-1) -j


                                if(j==0): # the inner most bin

                                        if(cn_PnPoly(lines[Contour_ind], sample[i,:])==1):
                                                bin_label[i]=j+1
                                                bin_pop[j] +=1


                                if(0<j<nBins-1): # intermediate contours

                                        if( (cn_PnPoly(lines[Contour_ind], sample[i,:]) ==  1) and (cn_PnPoly(lines[Contour_ind+1], sample[i,:])==0) ):
                                                bin_label[i]=j+1
                                                bin_pop[j] +=1

                                if(j==nBins-1): # the bin out of expanded contour
                                        if(cn_PnPoly(lines[0],sample[i,:])==0):
                                                bin_label[i]=j+1
                                                bin_pop[j] +=1
        print('plotting')
        # plt.clf()
        # for i in tqdm(range(len(sample))):
        #         for j in range(0,nBins):
        #                 plt.scatter(sample[i,0], sample[i,1], c=color[j])


        sample_in_traingle = sample[ind_in_triangle]
        plt.plot([center[0], A[0] ], [center[1], A[1]], ls='-', c='r', alpha=0.5, lw=1)
        plt.plot([center[0], B[0]  ], [center[1], B[1]], ls='-', c='r', alpha=0.5, lw=1)
        #plt.plot([A[0], B[0] ], [A[1], B[1]], c='r', alpha=0.5)
        # plt.plot([center[0], (A[0]+B[0]-2*center[0])*0.4], [center[1], (A[1]+B[1]-2*center[1])*0.4], lw=2)
        plt.plot([x_erq,max(sample[:,0])], [y_erq,y_erq], ls='--', c='black', lw=1)
        plt.plot([x_erq,x_erq], [y_erq, max(sample[:,1])], ls='--', c='black', lw=1)
        plt.scatter(sample[:,0], sample[:,1], s=0.1)
        # plt.xlim(-0.1, 1.1)
        # plt.ylim(-0.1, 1.1)
        for i in range(0, nContours):
                l = lines[i]

                if (i<len(expansion)):
                        plt.plot(l[:,0], l[:,1], c = 'black', ls='--', alpha=0.5, lw=1)
                # else:
                        # plt.plot(l[:,0], l[:,1], c = 'black', ls='-', alpha=0.5, lw=1)


    
        
        # # # plb.axis('equal')
        ticks_x = ticker.FuncFormatter(lambda x, 
                                   pos: '{0:g}'.format(round(x*rangeData[0]+ minData[0], 2)))
        ax.xaxis.set_major_formatter(ticks_x)

        ticks_y = ticker.FuncFormatter(lambda x, 
                                  pos: '{0:g}'.format(round(x*rangeData[1]+minData[1],2)))
        ax.yaxis.set_major_formatter(ticks_y)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(tit)
        # r = np.max(sample, axis=0) - np.min(sample, axis=0)
        # plt.xlim(min(sample[:,0]) - 0.1*r[0], max(sample[:,0]) +0.1*r[0])
        # plt.ylim(min(sample[:,1]) - 0.1*r[1], max(sample[:,1]) + 0.1*r[1])
        for j in tqdm(range(nBins)):  
                if(j>0):
                        x=[]; y=[]
                        for i in range(len(sample)):
                                if (bin_label[i]==j+1):
                                        x.append(sample[i,0])
                                        y.append(sample[i,1])

                        bin_med_position0 = np.median(x)
                        bin_med_position1 = np.median(y)
                        plt.text(bin_med_position0, bin_med_position1, str(j), fontsize=8, color='red')

                        ax.text(.77, 0.35-j*0.04, 'Bin-'+ str(j)+ ' : #'+str(int(bin_pop[j])),fontsize=7, color='black' )
                else:
                        ax.text(.77, .35, 'Bin-C' + ' : #'+str(tip_pop),fontsize=7, color='black' )
                        ax.text(np.median(sample[:,0]), np.median(sample[:,1]), 'C', fontsize=8, color='red')
        # if(expansion_handle==True): plt.title(str(round(expansion[0],1))+', '+ str(round(expansion[1],1))+ ', ' + str(round(expansion[2],1)))
        # plt.axis('equal')
        
        plt.savefig(path, bbox_inches='tight', format='png', dpi=200)
        # plt.show()
        plt.close()

        print('tip_pop', tip_pop )
        if(tip==True):
                return bin_label, bin_pop, tip_label
        else:
                return bin_label, bin_pop


def vectors_uniter(sampleVector):
        import numpy as np
        sampleVector_unit = np.zeros(np.shape(sampleVector))
        dim = np.shape(sampleVector)
        if(np.size(dim)>1):
                L = np.linalg.norm(sampleVector,axis=1)
                for i in range(len(sampleVector)):

                        sampleVector_unit[i,:] = sampleVector[i,:]/L[i]
        else:
                sampleVector_unit = sampleVector/np.linalg.norm(sampleVector)
        return sampleVector_unit


def opening_angle_finder(ERQ, Main_center,  enclosing_ratio, resolution):
        import numpy as np
        from scipy.stats import cumfreq
        ERQ_center = np.median(ERQ, axis=0)
        ERQ_vector = np.array(ERQ_center) - Main_center
        ERQ_vector_unit = vectors_uniter(ERQ_vector)
        ERQ = ERQ - Main_center
        ERQ_unit = vectors_uniter(ERQ)

        wedge_direction = modified_arctan(ERQ_vector[0], ERQ_vector[1])
        ERQ_angles=[]
        for i in range(len(ERQ)):
                ERQ_angles.append(np.arccos(np.dot(ERQ_vector_unit, ERQ_unit[i,:])))



        res = cumfreq(ERQ_angles, numbins=resolution)
        x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size,
                                        res.cumcount.size)
        count = res.cumcount
        # return np.rad2deg(np.arccos(x[np.where(count==int(enclosing_ratio*len(ERQ_dir)))]))
        angle = x[np.where(count>=int(enclosing_ratio*len(ERQ_angles)))]
        return angle[0] , wedge_direction



def opening_angle_finder_nd(ERQ_angles,  enclosing_ratio, resolution):
    from scipy.stats import cumfreq

    res = cumfreq(ERQ_angles, numbins=resolution)
    x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size,
                                 res.cumcount.size)
    count = res.cumcount
    # return np.rad2deg(np.arccos(x[np.where(count==int(enclosing_ratio*len(ERQ_dir)))]))
    angle = x[np.where(count>=int(enclosing_ratio*len(ERQ_angles)))]
    return angle[0]


def uniter(sampleVector):
        import numpy as np
        sampleVector_unit = np.zeros(np.shape(sampleVector))
        L = np.linalg.norm(sampleVector,axis=1)
        for i in range(len(sampleVector)):

                sampleVector_unit[i,:] = sampleVector[i,:]/L[i]
        return sampleVector_unit

def separation(X, vec):
        """
        returns the cos(dot(X,vec)  in nD
        if X is a set of points
        """
        from numpy.linalg import norm
        import numpy as np
        X_unit = uniter(X)
        vec_unit = vec/norm(vec)

        # return np.rad2deg(np.arccos(np.dot(X_unit, vec_unit)))
        return (np.dot(X_unit, vec_unit))



def rectangle_corner_finder(opening_distance, rectangle_length, MainCenter, rectangle_vector):
        import numpy as np
        """
        finds the corners of the rectangle given the rectangle length
        and direction of the length """

        U = rectangle_vector/np.linalg.norm(rectangle_vector)

        V1 = np.array([-U[1], U[0]])  # unit vector perpendicular to rectangle_vector
        V2 = np.array([U[1], -U[0]])  # unit vector perpendicular to rectangle vector opposite to V1

        A = V2*opening_distance + MainCenter
        B = V1*opening_distance + MainCenter
        C = A+ rectangle_length*U
        D = B+ rectangle_length*U
        return A, B, C, D


def rectangle_corner_finder(opening_distance, rectangle_length, MainCenter, rectangle_vector):
        import numpy as np
        """
        finds the corners of the rectangle given the rectangle length
        and direction of the length """

        U = rectangle_vector/np.linalg.norm(rectangle_vector)

        V1 = np.array([-U[1], U[0]])  # unit vector perpendicular to rectangle_vector
        V2 = np.array([U[1], -U[0]])  # unit vector perpendicular to rectangle vector opposite to V1

        A = V2*opening_distance + MainCenter
        B = V1*opening_distance + MainCenter
        C = A+ rectangle_length*U
        D = B+ rectangle_length*U
        return A, B, C, D

def radius_finder(CERQ, Main_center,  enclosing_ratio, resolution):

        import numpy as np
        from scipy.stats import cumfreq

        CERQ_center = np.median(CERQ, axis=0)
        CERQ_vector = np.array(CERQ_center) - Main_center
        # ERQ_vector_unit = vectors_uniter(ERQ_vector)
        CERQ1 = CERQ - Main_center
        # ERQ_unit = vectors_uniter(ERQ)
        CERQ_r = np.linalg.norm(CERQ1, axis=1) # distance of each ERQ point from the MCL
        maxR_CERQ = np.max(CERQ_r)
        CERQ_theta = np.arccos(separation(CERQ1, CERQ_vector))
        wedge_direction = modified_arctan(CERQ_vector[0], CERQ_vector[1])
#         ERQ distances
        CERQ_unit= uniter(CERQ1)
        CERQ_vector_unit=CERQ_vector/np.linalg.norm(CERQ_vector)

        CERQ_theta = np.arccos(dot(CERQ_unit, CERQ_vector_unit)) # angle of each point from the line
        CERQ_distances= CERQ_r*np.sin(CERQ_theta)

        res = cumfreq(CERQ_distances, numbins=resolution)
        x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size,
                                        res.cumcount.size)
        count = res.cumcount
        # return np.rad2deg(np.arccos(x[np.where(count==int(enclosing_ratio*len(ERQ_dir)))]))
        distance= x[np.where(count>=int(enclosing_ratio*len(CERQ_distances)))]
        return distance[0], wedge_direction, maxR_CERQ

