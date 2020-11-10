import numpy as np 
import scipy.special as sc
import scipy.integrate as si



#  integrate pdf once to get CDF
#  now you can interpolate CDF so you'll get infinite resolution 
#  then 
def ConfInter(PDF, x, level):
     """
     ConfInter: Computes confidence interval of a PDF given the
     confidence level and data points.
     """
     A = si.simps(PDF, even='avg')
     PDF = PDF/A

     def CDF(PDF, ind_a):

          """
          returns CDF of a PDF at point a iwith index=ind_a
          """

          return si.simps(PDF[0:ind_a], even='avg')
     # print(CDF(PDF, len(PDF)))
     eps=1e-2
     # peak_range = np.array(np.where(abs(PDF-PDF.max())<eps))
     # print(peak_range)
     # peak_ind = int(peak_range.shape[1]*0.5) # position of maximum
     # max_ind = peak_range[0,peak_ind]
     d_min=1e5
     for i in range(1, len(PDF)):
          for j in range(1, len(PDF)):
               if(j<i):
                    # print(CDF(PDF, i) - CDF(PDF,j))
                    if(abs(CDF(PDF, i) - CDF(PDF,j) -level)<eps):
                         d=abs(x[i] - x[j])
                         # print('d', d)
                         if(d<d_min):
                              l1 = x[j]
                              l2 = x[i]
                              d_min = d
                              # print(iStart)

     return [l1, l2]


