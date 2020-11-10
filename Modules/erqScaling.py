"""
Tools for normalizing the data-set
"""

def MinMaxScaler(X):
	""" Scaling to min-max range for each column of 
	the data set"""
	import numpy as np
	m = np.min(X, axis=0)
	r = np.max(X, axis=0) - m
	X1=(X - m)/r
	
	return X1,m,r


def scale(x):
        """
        Scaling to mean and variance
        output has zero mean and variance of 1
        x is an N-dim array and output will be overwritten 
        on x
        """
        import numpy as np


        dim = x.shape[1]
        y = np.empty([len(x), dim])
        mean= np.empty([dim])
        std= np.empty([dim])
        for i in range(0,dim):
                mean[i] = np.mean(x[:,i])
                std[i] = np.std(x[:,i])
                y[:,i]=(x[:,i] - mean[i])/std[i]

        return (y, mean, std)

def scale_inv(x, mean, var):
        """
        Scaling back to mean and variance
        input has zero mean and variance of 1
        but the output gets back to its initial 
        status x is an N-dim array and output will 
        be overwritten on x
        """
        import numpy as np
        y=np.empty([len(x)])
        for i in range(len(y)):
                y[i]=x[i] *var + mean

        return y
