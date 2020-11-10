import numpy
import scipy.stats

def binomial_confidence(n, k, q=0.6875):
    """
    Calculates binomial confidence intervals using the beta distribution.
    See e.g. Cameron 2011, PASA, 28, 128
    @param n: number of tries
    @type n: number/array
    @param k: number of successes
    @type k: number/array
    @param q: confidence interval (i.e. 68% would be a 1sigma confidence interval)
    @type q: float
    @return: confidence limit at q,p hat,confidence limit at 1-q
    @rtype: array
    """
    if type(n) is not list:
        n = [n]
    if type(k) is not list:
        k = [k]
    if len(n) != len(k):
        print "n and k must have same shape."
        raise ValueError
    full_n = numpy.array(n,dtype=float)
    full_k = numpy.array(k,dtype=float)
    conf_low = numpy.zeros_like(full_n)
    conf_high = numpy.zeros_like(full_n)
    p_hat = numpy.zeros_like(full_n)
    j = 0
    for n,k in zip(full_n, full_k):
        if k == 0:
            p_hat[j] = 0
            conf_low[j] = 0
            conf_high[j] = scipy.stats.beta.ppf(q, k+1, n-k+1)
        elif n == k:
            p_hat[j] = 1
            conf_low[j] = scipy.stats.beta.ppf((1.-q), k+1, n-k+1)
            conf_high[j] = 1
        else:
            conf_low[j] = scipy.stats.beta.ppf((1-q)/2., k+1, n-k+1)
            conf_high[j] = scipy.stats.beta.ppf(1-(1-q)/2., k+1, n-k+1)
            p_hat[j] = 1.*k/n
        j += 1
    return numpy.array([conf_low, p_hat, conf_high])