import astropy.io.fits as pyfits
from numpy import *
import scipy

def readDR10spec(input):
    """
    Load a 1D SDSS DR10 spectrum and extract wavelength solution.
    @param input: input file name
    @type input: string
    """
    dat = pyfits.open(input)
    wl = 10.0**(dat[1].data['loglam'])
    temperr = dat[1].data['ivar']
    air = mean(dat[2].data['airmass'])
    mask = temperr <= 0.0
    temperr[mask] = 1.0e-5
    err = 1./sqrt(temperr)
    return {'wl':wl,'flux':dat[1].data['flux'],'error':err,'z':dat[2].data['Z'],'air':air}

def readDR10specnoerr(input):
    """
    Load a 1D SDSS DR10 spectrum and extract wavelength solution.
    @param input: input file name
    @type input: string
    """
    dat = pyfits.open(input)
    wl = 10.0**(dat[1].data['loglam'])
#    temperr = dat[1].data['ivar']
    air = mean(dat[2].data['airmass'])
#    mask = temperr <= 0.0
#    temperr[mask] = 1.0e-5
#    err = 1./sqrt(temperr)
    return {'wl':wl,'flux':dat[1].data['flux'],'air':air}

def readDR7spec(input):
    """
    Load a 1D SDSS DR7 spectrum and extract wavelength solution.
    @param input: input file name
    @type input: string
    """
    dat = pyfits.open(input)
    w1 = dat[0].header['CRVAL1']
    wp = w1
    dw = dat[0].header['CD1_1']
    z = dat[0].header['z']
    wl = ones(n)
    for k in range(1,n):
       wp = wp + dw
       wl[k] = 10.0**wp
    return {'wl':wl,'wlz':wlz,'flux':dat[0].data[0],'error':dat[0].data[1],'z':z}
