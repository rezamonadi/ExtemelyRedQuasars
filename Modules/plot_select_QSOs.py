import pyfits
import numpy
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import scipy.ndimage
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from readSDSSspectrafast import *

# some line wavelengths
s6 = array([933.862,944.523])
c3 = array([977.020,1175.71,1908.73])
n3 = array([989.799, 1750.0])
lyb = 1025.722
o6 = array([1031.926,1037.617])
s4 = array([1062.664,1073.518])
n2 = array([1083.993,1084.580,1085.701])
p5 = array([1117.977,1128.008])
si3 = array([1206.500,1892.03])
lya = 1215.670
he2 = 1640.4
mg2 = array([2796.35, 2803.53])
n4 = 1486.5
o4 = array([1397.21,1399.78,1404.79,1407.39])
n5 = array([1238.821,1242.804])
si2a = array([1190.416,1193.290,1194.500])
si2b = array([1260.422,1264.738,1526.7066,1533.4310])
#o1 = array([1302.1685,1305.0])
o1 = 1304.
#c2 = array([1334.532,1335.708])
c2 = 1335.
si4 = array([1393.755,1402.770])
c4 = array([1548.195,1550.770])
fe2 = array([1608.4511, 1786.7])
fe3 = array([1122.52,1124.87,1128.72,1131.908,1895.46,1914.06,1926.30])
o3 = 1664.
al2 = array([1670.7874])
al3 = array([1854.716,1862.7895])

c = 2.9979e5

# *approx* ratio of NV to CIV wavelengths
ncrat = 0.8005
# amount of additive increase to NV b values in angstroms = difference between NV and CIV doublet separations
ncb = 1.9
# use this for plotting NV
def func(x, scle, wing, A0, x0, b0, A1, x1):
   b1 = wing*b0
   b0n = b0*ncrat+ncb
   b1n = b1*ncrat+ncb
   x0n = x0*ncrat
   x1n = x1*ncrat
######## This next line is a fudge
   if scle == 1.0: scle = 0.0
   return 1.0 + abs(scle*A0)*exp(-(((x-x0n)/b0n)**2.)) + abs(scle*A1)*exp(-(((x-x1n)/b1n)**2.)) + abs(A0)*exp(-(((x-x0)/b0)**2.)) + abs(A1)*exp(-(((x-x1)/b1)**2.))

# output pdf file name:
pdf = PdfPages('/home/reza/erq/fred/plot_select_QSOs.pdf')

# input fits file
thdu = pyfits.open('/home/reza/erq/fred/C4N5REWs_DR12v11.fits')
tdata = thdu[1].data
print(shape(tdata))
name = tdata.field('sdss_name')
plate = tdata.field('Plate')
mjd = tdata.field('MJD')
fiberid = tdata.field('FiberID')
nqsos = len(plate)

dec = zeros(nqsos)
ra = zeros(nqsos)
for k in range(nqsos):
   dec[k] = name[k][9:12]
   ra[k] = name[k][0:2]

zdr12 = tdata.field('z_dr12')  # This is the 'best' redshift available: PCA or z_vi if PCA not good
zp1 = zdr12 + 1.0
bal_vi = tdata.field('bal_flag_vi')     # this is the BAL visual inspection flag
bi_civ = tdata.field('bi_civ')
err_bi_civ = tdata.field('err_bi_civ')
qflag = tdata.field('qflag')
fwhm = tdata.field('fwhm_gf')
fwhm_err = tdata.field('fwhme_gf')
fwhm_snr = zeros(nqsos)
mask = fwhm_err > 0.0
fwhm_snr[mask] = fwhm[mask]/fwhm_err[mask]
fwhm_gfc = tdata.field('fwhm_gfc')
fwhm_gfw = tdata.field('fwhm_gfw')
shift_gf = tdata.field('shift_gf')
rat = tdata.field('rat_gf')
sigma = tdata.field('sigma_gf')
kurt75 = tdata.field('kurt75_gf')
kurt80 = tdata.field('kurt80_gf')
asy = tdata.field('asy_gf')
rew = tdata.field('rew_gf')
rew_gfc = tdata.field('rew_gfc')
rew_gfw = tdata.field('rew_gfw')
rew_err = tdata.field('rewe_gf')
rew_snr = tdata.field('rewsnr_gf')
peak = tdata.field('peak_gf')
peak_snr = tdata.field('peaksnr_gf')
wcore = tdata.field('wcore_gf')
rew_nv = tdata.field('rew_nv')
frat_nvciv = tdata.field('frat_nv/civ')
scale = tdata.field('scale')
f1450 = tdata.field('f1450')
alpha_nv = tdata.field('alpha_nv')
alpha_civ = tdata.field('alpha_civ')
alpha_all = tdata.field('alpha_all')
alpha_allc = tdata.field('alpha_allc')
ABr = tdata.field('ABr')
ABi = tdata.field('ABi')
cc_flags = tdata.field('cc_flags')
ABw1 = tdata.field('ABw1')
ABw1_snr = tdata.field('ABw1_snr')
ABw2 = tdata.field('ABw2')
ABw2_snr = tdata.field('ABw2_snr')
ABw3 = tdata.field('ABw3')
ABw3_snr = tdata.field('ABw3_snr')
ABw4 = tdata.field('ABw4')
ABw4_snr = tdata.field('ABw4_snr')
rmz = tdata.field('r-z')
imw4 = tdata.field('i-w4')
rmw4 = tdata.field('r-w4')
imw3 = tdata.field('i-w3')
imw1 = tdata.field('i-w1')
w3mw4 = tdata.field('w3-w4')
snr_spec = tdata.field('snr_spec')
snr_1700 = tdata.field('snr_1700')

w1mw2 = ABw1 - ABw2
w2mw3 = ABw2 - ABw3
imw2 = ABi - ABw2

mask = (ABw1 < 0.0) | (ABw1_snr < 3.0) | (cc_flags != '0000')
ABw1[mask] = -99.0
w1mw2[mask] = -99.0

mask = (ABw2 < 0.0) | (ABw2_snr < 3.0) | (cc_flags != '0000')
ABw2[mask] = -99.0
w1mw2[mask] = -99.0
w2mw3[mask] = -99.0
imw2[mask] = -99.0

mask = (ABw3 < 0.0) | (ABw3_snr < 3.0) | (cc_flags != '0000')
ABw3[mask] = -99.0
w2mw3[mask] = -99.0
w3mw4[mask] = -99.0
imw3[mask] = -99.0

mask = (ABw4 < 0.0) | (ABw4_snr < 2.5) | (cc_flags != '0000')
ABw4[mask] = -99.0
w3mw4[mask] = -99.0
imw4[mask] = -99.0
rmw4[mask] = -99.0

# read in a list of quasar names from an ascii table
qname = loadtxt("/home/reza/erq/fred/plot_select_QSOs.txt", dtype=str)

wind = -1
count = 0
nerqs = 0
ncore = 0
nbcore = 0
nT2core = 0
nbals = 0
nbalsall = 0
nT2 = 0
nT2all = 0
ind = []
redshift = []
colorcore = []
colorerqs = []
colorall = []
color = []
nqsos = len(zdr12)
for j in range(nqsos):
   if name[j] not in qname: continue
   count = count + 1
   wind = wind + 1
# Sample statistics for all selected, then ERQs, then core ERQs:
#   if imw3[j] > 1.0: print imw3[j]
   if bal_vi[j] > 0: nbalsall = nbalsall + 1
   elif fwhm[j] < 2000:
         nT2all = nT2all+1
   if (zdr12[j] >= 2.0) and (zdr12[j] <= 3.4) and imw3[j] >= 4.6 and rew_snr[j] >= 4.0 and fwhm_snr[j] >= 4.0 and qflag[j] == 0:
      nerqs = nerqs+1
      if imw3[j] > 0.0: colorerqs.append(imw3[j])
      if bal_vi[j] > 0:
         nbals = nbals+1
      elif fwhm[j] < 2000:
         nT2 = nT2+1
      if rew[j] >= 100.0:
         ncore = ncore+1
         if imw3[j] > 0.0: colorcore.append(imw3[j])
         if bal_vi[j] > 0:
            nbcore = nbcore+1
         elif fwhm[j] < 2000:
            nT2core = nT2core+1
   ind.append(j)
   if imw3[j] > 0.0: colorall.append(imw3[j])
   redshift.append(zdr12[j])
   file = '/home/reza/erq/fred/sdss/spec-%d-%d-%04d.fits' % (plate[j],mjd[j],fiberid[j])
   file2 = '/home/reza/erq/fred/sdss/%d/spec-%d-%d-%04d.fits' % (plate[j],plate[j],mjd[j],fiberid[j])
#   print count,j,'retrieving file / quasar:', file, name[j]
# read the spectrum file
   spec = readDR10spec(file2)
   wz = spec['wl']/zp1[j]
   flux = spec['flux']
   fluxe = spec['error']
   sclean = flux.copy()
   scleane = fluxe.copy()
# calculate the continua define for CIV and NV
   cont = f1450[j] * (wz/1450.0)**alpha_civ[j]
   cont_tot = cont.copy()
   if rew_nv[j] > 0.0:
      mask = wz < 1450.0
      wzn = wz[mask]
      contnv = f1450[j] * (wzn/1450.0)**alpha_nv[j]
      cont_tot[mask] = contnv
   if alpha_all[j] > -5.0:
      mask = wz > 1600.0
      cont_all = f1450[j] * (wz/1450.0)**alpha_all[j]
      cont_tot[mask] = cont_all[mask]
   
# clean CRs off the flux array
   for k in range(10,len(wz)-10,1):
      mask = (wz > wz[k]-2.5) & (wz < wz[k]+2.5)
      medf = median(flux[mask])
      medfe = median(fluxe[mask])
      if abs(flux[k]-medf) > 5.0*medfe: sclean[k] = medf
      if abs(fluxe[k]-medfe)/medfe > 4.0: scleane[k] = medfe

# Calculate the profile fits if REW(CIV) > 0.0
   if rew[j] > 0.0:
      x0 = wcore[j]
      x1 = x0 + shift_gf[j]
      wing = fwhm_gfw[j]/fwhm_gfc[j]
      b0 = 1549.05*fwhm_gfc[j]/(c*1.665)
      A0 = rew_gfc[j]/(1.772*b0)
      if wing < 0.01:
         A1 = 0.0
         wing = 0.01
      else:
         A1 = rew_gfw[j]/(1.772*b0*wing)
      lfits = cont_tot*func(wz, scale[j], wing, A0, x0, b0, A1, x1)

#smooth and find max flux within the plot wavelength range
#   sm_raw = scipy.ndimage.filters.gaussian_filter1d(flux,1.5)
   sm_spec = scipy.ndimage.filters.gaussian_filter1d(sclean,1.2)
   sm_spec5 = scipy.ndimage.filters.gaussian_filter1d(sclean,1.7)
   sm_spece = scipy.ndimage.filters.gaussian_filter1d(scleane,1.2)
   wmin = max(1180.0,wz[20])
   mask = logical_and(wz > wmin, wz < 1600)
   ymax = max(sm_spec[mask].max(),lfits[mask].max()) * 1.14
   if wz[0] > 1250.0: ymax = ymax*1.12
##### The next line set the scale to look for BALs. Comment out otherwise:
#   ymax = 3.5*f1450[j]
#   mask = (wz > 2500) & (wz < 2830)
#   ymax = max(ymax,sm_spec[mask].max())
#   mask = logical_and(wz > 1720, wz < 1820)
   med = median(sm_spec[mask])
   if wind == 0:
      fig = plt.figure(figsize=(8.0,12.0))
      gs = gridspec.GridSpec(5,1)
      gs.update(left=0.12, right=0.90, top=0.98, bottom=0.08, hspace=0.0)

   print('J%s  %.2f  %.1f   %.1f   %.0f+/-%.0f   %.0f+/-%.0f   %.1f   %.1f   %.2f   %.2f   %d' 
   % (name[j],zdr12[j],ABi[j],imw3[j],rew[j],rew_err[j],fwhm[j],fwhm_err[j],frat_nvciv[j],snr_1700[j],w1mw2[j],w2mw3[j],j) )

   ax1 = plt.subplot(gs[wind,:])
   plt.ylim(-0.05*ymax, ymax)
   plt.xlim(1140,1980)
   if rew[j] > 0.0:
      mask = wz < 1600.0
      plt.plot(wz[mask],lfits[mask],c='m',lw=0.7,ls='-',alpha=0.7,zorder=1)
   plt.plot(wz,cont_tot,c='orange',lw=0.7,ls='-',alpha=1.0,zorder=2)
   plt.plot(wz,sm_spec,c='k',lw=1.0,ls='steps-mid')
   plt.plot(wz,sm_spece,c='limegreen',lw=0.4,ls=':')
   plt.axvline(lya,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(lyb,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(he2,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(o6[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(o6[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si3[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(n5[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(n5[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(mg2[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(mg2[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(n3[1],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(c3[0],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(c3[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(c3[2],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(n2[0],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(n2[1],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(n2[2],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si2b[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si2b[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si2b[2],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si2b[3],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si4[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(si4[1],c='c',lw=0.6,ls='--', alpha=0.4)
#plt.axvline(o4[2],c='c',lw=0.6,ls='--', alpha=0.4)
#plt.axvline(o4[3],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(s4[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(s4[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(c4[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(c4[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(c2,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(fe2[1],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(fe3[2],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(fe3[3],c='c',lw=0.6,ls='--', alpha=0.4)
#   plt.axvline(fe3[4],c='indigo',lw=0.6,ls='--',alpha=0.4)
#   plt.axvline(fe3[5],c='indigo',lw=0.6,ls='--',alpha=0.4)
#   plt.axvline(fe3[6],c='indigo',lw=0.6,ls='--',alpha=0.4)
   plt.axvline(o1,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(n4,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(o3,c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(al3[0],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axvline(al3[1],c='c',lw=0.6,ls='--', alpha=0.4)
   plt.axhline(0.0,c='k',lw=0.6,ls='-')
#plt.text(1034,  0.93*ymax,'OVI',ha='center',va='center',fontsize=9.0)
#plt.text(1069,  0.88*ymax,'SIV',ha='center',va='center',fontsize=9.0)
#plt.text(1128,  0.93*ymax,'FeIII',ha='center',va='center',fontsize=9.0)
   plt.text(1197,  0.93*ymax,'Ly$\\alpha$',ha='center',va='center',fontsize=9.0)
   plt.text(1240,  0.93*ymax,'NV',ha='center',va='center',fontsize=9.0)
   plt.text(1264,  0.89*ymax,'SiII',ha='center',va='center',fontsize=9.0)
   plt.text(1303,  0.93*ymax,'OI',ha='center',va='center',fontsize=9.0)
   plt.text(1335,  0.93*ymax,'CII',ha='center',va='center',fontsize=9.0)
   plt.text(1400,  0.93*ymax,'SiIV',ha='center',va='center',fontsize=9.0)
   plt.text(1486,  0.93*ymax,'NIV]',ha='center',va='center',fontsize=9.0)
   plt.text(1527,  0.87*ymax,'SiII',ha='center',va='center',fontsize=9.0)
   plt.text(1549,  0.93*ymax,'CIV',ha='center',va='center',fontsize=9.0)
   plt.text(1743,  0.93*ymax,'NIII]',ha='center',va='center',fontsize=9.0)
   plt.text(1790,  0.93*ymax,'FeII',ha='center',va='center',fontsize=9.0)
   plt.text(1625,  0.93*ymax,'HeII',ha='center',va='center',fontsize=9.0)
   plt.text(1672,  0.93*ymax,'OIII]',ha='center',va='center',fontsize=9.0)
   plt.text(1859,  0.93*ymax,'AlIII',ha='center',va='center',fontsize=9.0)
   plt.text(1888,  0.87*ymax,'SiIII]',ha='center',va='center',fontsize=9.0)
   plt.text(1909,  0.93*ymax,'CIII]',ha='center',va='center',fontsize=9.0)
#   plt.text(2798,  0.93*ymax,'MgII',ha='center',va='center',fontsize=9.0)
#   plt.text(1926,  0.88*ymax,'FeIII',ha='center',va='center',fontsize=9.0)
   if (wind <= 3) and j != nqsos-1: plt.setp(ax1.get_xticklabels(), visible=False)
   if wind == 4 or j == nqsos-1: plt.xlabel('Rest Wavelength (A)')
   plt.ylabel('Flux')
   xloc = 1790
   if imw3[j] < 0.0:
      plt.text(xloc, 0.68*ymax, 'J%s   z$_e$=%.2f\n  i=%.1f  i-w3=---  w3-w4=---  r-z=%.1f\n REW=%.0f+/-%.0f  FWHM=%.0f+/-%.0f\n nv/civ=%.1f  kt=%.2f  ph=%.1f  bal_vi=%.1f' % (name[j],zdr12[j],ABi[j],rmz[j],rew[j],rew_err[j],fwhm[j],fwhm_err[j],frat_nvciv[j],kurt80[j],peak[j],bal_vi[j]), ha='center', va='center', fontsize=8.7)
   elif w3mw4[j] < 0.0:
      plt.text(xloc, 0.68*ymax, 'J%s   z$_e$=%.2f\n  i=%.1f  i-w3=%.1f  w3-w4=---  r-z=%.1f\n REW=%.0f+/-%.0f  FWHM=%.0f+/-%.0f\n nv/civ=%.1f  kt=%.2f  ph=%.1f  bal_vi=%.1f' % (name[j],zdr12[j],ABi[j],imw3[j],rmz[j],rew[j],rew_err[j],fwhm[j],fwhm_err[j],frat_nvciv[j],kurt80[j],peak[j],bal_vi[j]), ha='center', va='center', fontsize=8.7)
   else:
      plt.text(xloc, 0.68*ymax, 'J%s   z$_e$=%.2f\n  i=%.1f  i-w3=%.1f  w3-w4=%.1f  r-z=%.1f\n REW=%.0f+/-%.0f  FWHM=%.0f+/-%.0f\n nv/civ=%.1f  kt=%.2f  ph=%.1f  bal_vi=%.1f' % (name[j],zdr12[j],ABi[j],imw3[j],w3mw4[j],rmz[j],rew[j],rew_err[j],fwhm[j],fwhm_err[j],frat_nvciv[j],kurt80[j],peak[j],bal_vi[j]), ha='center', va='center', fontsize=8.7)

   if wind == 4:
      plt.savefig(pdf, format='pdf', bbox_inches='tight')
      plt.close()
      wind = -1

if wind != -1:
   if (wind < 4): plt.setp(ax1.get_xticklabels(), visible=True)
   plt.xlabel('Rest Wavelength (A)')
   plt.savefig(pdf, format='pdf', bbox_inches='tight')
   plt.close()
pdf.close()

print('# All, # BALs_VI, # T2, <imw3>, #color =', count, nbalsall, nT2all, median(colorall), len(colorall) )
print('# ERQs, # BALs_VI, # T2, <imw3>, #color =', nerqs, nbals, nT2, median(colorerqs), len(colorerqs))
print('# core, # BALs_VI core, # T2 core, <imw3>, #color, median(redshift) =', ncore, nbcore, nT2core, median(colorcore), len(colorcore), median(redshift))
"""
newtdata = tdata[ind]
hdu = pyfits.BinTableHDU(newtdata)
hdu.writeto("Type2_all_v11.fits")
"""
