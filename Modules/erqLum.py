# Luminosity of Quasars using Flat Lambda Cold Dark Matter model 
def Lum(m_AB, z_dr12, lambda_obs, lambda_nominal):
    import numpy as np
    import astropy.cosmology as cosmo
    row = len(z_dr12)
    F_nu = 10**(-(48.6+m_AB)*0.4)
    c = 2.998e10 # cm/sec
    H0 = 67.3 # km/s MP^-1
    omega_m0 = 0.315
    omega_gamma0=0.685
    D_p=np.empty([row])
    D_Lum=np.empty([row])
    model = cosmo.FlatLambdaCDM(H0=H0, Om0=omega_m0)
    for i in range(row):
        D_Lum[i]= model.luminosity_distance(z_dr12[i]).value
        D_Lum[i]*=3085677600000000000000000 # Mpc to cm
    lambda_rest = lambda_obs/(1+z_dr12)
    F_obs = (c/lambda_obs**2)*F_nu
    Lum_rest = F_obs*4*np.pi*D_Lum**2*(1+z_dr12)
    Lum_nominal = Lum_rest*(lambda_nominal/lambda_rest)**(-0.65)
    return np.log10(8*Lum_rest*lambda_rest)
