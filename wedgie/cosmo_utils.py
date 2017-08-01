import numpy as np
import gen_utils as gu
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from astropy import constants as c
import aipy

# options for cosmologies: ["FLRW", "LambdaCDM", "FlatLambdaCDM", "wCDM", "FlatwCDM",
#           "Flatw0waCDM", "w0waCDM", "wpwaCDM", "w0wzCDM", "WMAP5", "WMAP7",
#           "WMAP9", "Planck13", "Planck15", "default_cosmology"]

# Temperature conversion constants

PAPER_BEAM_POLY = np.array([ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03])

HERA_BEAM_POLY = np.array([  8.07774113e+08,  -1.02194430e+09,   
    5.59397878e+08,  -1.72970713e+08, 3.30317669e+07,  -3.98798031e+06,   
    2.97189690e+05,  -1.24980700e+04, 2.27220000e+02]) # See HERA Memo #27

HERA_OP_OPP = 2.175 # See HERA Memo #27
PAPER_OP_OPP = 2.35

NOISE_EQUIV_BW = {
    'blackman': 1.73,
    'blackman-harris': 2.01,
    'gaussian0.4': 1.45,
    'hamming': 1.37,
    'hanning': 1.50,
    'kaiser2': 1.50,
    'kaiser3': 1.80,
    'parzen': 1.33,
    'none': 1,
}

# Cosmology routines

def f2z(fq):
    """
    Convert frequency to redshift for 21cm line.
    """
    F21 = 1.42040575177*u.GHz
    return (F21 / fq - 1.)

def f2eta(f):
    """
    Convert an array of frequencies to an array of etas (freq^-1) 
    corresponding to the bins that an FFT generates.
    
    Expects GIGAHERTZ frequencies
    """
    return np.fft.fftfreq(f.shape[-1], f[1]-f[0])*u.ns

def dL_df(z):
    """
    [h^-1 Mpc]/GHz, from Furlanetto et al. (2006)
    """
    return (1.7 / 0.1) * (np.sqrt((1+z)/10.) * (1./np.sqrt(cosmo.Om0/0.15)) * 1e3 ) * u.Mpc/(u.GHz)

def dL_dth(z):
    """
    [h^-1 Mpc]/radian, from Furlanetto et al. (2006)
    """
    #0.00291 == aipy.const.arcmin
    return 1.9 * (1./0.000291) * np.power((1+z)/10.,0.2) *u.Mpc/(u.radian)

def X2Y(z):
    """
    [h^-3 Mpc^3] / [str * GHz] scalar conversion between observing and cosmological coordinates
    """
    return np.power(dL_dth(z),2.) * dL_df(z)

def dk_du(z):
    """
    2pi * [h Mpc^-1] / [wavelengths], valid for u >> 1.
    """
    # from du = 1/dth, which derives from du = d(sin(th)) using the small-angle approx
    return 2*np.pi / dL_dth(z)

def dk_deta(z):
    """
    2pi * [h Mpc^-1] / [GHz^-1]
    """
    return (2*np.pi / dL_df(z)).to(u.GHz/u.Mpc)

def eta2kpl(etas,z):
    """
    Convert an array of etas (freq^-1) to k_parallel (h/Mpc)
    """
    return dk_deta(z) * etas


def jy2T(f, bm_poly=HERA_BEAM_POLY):
    """
    Return [mK] / [Jy] for a beam size vs. frequency (in GHz) defined by the
    polynomial bm_poly.  Default is a poly fit to the HERA primary beam.
    """
    lam = c.c/ (f*u.GHz)
    bm = np.polyval(bm_poly, f)
    return 1e-23 * np.power(lam,2.) / (2 * c.k_B * bm) * 1e3

def freq2kpl(freqs,fold=False):
    """
    Convert an array of frequencies to k_parallel in h/Mpc.
    Assumes redshift based on central frequency, and equally-spaced frequencies.
    
    fold=True returns the positive-half of the k_parallel array.
    """
    cen_fq = gu.findMiddle(freqs)
    z = f2z(cen_fq)
    etas = f2eta(freqs)
    kpl = np.fft.fftshift(eta2kpl(etas,z))
    if fold:
        return kpl[kpl.shape[0]/2:]*(1./u.Mpc)
    else:
        return kpl*(1./u.Mpc)

def uv2kpr(blmag,cen_fq):
    """
    Compute k_perpendicular from the magnitude of the baseline vector and the
    central frequency of observation.
    
    blmag is an np.array of baseline vector magnitudes, units of length
    cen_fq is the central frequency.
    
    Returns k_perpendicular in units of h/Mpc
    """
    z = f2z(cen_fq)
    lam = c.c/(cen_fq)
    uvmag = blmag/lam
    kpr = 2*np.pi*uvmag/(cosmo.comoving_transverse_distance(z)*cosmo.h)
    return kpr.to(1./u.Mpc)

def scale2Pk_crude(sq_arr,cen_fq,dfq,avg_scale=0.03):
    """
    Crudely scales squared HERA visibilities to an approximate Kelvin^2
    level using the Beardesley et al. and Parsons et al. HERA Memos:
    
    reionization.org/wp-content/uploads/2017/04/HERA19_Tsys_3April2017.pdf,  reionization.org/wp-content/uploads/2013/03/Power_Spectrum_Normalizations_for_HERA.pdf
    
    Input: sq_arr: gridded power spectra in raw visibility^2 units
           cen_fq: central frequency IN GIGAHERTZ
           dfq: channel width IN GIGAHERTZ
           avg_scale: divide this number from every antenna to go do Jy
    """
    jy2_arr = sq_arr/np.power(avg_scale,2.)
    fq0 = .150 # GHz
    z = f2z(fq0)
    B = dfq * NOISE_EQUIV_BW['none'] #proper normaliztion
    bm_HERA = np.polyval(HERA_BEAM_POLY, cen_fq) * 2.35
    scalar_HERA = X2Y(z) * bm_HERA * B #[h^-3 Mpc^3] / [sr * GHz]
    #XXX I need to double-check the below factors
    jy_to_mK2_HERA = scalar_HERA * jy2T(cen_fq)**2
    return jy2_arr*np.power(jy_to_mK2_HERA,2.)
    