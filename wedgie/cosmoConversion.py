#!/usr/bin/env python

import numpy as np
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from astropy import constants as c

F21 = 1.42040575177*u.GHz

def findMiddle(input_list):
    """
    This method returns the central value of a list.
    If there is an even number of values, it returns
    the average of the central two. 
    
    i.e. this should be used for constantly-spaced-value lists.
    
    It's kinda stupid.
    """
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return input_list[int(middle - .5)]
    else:
        return np.mean([input_list[int(middle)], input_list[int(middle-1)]])

class Visibility(object):
    """
    Interferometric observable V(nu,vec(b)).
    
    Cosmological conversions rely on universal isotropy, so we only care about
    the magnitude of the baseline vector, not its direction.
    """
    def __init__(self,freq=None,blmag=None):
        self._freq = freq
        self._blmag = blmag
        
    @property
    def freq(self):
        return self._freq
    
    @freq.setter
    def freq(self,newredshift):
        """
        Converts redshift of 21cm emission to frequency
        """
        self._freq = F21/(1+newredshift)
     
    @property
    def blmag(self):
        return self._blmag
    
    @blmag.setter
    def blmag(self,new_blmag):
        self._blmag = new_blmag 
     
    @property
    def redshift(self):
        return (F21/self.freq) - 1.
        
    @property
    def eta(self):
        return np.fft.fftfreq(self.freq.shape[-1], self.freq[1]-self.freq[0]) * (1/self.freq.unit)
    
    def dL_df(self):
        """
        [h^-1 Mpc]/GHz, from Furlanetto et al. (2006)
        """
        return (1.7 / 0.1) * (np.sqrt((1+self.redshift)/10.) * (1./np.sqrt(cosmo.Om0/0.15)) * 1e3 ) * u.Mpc/(u.GHz)
    
    def dL_dth(self):
        """
        [h^-1 Mpc]/radian, from Furlanetto et al. (2006)
        """
        #0.00291 == aipy.const.arcmin
        return 1.9 * (1./0.000291) * np.power((1+self.redshift)/10.,0.2)*u.Mpc/(u.radian)

    def dk_du(self):
        """
        2pi * [h Mpc^-1] / [wavelengths], valid for u >> 1.
        """
        # from du = 1/dth, which derives from du = d(sin(th)) using the small-angle approx
        return (2*np.pi / self.dL_dth()).decompose()

    def dk_deta(self):
        """
        2pi * [h Mpc^-1] / [GHz^-1]
        """
        return (2*np.pi / self.dL_df()).to(u.GHz/u.Mpc)
    
    def eta2kpl(self):
        """
        Convert an array of etas (freq^-1) to k_parallel (h/Mpc)
        """
        return (self.dk_deta()*self.eta).decompose()
    
    def freq2kpl(self,fold=False):
        """
        Convert an array of frequencies to k_parallel in h/Mpc.
        Assumes redshift based on central frequency, and equally-spaced frequencies.
    
        fold=True returns the positive-half of the k_parallel array.
        """
        cen_fq = findMiddle(self.freq)
        kpl = np.fft.fftshift(self.eta2kpl())
        if fold:
            return kpl[kpl.shape[0]/2:]*(1./u.Mpc)
        else:
            return kpl*(1./u.Mpc)
    
    def uv2kpr(self):
        """
        Compute k_perpendicular from the magnitude of the baseline vector and the
        central frequency of observation.
    
        blmag is an np.array of baseline vector magnitudes, units of length
        cen_fq is the central frequency.
    
        Returns k_perpendicular in units of h/Mpc
        """
        lam = c.c/self.freq
        uvmag = self.blmag/lam
        xi = np.zeros_like(self.redshift)
        for i,z in enumerate(self.redshift):
            xi[i] = cosmo.comoving_transverse_distance(z)
        kpr = 2*np.pi*uvmag/(xi*cosmo.h)
        return kpr.to(1./u.Mpc)
    