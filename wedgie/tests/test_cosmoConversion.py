#!/usr/bin/env python

"""
Tests for cosmoConversion.py
"""

import nose.tools as nt
from wedgie import cosmoConversion as cc
from astropy import units as u

class TestMethods(object):
    
    def setUp(self):
        self.testarray_odd = [0,1,2]
        self.testarray_even= [0,1,2,3]
        
    def test_findMiddle(self):
        nt.assert_equal(cc.findMiddle(self.testarray_odd), 1)
        nt.assert_equal(cc.findMiddle(self.testarray_even), 1.5)
    
class Test_Visibility(object):
    """
    Tests for the Visibility object.
    Assumes astropy's Planck15 cosmology.
    """
    def setUp(self):
        self.frequencies = [100.,120.,160.,180.,200.]*u.MHz
        self.blmag= 30.*u.m
        self.mainvis = cc.Visibility(freq=self.frequencies, blmag=self.blmag)
        
    def test_freq(self):
        vis = cc.Visibility(freq=self.frequencies)
        nt.assert_true(vis.freq[1] == self.frequencies[1])
    
    def test_blmag(self):
        vis = cc.Visibility(blmag=self.blmag)
        nt.assert_equal(vis.blmag, self.blmag)
        vis.blmag = 1*u.m
        nt.assert_equal(vis.blmag, 1*u.m)
    
    def test_redshift(self):
        z = self.mainvis.redshift
        nt.assert_almost_equals(z[0].value,13.20405752,places=3)
    
    def test_eta(self):
        etas = self.mainvis.eta
        nt.assert_equal(etas[0],0.)
    
    #XXX NEED NUMERICAL VALUES TO TEST AGAINST FOR THE REST OF THESE
    def test_dL_df(self):
        dLdf = self.mainvis.dL_df()
    
    def test_dL_dth(self):
        dLdth = self.mainvis.dL_dth()
    
    def test_dk_du(self):
        dkdu = self.mainvis.dk_du()
    
    def test_dk_deta(self):
        dkdeta = self.mainvis.dk_deta()
    
    def test_eta2kpl(self):
        kpl = self.mainvis.eta2kpl()
    
    def test_freq2kpl(self):
        kpl_all = self.mainvis.freq2kpl()
        kpl_fld = self.mainvis.freq2kpl(fold=True)
    
    # THIS FAILS
    def test_uv2kpr(self):
        kpr = self.mainvis.uv2kpr()
    
        
    
    