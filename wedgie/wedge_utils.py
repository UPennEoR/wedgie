"""
This module contains functions to create and plot wedges/pitchforks from HERA data files.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
Paul La Plante <plaplant_at_sas.upenn.edu>
"""

import aipy
from pyuvdata import UVData
import pyuvdata.utils as uvutils

import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import scipy.constants as sc

import gen_utils as gu
import cosmo_utils as cu

import decimal

import hera_cal

# For Interactive Development
from IPython import embed


class Wedge(object):
    def __init__(self, args, files, calfile, pol, ex_ants, freq_range, history):
        self.args = args
        self.files = files
        self.pol = pol
        self.calfile = calfile
        self.history = history
        self.freq_range = freq_range
        self.ex_ants = ex_ants

        self.npz_name = str()
        self.caldata = tuple()
        self.info = None
        self.data = None
        self.flags = None
        self.aa = None
        self.bl_length = int()
        self.vis_sq_bl = None
        self.vis_sq_slope = None
        self.vis_sq_antpair = None
        self.fft_2Ddata = None
        self.zenith = None

        self.wedgeslices = list()
        self.cwedgeslices = list()
        self.antpairslices = list()

        self.lst = int()
        self.lst_range_str = list()
        self.lst_range_num = list()
        self.times = list()

        self.delays = list()

        decimal.getcontext().prec = 6
        if self.calfile is not None:
            self.calculate_caldata()

    # Methods common throughout Wedge Creation
    def name_npz(self, tag):
        file_start = self.files[self.files.keys()[0]][0].split('/')[-1].split('.')
        file_end = self.files[self.files.keys()[0]][-1].split('/')[-1].split('.')

        zen_day = file_start[:2]
        time_range = ['{start}_{end}'.format(start=file_start[2], end=file_end[2])]
        HH = [file_start[4]]
        ext = [file_start[5]]
        pol = [self.pol]
        freq_range = ['{start}_{end}'.format(start=self.freq_range[0], end=self.freq_range[1])]
        tag = [tag]

        npz_name = zen_day + time_range + pol + freq_range + HH + ext + tag + ["npz"]

        if self.ex_ants:
            ex_ants = [str(ant) for ant in self.ex_ants]
            ex_ants = "_".join(ex_ants)
            npz_name.insert(4, ex_ants)
        else:
            ex_ants = ["None"]
            npz_name.insert(4, ex_ants)

        npz_name = '.'.join(npz_name)

        self.npz_name = self.args.path + npz_name

        print(self.npz_name)

    def calculate_caldata(self, uvd=None):
        """
        Returns a dictionary of baseline lengths and the corresponding pairs. The data is based
        on a calfile. ex_ants is a list of integers that specify antennae to be exlcuded from
        calculation.

        Requires cal file to be in PYTHONPATH.
        """
        if self.calfile is not None:
            try:
                print('Reading calfile: {cfile}'.format(cfile=self.calfile))
                exec("import {cfile} as cal".format(cfile=self.calfile))
                antennae = cal.prms['antpos_ideal']
            except ImportError:
                raise Exception("Unable to import {cfile}.".format(cfile=self.calfile))
        else:
            # build antenna positions from data file itself
            print('Generating calibration information from data file')
            antennae = {}
            lat, lon, alt = uvd.telescope_location_lat_lon_alt
            for i, antnum in enumerate(uvd.antenna_numbers):
                pos = uvd.antenna_positions[i, :] + uvd.telescope_location
                xyz = uvutils.ENU_from_ECEF(pos, latitude=lat, longitude=lon, altitude=alt)
                antennae[antnum] = {'top_x': xyz[0], 'top_y': xyz[1], 'top_z': xyz[2]}

        # Remove all placeholder antennae from consideration
        # Remove all antennae from ex_ants from consideration
        ants = []
        for ant in antennae.keys():
            if (not antennae[ant]['top_z'] < 0) and (ant not in self.ex_ants):
                ants.append(ant)

        # Form pairs of antennae
        # Store unique baselines and slopes for later use
        pairs, baselines, slopes = {}, [], []
        for ant_i in ants:
            for ant_j in ants:
                if (ant_i >= ant_j):
                    continue
                pair = (ant_i, ant_j)

                # Calculate baseline length
                dx = decimal.Decimal(antennae[pair[0]]['top_x'] - antennae[pair[1]]['top_x'])
                dy = decimal.Decimal(antennae[pair[0]]['top_y'] - antennae[pair[1]]['top_y'])
                baseline = float((dx**2 + dy**2).sqrt())
                baselines.append(baseline)

                # Calculate slope between an antenna pair
                dy = decimal.Decimal(antennae[pair[1]]['top_y'] - antennae[pair[0]]['top_y'])
                dx = decimal.Decimal(antennae[pair[1]]['top_x'] - antennae[pair[0]]['top_x'])
                if dx != 0:
                    slope = float(dy / dx)
                else:
                    slope = np.inf
                slopes.append(slope)

                pairs[pair] = (baseline, slope)

        # Remove duplicates baseline and slope values
        baselines = set(baselines)
        slopes = set(slopes)

        # Initalize antdict with baselines as keys and empty lists as values
        antdict = {baseline: [] for baseline in baselines}

        # Add pairs to the list of their respective baseline
        for pair in pairs:
            baseline = pairs[pair][0]
            antdict[baseline].append(pair)

        # Initialize slopedict with baselines for keys and the dictionary of slopes for each value
        slopedict = {}
        for baseline in baselines:
            slopedict[baseline] = {slope: [] for slope in slopes}

        # Add pairs to their respective slope within their respective baseline
        for pair in pairs:
            baseline = pairs[pair][0]
            slope = pairs[pair][1]
            slopedict[baseline][slope].append(pair)

        for baseline in slopedict.copy():
            for slope in slopedict[baseline].copy():
                if slopedict[baseline][slope] == []:
                    del slopedict[baseline][slope]

        self.caldata = (antdict, slopedict, pairs, sorted(list(baselines)), sorted(list(slopes)))

    def cleanfft(self, pair, clean=1e-3):
        # FFT
        w = aipy.dsp.gen_window(self.data[pair][self.pol].shape[-1], window='blackman-harris')
        _dw = np.fft.ifft(self.data[pair][self.pol] * w)

        # CLEAN
        _ker = np.fft.ifft(self.flags[pair][self.pol] * w)
        gain = aipy.img.beam_gain(_ker)
        for time in range(_dw.shape[0]):
            _dw[time, :], info = aipy.deconv.clean(_dw[time, :], _ker[time, :], tol=clean)
            _dw[time, :] += info['res'] / gain

        self.fft_2Ddata = np.ma.array(_dw)

    def phasing(self, ntimes, antpair):
        for i in range(ntimes):
            if i != 0:
                old_zenith = self.zenith
            time = self.info['times'][i]
            self.aa.set_jultime(time)

            self.lst = self.aa.sidereal_time()
            self.lst_range_str.append(str(self.lst))
            self.lst_range_num.append(self.lst)

            self.zenith = aipy.phs.RadioFixedBody(self.lst, self.aa.lat)
            self.zenith.compute(self.aa)

            if i % 2:
                v1 = self.fft_2Ddata[i - 1, :]
                phase_correction = np.conj(self.aa.gen_phs(self.zenith, antpair[0], antpair[1])) * self.aa.gen_phs(old_zenith, antpair[0], antpair[1])
                v2 = self.fft_2Ddata[i, :] * phase_correction[self.freq_range[0]:self.freq_range[1]]
                self.vis_sq_antpair[i // 2, :] = np.conj(v1) * v2

    # Forming Stokes Parameters:
    def form_stokesI(self):
        """Calculate I (VI = Vxx + Vyy)"""
        uvxx = UVData()
        uvxx.read_miriad(self.files['xx'])
        uvyy = UVData()
        uvyy.read_miriad(self.files['yy'])

        # compute redundancy information if it hasn't been done yet
        if self.caldata == tuple():
            self.calculate_caldata(uvxx)

        # get metadata
        info = {}
        # convert from Hz -> GHz
        info['freqs'] = uvxx.freq_array[0, :] / 1e9
        info['times'] = np.unique(uvxx.time_array)

        # get data and flags
        dI = {}
        fI = {}
        for key, dxx in uvxx.antpairpol_iter():
            # ignore auto correlations
            if key[0] == key[1]:
                continue
            ind1, ind2, ipol = uvxx._key2inds(key)
            for ind in [ind1, ind2]:
                if len(ind) == 0:
                    continue
                # assumes that xx and yy files only have a single polarization in them
                dyy = uvyy.data_array[ind, 0, :, 0]
                fxx = uvxx.flag_array[ind, 0, :, 0]
                fyy = uvyy.flag_array[ind, 0, :, 0]

                antkey = key[:2]
                dI[antkey] = {'I': dxx + dyy}
                fI[antkey] = {'I': fxx + fyy}

        # clean up
        del uvxx, uvyy

        # assign to object variables
        self.info = info
        self.data = dI
        self.flags = fI


    def form_stokesQ(self):
        """Calculate Q (VQ = Vxx - Vyy)"""
        uvxx = UVData()
        uvxx.read_miriad(self.files['xx'])
        uvyy = UVData()
        uvyy.read_miriad(self.files['yy'])

        # compute redundancy information if it hasn't been done yet
        if self.caldata == tuple():
            self.calculate_caldata(uvxx)

        # get metadata
        info = {}
        # convert from Hz -> GHz
        info['freqs'] = uvxx.freq_array[0, :] / 1e9
        info['times'] = np.unique(uvxx.time_array)

        # get data and flags
        dQ = {}
        fQ = {}
        for key, dxx in uvxx.antpairpol_iter():
            # ignore auto correlations
            if key[0] == key[1]:
                continue
            ind1, ind2, ipol = uvxx._key2inds(key)
            for ind in [ind1, ind2]:
                if len(ind) == 0:
                    continue
                # assumes that xx and yy files only have a single polarization in them
                dyy = uvyy.data_array[ind, 0, :, 0]
                fxx = uvxx.flag_array[ind, 0, :, 0]
                fyy = uvyy.flag_array[ind, 0, :, 0]

                antkey = key[:2]
                dQ[antkey] = {'Q': dxx - dyy}
                fQ[antkey] = {'Q': fxx + fyy}

        # clean up
        del uvxx, uvyy

        # assign to object variables
        self.info = info
        self.data = dQ
        self.flags = fQ

    def form_stokesU(self):
        """Calculate U (VU = Vxy + Vyx)"""
        uvxy = UVData()
        uvxy.read_miriad(self.files['xy'])
        uvyx = UVData()
        uvyx.read_miriad(self.files['yx'])

        # compute redundancy information if it hasn't been done yet
        if self.caldata == tuple():
            self.calculate_caldata(uvxy)

        # get metadata
        info = {}
        # convert from Hz -> GHz
        info['freqs'] = uvxy.freq_array[0, :] / 1e9
        info['times'] = np.unique(uvxy.time_array)

        # get data and flags
        dU = {}
        fU = {}
        for key, dxy in uvxy.antpairpol_iter():
            # ignore auto correlations
            if key[0] == key[1]:
                continue
            ind1, ind2, ipol = uvxy._key2inds(key)
            for ind in [ind1, ind2]:
                if len(ind) == 0:
                    continue
                # assumes that xy and yx files only have a single polarization in them
                dyx = uvyx.data_array[ind, 0, :, 0]
                fxy = uvxy.flag_array[ind, 0, :, 0]
                fyx = uvyx.flag_array[ind, 0, :, 0]

                antkey = key[:2]
                dU[antkey] = {'U': dxy + dyx}
                fU[antkey] = {'U': fxy + fyx}

        # clean up
        del uvxy, uvyx

        # assign to object variables
        self.info = info
        self.data = dU
        self.flags = fU

    def form_stokesV(self):
        """Calculate V (VV = -i*Vxy + i*Vyx)"""
        uvxy = UVData()
        uvxy.read_miriad(self.files['xy'])
        uvyx = UVData()
        uvyx.read_miriad(self.files['yx'])

        # compute redundancy information if it hasn't been done yet
        if self.caldata == tuple():
            self.calculate_caldata(uvxy)

        # get metadata
        info = {}
        # convert from Hz -> GHz
        info['freqs'] = uvxy.freq_array[0, :] / 1e9
        info['times'] = np.unique(uvxy.time_array)

        # get data and flags
        dV = {}
        fV = {}
        for key, dxy in uvxy.antpairpol_iter():
            # ignore auto correlations
            if key[0] == key[1]:
                continue
            ind1, ind2, ipol = uvxy._key2inds(key)
            for ind in [ind1, ind2]:
                if len(ind) == 0:
                    continue
                # assumes that xy and yx files only have a single polarization in them
                dyx = uvyx.data_array[ind, 0, :, 0]
                fxy = uvxy.flag_array[ind, 0, :, 0]
                fyx = uvyx.flag_array[ind, 0, :, 0]

                antkey = key[:2]
                dV[antkey] = {'V': -1j * dxy + 1j * dyx}
                fV[antkey] = {'V': fxy + fyx}

        # clean up
        del uvxy, uvyx

        # assign to object variables
        self.info = info
        self.data = dV
        self.flags = fV

    # Load data of one polarization:
    def load_file(self):
        """Loads data with given polarization, self.pol, from files, self.files"""
        uv = UVData()
        uv.read_miriad(self.files[self.pol])

        # compute redundancy information if it hasn't been done yet
        if self.caldata == tuple():
            self.calculate_caldata(uv)

        # get metadata
        info = {}
        # convert from Hz -> GHz
        info['freqs'] = uv.freq_array[0, :] / 1e9
        info['times'] = np.unique(uv.time_array)

        # get data and flags
        d = {}
        f = {}
        for key, data in uv.antpairpol_iter():
            # ignore auto correlations
            if key[0] == key[1]:
                continue
            ind1, ind2, ipol = uv._key2inds(key)
            for ind in [ind1, ind2]:
                if len(ind) == 0:
                    continue
                flags = uv.flag_array[ind, 0, :, ipol]
                antkey = key[:2]
                polkey = key[2].lower()
                d[antkey] = {polkey: data}
                f[antkey] = {polkey: flags}

        # clean up
        del uv

        # assign to object variables
        self.info = info
        self.data = d
        self.flags = f

    # Format flags for correct application:
    def apply_flags(self):
        """Turns flags from False/True --> 1/0."""
        for i in self.flags:
            for j in self.flags[i]:
                    self.flags[i][j] = np.logical_not(self.flags[i][j]).astype(int)

        for pair in self.caldata[2]:
            self.data[pair][self.pol] *= self.flags[pair][self.pol]

    # Types of Wedge creating methods:
    def timeavg(self):
        """Take the log10 of the fftshift of the absolute value of the time averaged visibility data."""
        self.wedgeslices = np.log10(np.fft.fftshift(np.abs(np.nanmean(self.cwedgeslices, axis=1)), axes=1))

    def blavg(self):
        """Take the log10 of the fftshift of the absolute value of the visibility data."""
        self.wedgeslices = np.log10(np.fft.fftshift(np.abs(self.cwedgeslices), axes=1))

    def flavors(self):
        """Create a pitchfork averaged over baseline orientation and in time"""
        self.info['freqs'] = self.info['freqs'][self.freq_range[0]:self.freq_range[1]]
        for antpair in self.data.keys():
            self.data[antpair][self.pol] = self.data[antpair][self.pol][:, self.freq_range[0]:self.freq_range[1]]
            self.flags[antpair][self.pol] = self.flags[antpair][self.pol][:, self.freq_range[0]:self.freq_range[1]]

        ntimes = len(self.info['times'])
        nchan = len(self.info['freqs'])

        if self.calfile is not None:
            uv = aipy.miriad.UV(self.files[self.files.keys()[0]][0])
            self.aa = aipy.cal.get_aa(self.calfile, uv['sdf'], uv['sfreq'], uv['nchan'])
        else:
            uv = UVData()
            uv.read_miriad(self.files[self.files.keys()[0]][0])
            # convert from Hz -> GHz
            freqs = uv.freq_array[0, :] / 1e9
            self.aa = hera_cal.utils.get_aa_from_uv(uv, freqs)
        del uv
        self.aa.set_active_pol(self.pol)

        for baselength in sorted(self.caldata[1]):
            self.vis_sq_bl = np.zeros((ntimes // 2, nchan), dtype=np.complex128)
            for slope in sorted(self.caldata[1][baselength]):
                self.vis_sq_slope = np.zeros((ntimes // 2, nchan), dtype=np.complex128)
                for antpair in self.caldata[1][baselength][slope]:
                    self.vis_sq_antpair = np.zeros((ntimes//2, nchan), dtype=np.complex128)
                    self.cleanfft(antpair)
                    self.phasing(ntimes, antpair)
                    self.vis_sq_slope += self.vis_sq_antpair
                self.vis_sq_slope /= len(self.caldata[1][baselength])
                self.wedgeslices.append(np.log10(np.fft.fftshift(np.abs(np.nanmean(self.vis_sq_slope, axis=0)))))
                print('Wedgeslice for baseline {} and slope {} complete.'.format(baselength, slope))

    def bltype(self):
        """Create a plot with a wedgeslice for each antenna pair from a given baselength"""
        self.info['freqs'] = self.info['freqs'][self.freq_range[0]:self.freq_range[1]]
        for antpair in self.data.keys():
            self.data[antpair][self.pol] = self.data[antpair][self.pol][:, self.freq_range[0]:self.freq_range[1]]
            self.flags[antpair][self.pol] = self.flags[antpair][self.pol][:, self.freq_range[0]:self.freq_range[1]]

        ntimes = len(self.info['times'])
        nchan = len(self.info['freqs'])

        if self.calfile is not None:
            uv = aipy.miriad.UV(self.files[self.files.keys()[0]][0])
            self.aa = aipy.cal.get_aa(self.calfile, uv['sdf'], uv['sfreq'], uv['nchan'])
        else:
            uv = UVData()
            uv.read_miriad(self.files[self.files.keys()[0]][0])
            # convert from Hz -> GHz
            freqs = uv.freq_array[0, :] / 1e9
            self.aa = hera_cal.utils.get_aa_from_uv(uv, freqs)
        del uv
        self.aa.set_active_pol(self.pol)

        bl_type = self.args.bl_type - 1
        self.bl_length = self.caldata[3][bl_type]
        antpairs = self.caldata[0][self.bl_length]

        for antpair in antpairs:
            self.vis_sq_antpair = np.zeros((ntimes // 2, nchan), dtype=np.complex128)
            self.cleanfft(antpair)
            self.phasing(ntimes, antpair)

            self.antpairslices.append(np.log10(np.fft.fftshift(np.abs(self.vis_sq_antpair), axes=1)))
            print('Wedgeslice for antenna pair {} complete.'.format(antpair))

    # Wedge Creation
    def starter(self):
        self.info['freqs'] = self.info['freqs'][self.freq_range[0]:self.freq_range[1]]
        for antpair in self.data.keys():
            self.data[antpair][self.pol] = self.data[antpair][self.pol][:, self.freq_range[0]:self.freq_range[1]]
            self.flags[antpair][self.pol] = self.flags[antpair][self.pol][:, self.freq_range[0]:self.freq_range[1]]

        ntimes = len(self.info['times'])
        nchan = len(self.info['freqs'])

        if self.calfile is not None:
            uv = aipy.miriad.UV(self.files[self.files.keys()[0]][0])
            self.aa = aipy.cal.get_aa(self.calfile, uv['sdf'], uv['sfreq'], uv['nchan'])
        else:
            uv = UVData()
            uv.read_miriad(self.files[self.files.keys()[0]][0])
            # convert from Hz -> GHz
            freqs = uv.freq_array[0, :] / 1e9
            self.aa = hera_cal.utils.get_aa_from_uv(uv, freqs)
        del uv
        self.aa.set_active_pol(self.pol)

        for baselength in sorted(self.caldata[0]):
            self.vis_sq_bl = np.zeros((ntimes//2, nchan), dtype=np.complex128)
            for antpair in self.caldata[0][baselength]:
                self.vis_sq_antpair = np.zeros((ntimes//2, nchan), dtype=np.complex128)
                self.cleanfft(antpair)
                self.phasing(ntimes, antpair)
                self.vis_sq_bl += self.vis_sq_antpair
            self.vis_sq_bl /= len(self.caldata[0][baselength])
            self.cwedgeslices.append(self.vis_sq_bl)

    def form_times(self):
        self.lst_range_str = np.unique(self.lst_range_str)
        self.lst_range_num = np.unique(self.lst_range_num)
        self.lst_range = np.vstack((self.lst_range_num, self.lst_range_str))
        self.times = np.vstack((self.lst_range, self.info['times']))

    def form_delays(self):
        channel_width = (self.info['freqs'][1] - self.info['freqs'][0]) * (10**3)
        num_bins = len(self.info['freqs'])
        self.delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))

    def savenpz(self):
        np.savez(
            self.npz_name,
            pol=self.pol,
            caldata=self.caldata,
            wslices=self.wedgeslices,
            cwslices=self.cwedgeslices,
            aslices=self.antpairslices,
            delays=self.delays,
            times=self.times,
            bl_length=self.bl_length,
            hist=self.history)


# Create Wedge from Pitchfork
def wedge_delayavg(npz_name, multi=False):

    plot_data = np.load(npz_name)
    delays, wedgevalues, baselines = plot_data['dlys'], plot_data['cldt'][3], plot_data['cldt'][3]
    # d_start = plot_data['dlys'][0]
    # d_end = plot_data['dlys'][-1]
    split = (len(wedgevalues[0, :])/2)

    wedgevalues2 = np.zeros((len(wedgevalues), len(delays)))

    for baselength in range(len(wedgevalues)):
        for i in range(split):
            avg = ((wedgevalues[baselength, (split-1+i)]+wedgevalues[baselength, (split+i)])/2)
            wedgevalues2[baselength][split-i] = avg
    delayavg_wedgevalues = wedgevalues2.T.T.T
    npz_delayavg = (npz_name[:-11] + 'delayavg.npz')
    np.savez(npz_delayavg, wdgslc=delayavg_wedgevalues, dlys=delays, bls=baselines)
    # saving to longer arrary fml??? idk
    print("got here!!!")
    return npz_delayavg


# Plotting Functions:
def plot_timeavg(npz_name, path):
    # load and format data
    data = np.load(npz_name)
    npz_name = npz_name.split('/')[-1]
    delays = data['delays']
    wedgeslices = data['wslices']
    caldata = data['caldata']
    times = data['times']

    # sets max/min values to plot depending on if its simulated data
    if 'SIM' in npz_name:
        vmax = -2
        vmin = -12
    else:
        vmax = -1
        vmin = -5

    # format data so y axis properly scaled
    plotindeces = [int(round(i*10)) for i in caldata[3]]
    plotdata = np.zeros((plotindeces[-1], wedgeslices.shape[-1]), dtype=np.float64)
    j = 0
    for i in range(len(plotindeces)):
        plotdata[j:plotindeces[i]] = wedgeslices[i]
        j = plotindeces[i]

    # plot the data
    plt.imshow(
        plotdata,
        aspect='auto',
        interpolation='nearest',
        extent=[delays[0], delays[-1], plotindeces[-1], 0],
        vmin=vmin,
        vmax=vmax)

    # label colorbar and axes, format axes
    plt.colorbar().set_label(r'$\log_{10}({\rm mK}^2)$')
    plt.xlabel("Delay [ns]")
    plt.xlim((-450, 450))
    plt.ylabel("Baseline Length [m]")
    plt.yticks(plotindeces, [round(n, 1) for n in caldata[3]])

    # set titles
    plt.suptitle("JD: {JD}; LST {start} to {end}".format(JD=npz_name.split('.')[1], start=times[1][0][:-6], end=times[1][-1][:-6]))
    plt.title(npz_name.split('.')[3])

    # plot center line
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

    # plot horizon lines
    horizons = []
    for length in caldata[3]:
        horizons.append(length / sc.c * 10**9)
    j = 0
    for i in range(len(horizons)):
        x1, y1 = [horizons[i], horizons[i]], [j, plotindeces[i]]
        x2, y2 = [-horizons[i], -horizons[i]], [j, plotindeces[i]]
        plt.plot(x1, y1, x2, y2, color='white', linestyle='--', linewidth=.75)
        j = plotindeces[i]

    # save and close
    plt.savefig(path + npz_name[:-4] + '.png')
    plt.show()
    plt.close()
    plt.clf()


def plot_timeavg_multi(npz_names, path):
    nplots = len(npz_names)

    # Form subplots
    f, axes = plt.subplots(1, nplots, sharex=True, sharey=True, figsize=(12, 5))

    # Loop through each subplot and plot
    for i, ax in enumerate(axes):
        # Get data from npz_files
        npz_name = npz_names[i]
        data = np.load(npz_name)
        npz_name = npz_name.split('/')[-1]
        delays = data['delays']
        wedgeslices = data['wslices']
        caldata = data['caldata']
        times = data['times']

        ax.set_title(npz_name.split('.')[3])

        if 'SIM' in npz_name:
            vmax = -2
            vmin = -12
        else:
            vmax = -1
            vmin = -5

        # format data so y axis properly scaled
        plotindeces = [int(round(i*10)) for i in caldata[3]]
        plotdata = np.zeros((plotindeces[-1], wedgeslices.shape[-1]), dtype=np.float64)
        j = 0
        for i in range(len(plotindeces)):
            plotdata[j:plotindeces[i]] = wedgeslices[i]
            j = plotindeces[i]

        # plot the data
        plot = ax.imshow(
            plotdata,
            aspect='auto',
            interpolation='nearest',
            extent=[delays[0], delays[-1], plotindeces[-1], 0],
            vmin=vmin,
            vmax=vmax)

        # Plot center line to easily see peak offset
        ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

        # plot horizon lines
        horizons = []
        for length in caldata[3]:
            horizons.append(length / sc.c * 10**9)
        j = 0
        for i in range(len(horizons)):
            x1, y1 = [horizons[i], horizons[i]], [j, plotindeces[i]]
            x2, y2 = [-horizons[i], -horizons[i]], [j, plotindeces[i]]
            ax.plot(x1, y1, x2, y2, color='white', linestyle='--', linewidth=.75)
            j = plotindeces[i]

        del npz_name

    plt.suptitle("JD: {JD}; LST {start} to {end}".format(JD=npz_names[0].split('/')[-1].split('.')[1], start=times[1][0][:-6], end=times[1][-1][:-6]))

    # Smooshes the plots together
    f.subplots_adjust(wspace=0)

    plt.xlim((-450, 450))
    plt.yticks(plotindeces, [round(n, 1) for n in caldata[3]])

    # X and Y axis labels
    f.text(0.5, 0.025, 'Delay [ns]', ha='center')
    f.text(0.075, 0.5, 'Baseline Length [m]', va='center', rotation='vertical')

    # Colorbar
    cbar_ax = f.add_axes([0.9125, 0.25, 0.025, 0.5])
    cbar = f.colorbar(plot, cax=cbar_ax)
    cbar.set_label(r'$\log_{10}({\rm mK}^2)$')

    # Naming and saving
    npz_names = [npz_name.split('/')[-1] for npz_name in npz_names]
    pols = [npz_name.split('.')[3] for npz_name in npz_names]
    pols = ["".join(pols)]
    f1, f2 = npz_names[0].split('/')[-1].split('.')[:3], npz_names[0].split('/')[-1].split('.')[4:-1]
    npz_name = ".".join(f1 + pols + f2)
    plt.savefig(path + npz_name + '.png')

    # plt.show()
    # plt.close()
    plt.clf()


def plot_blavg(npz_name, path):
    data = np.load(npz_name)
    npz_name = npz_name.split('/')[-1]
    delays = data['delays']
    wedgeslices = data['wslices']
    caldata = data['caldata']
    times = data['times']

    # Format time scale in minutes
    time_start = (float(times[2][-1]) - float(times[2][0])) * 24*60
    time_end = 0

    # create subplot to plot data
    f, axes = plt.subplots(len(wedgeslices), 1, sharex=True, sharey=True)

    # Create subplot to make titles and lables look right
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel("Delay [ns]")
    plt.ylabel("Time [min]")
    plt.suptitle("JD: {JD}; LST {start} to {end}".format(JD=npz_name.split('.')[1], start=times[1][0][:-6], end=times[1][-1][:-6]))
    plt.title(npz_name.split('.')[3])

    # Loop through each wedgeslice, and imshow each in its own subplot
    for i in range(len(wedgeslices)):
        im = axes[i].imshow(
            wedgeslices[i],
            aspect='auto',
            interpolation='nearest',
            extent=[delays[0], delays[-1], time_start, time_end],
            vmax=-1,
            vmin=-5)

        horizon = caldata[3][i] / sc.c * 10**9
        x1, y1 = [horizon, horizon], [time_start, time_end]
        x2, y2 = [-horizon, -horizon], [time_start, time_end]
        axes[i].plot(x1, y1, x2, y2, color='white')

    cax, kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    cbar = plt.colorbar(im, cax=cax, **kw)
    cbar.set_label(r'$\log_{10}({\rm mK}^2)$')

    axes[0].set_xlim(-450, 450)

    plt.savefig(path + npz_name[:-4] + '.png')
    plt.show()


def plot_flavors(npz_name, path):
    data = np.load(npz_name)
    npz_name = npz_name.split('/')[-1]
    delays = data['delays']
    wedgeslices = data['wslices']
    caldata = data['caldata']
    times = data['times']

    if 'SIM' in npz_name:
        vmax = -2
        vmin = -12
    else:
        vmax = -1
        vmin = -5

    plt.imshow(
        wedgeslices,
        aspect='auto',
        interpolation='nearest',
        extent=[delays[0], delays[-1], len(wedgeslices), 0],
        vmin=vmin,
        vmax=vmax)

    plt.colorbar().set_label(r'$\log_{10}({\rm mK}^2)$')

    plt.xlabel("Delay [ns]")
    plt.xlim((-450, 450))

    plt.ylabel("Baseline Length [m]: Orientation")

    plt.suptitle("JD: {JD}; LST {start} to {end}".format(JD=npz_name.split('.')[1], start=times[1][-1][:-6], end=times[1][-1][:-6]))
    plt.title(npz_name.split('.')[3])

    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

    ticks = []
    slopedict = caldata[1]
    for baseline in sorted(slopedict.keys()):
        for slope in sorted(slopedict[baseline].keys()):
            ticks.append("{:.3}: {:8.3}".format(baseline, slope))
    plt.yticks(np.arange(len(ticks)), ticks)

    horizons = []
    for baseline in sorted(slopedict.keys()):
        for slope in sorted(slopedict[baseline].keys()):
            horizons.append(baseline / sc.c * 10**9)

    for j in range(len(horizons)):
        x1, y1 = [horizons[j], horizons[j]], [j, j+1]
        x2, y2 = [-horizons[j], -horizons[j]], [j, j+1]
        plt.plot(x1, y1, 'w', x2, y2, 'w')

    plt.savefig(path + npz_name[:-4] + '.png')
    plt.show()


def plot_flavors_multi(npz_names, path):
    nplots = len(npz_names)

    # Form subplots
    f, axes = plt.subplots(1, nplots, sharex=True, sharey=True, figsize=(12, 5))

    # Loop through each subplot and plot
    for i, ax in enumerate(axes):
        npz_name = npz_names[i]
        data = np.load(npz_names[i])
        npz_name = npz_name.split('/')[-1]
        delays = data['delays']
        wedgeslices = data['wslices']
        caldata = data['caldata']
        times = data['times']

        ax.set_title(npz_name.split('.')[3])

        if 'SIM' in npz_name:
            vmax = -2
            vmin = -12
        else:
            vmax = -1
            vmin = -5

        plot = ax.imshow(
            wedgeslices,
            aspect='auto',
            interpolation='nearest',
            extent=[delays[0], delays[-1], len(wedgeslices), 0],
            vmin=vmin,
            vmax=vmax)

        ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

        horizons = []
        for baseline in sorted(caldata[1].keys()):
            for slope in sorted(caldata[1][baseline].keys()):
                horizons.append(baseline / sc.c * 10**9)

        for j in range(len(horizons)):
            x1, y1 = [horizons[j], horizons[j]], [j, j+1]
            x2, y2 = [-horizons[j], -horizons[j]], [j, j+1]
            ax.plot(x1, y1, 'w', x2, y2, 'w')

    plt.suptitle("JD: {JD}; LST {start} to {end}".format(JD=npz_name.split('.')[1], start=times[1][0][:-6], end=times[1][-1][:-6]))

    # Smooshes the plots together
    f.subplots_adjust(wspace=0)

    plt.xlim((-450, 450))
    ticks = []
    for baseline in sorted(caldata[1].keys()):
        for slope in sorted(caldata[1][baseline].keys()):
            ticks.append("{:.3}: {:8.3}".format(baseline, slope))
    plt.yticks(np.arange(len(ticks)), ticks)

    # X and Y axis labels
    f.text(0.5, 0.025, 'Delay [ns]', ha='center')
    f.text(0.025, 0.5, 'Baseline Length [m]: Orientation', va='center', rotation='vertical')

    # Colorbar
    cbar_ax = f.add_axes([0.9125, 0.25, 0.025, 0.5])
    cbar = f.colorbar(plot, cax=cbar_ax)
    cbar.set_label(r'$\log_{10}({\rm mK}^2)$')

    # Naming and saving
    pols = [npz_name.split('.')[3] for npz_name in npz_names]
    pols = ["".join(pols)]
    f1, f2 = npz_names[0].split('/')[-1].split('.')[:3], npz_names[0].split('/')[-1].split('.')[4:-1]
    npz_name = ".".join(f1 + pols + f2)
    print(path + npz_name + '.png')
    plt.savefig(path + npz_name + '.png')

    plt.show()


def plot_bltype(npz_name):
    plot_data = np.load(npz_name)

    d_start = plot_data['delays'][0]
    d_end = plot_data['delays'][-1]
    t_start = plot_data['aslices'][0].shape[0]
    length = float(plot_data['bl_length'])

    # create subplot to plot data
    f, axarr = plt.subplots(len(plot_data['aslices']), 1, sharex=True, sharey=True)

    # add axes labels
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel("Delay (ns)")
    plt.ylabel("Time", labelpad=15)

    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2]+'.'+npz_name.split('.')[3]+'.baseline'+npz_name.split('.')[7])
    lengthstr = str(plot_data['bl_length'])
    plt.title('baseline length:'+lengthstr)

    # plot individual wedge slices
    for i in range(len(plot_data['aslices'])):
        # plot the graph
        im = axarr[i].imshow(plot_data['aslices'][i], aspect='auto', interpolation='nearest', vmin=-9, vmax=1, extent=[d_start, d_end, t_start, 0])
        # plot light delay time lines
        light_time = (plot_data['bl_length'])/sc.c*10**9
        x1, y1 = [light_time, light_time], [0, np.shape(plot_data['aslices'][i])[0]]
        x2, y2 = [-light_time, -light_time], [0, np.shape(plot_data['aslices'][i])[0]]
        axarr[i].plot(x1, y1, x2, y2, color='white')
        axarr[i].set_ylabel(plot_data['caldata'][0][length][i], fontsize=6)

    cax, kw = mpl.colorbar.make_axes([ax for ax in axarr.flat])
    plt.colorbar(im, cax=cax, **kw)

    # scale x axis to the significant information
    axarr[0].set_xlim(-450, 450)

    f.set_size_inches(6, 9, forward=True)
    plt.savefig(npz_name[:-3] + 'png')
    plt.show()


def plot_delayavg(npz_delayavg):

    plot_data = np.load(npz_delayavg)
    wedgevalues, baselines = plot_data['wslices'], plot_data['caldata'][3]
    d_start = plot_data['delays'][0]
    d_end = plot_data['delays'][-1]
    plt.imshow(wedgevalues, aspect='auto', interpolation='nearest', extent=[0, len(npz_delayavg), d_start, d_end], vmin=-3.0, vmax=1.0)
    # plot = plt.imshow(npz_delayavg, aspect='auto', interpolation='nearest',extent=[0,len(wedgevalues),d_start,d_end], vmin=-3.0, vmax=1.0)

    plt.xlabel("Baseline length (short to long)")
    plt.ylabel("Delay (ns)")
    cbar = plt.colorbar()
    cbar.set_label("log10((mK)^2)")
    plt.xlim((0, len(baselines)))
    plt.ylim(0, 450)
    plt.title(npz_delayavg.split('.')[1]+'.'+npz_delayavg.split('.')[2]+'.'+npz_delayavg.split('.')[3])

    # calculate light travel time for each baselength
    light_times = []
    for length in plot_data['caldata'][3]:
        light_times.append(length/sc.c*10**9)

    # plot lines on plot using the light travel time
    for i in range(len(light_times)):
        y1, x1 = [light_times[i], light_times[i]], [i, i+1]
        y2, x2 = [-light_times[i], -light_times[i]], [i, i+1]
        plt.plot(x1, y1, x2, y2, color='white')

    print("got here1")
    plt.savefig(npz_delayavg[:-12]+'delayavg.png')
    plt.show()


# 1D Plotting for timeavg wedges
def plot_1D(npz_name, baselines=[]):

    """
    Plots all baselines overlapped on a 1D plot.
    If baselines is a specified argument (start indexing with baseline length #1),
    then only plots the provided baselines.
    """

    plot_data = np.load(npz_name)

    if len(baselines):
        baselines = [i-1 for i in baselines]
    else:
        baselines = range(len(plot_data['wslices']))

    plt.figure(figsize=(12, 6))
    G = gridspec.GridSpec(2, 9)

    plt.subplot(G[:, 0:4])
    for i in baselines:
        plt.plot(plot_data['delays'], plot_data['wslices'][i], label='bl len '+str(plot_data['caldata'][3][i]))
    plt.xlabel('Delay (ns)')
    plt.ylabel('log10((mK)^2)')
    plt.legend(loc='upper left')
    plt.ylim((-3.5, 2.0))

    plt.subplot(G[:, 5:9])
    for i in baselines:
        plt.plot(plot_data['delays'], plot_data['wslices'][i])
    if len(baselines) == 1:
        light_time = plot_data['caldata'][3][baselines[0]]/sc.c*10**9
        plt.axvline(light_time, color='#d3d3d3', linestyle='--')
        plt.axvline(-1*light_time, color='#d3d3d3', linestyle='--')
        plt.axvline(0, color='#d3d3d3', linestyle='--')
    plt.xlim((-450, 450))
    plt.ylim((-3.5, 2.0))

    plt.xlabel('Delay (ns)')
    plt.ylabel('log10((mK)^2)')
    npz_name = npz_name.split('/')[-1]
    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2]+'.'+npz_name.split('.')[3])

    plt.show()


def plot_multi_1D(npz_names, baselines=[]):
    """Plots four 1D plots next to each other.
    If baselines is a specified argument (start indexing with baseline lengt #1),
    then only plots the the provided baselines."""

    plot_data = np.load(npz_names[0])
    # set up baselines
    if len(baselines):
        baselines = [i-1 for i in baselines]
    else:
        baselines = range(len(plot_data['wslices']))

    # set up the plotting space
    plt.figure(figsize=(18, 4))
    G = gridspec.GridSpec(1, 4)

    # plot each 1D plot
    polorder = ''
    for n in range(len(npz_names)):

        # load data, format plotting section
        plot_data = np.load(npz_names[n])
        plt.subplot(G[:, n:n+1])

        # plot the data
        for i in baselines:
            plt.plot(plot_data['delays'], plot_data['wslices'][i], label='bl len '+str(plot_data['caldata'][3][i]))
        if len(baselines) == 1:
            light_time = plot_data['caldata'][3][baselines[0]]/sc.c*10**9
            plt.axvline(light_time, color='#d3d3d3', linestyle='--')
            plt.axvline(-1*light_time, color='#d3d3d3', linestyle='--')
        plt.axvline(0, color='#d3d3d3', linestyle='--')
        plt.xlim((-450, 450))
        plt.ylim((-3.0, 2.0))

        if n == 0:
            plt.legend(loc='upper left')
        plt.xlabel('Delay (ns)')
        plt.ylabel('log10((mK)^2)')
        pol = npz_names[n].split('/')[-1].split('.')[3]
        plt.title(pol)
        polorder += pol

    npz_name = npz_names[0].split('/')[-1]
    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2])

    if len(baselines) == len(plot_data['wslices']):
        blstr = 'allbls'
    else:
        blstr = 'bl'
        for bl in baselines:
            blstr += str(bl+1)

    plt.tight_layout()
    plt.savefig(npz_name.split(polorder[0])[0]+polorder+npz_name.split(polorder[0])[-1][:-3] + "multi1D." + blstr + ".png")
    plt.show()


# Inside/Outside Averaging
def in_out_avg(npz_name):
    data = np.load(npz_name)
    history = data['hist'].tolist()
    num_files = len(history['filenames'])

    light_times = []
    for length in data['caldata'][3]:
        light_times.append(length / (sc.c * (10**9)))

    total_in, total_out = 0, 0
    total_in_count, total_out_count = 0, 0
    for i in range(len(data['caldata'][3])):
        for index, delay in enumerate(data['delays']):
            if abs(delay) >= light_times[i]:
                total_out += data['wslices'][i][index]
                total_out_count += 1
            else:
                total_in += data['wslices'][i][index]
                total_in_count += 1

    avg_in = total_in / total_in_count
    avg_out = total_out / total_out_count

    return (avg_in, avg_out, num_files)


def plot_avgs(npz_names, rng):
    name = ".".join(npz_names[-1].split('/')[-1].split('.')[:7])

    total_files = []
    avgs_in = []
    avgs_out = []
    for npz_name in npz_names:
        avgs_in_out = in_out_avg(npz_name)
        total_files.append(avgs_in_out[2])
        avgs_in.append(avgs_in_out[0])
        avgs_out.append(avgs_in_out[1])

    plt.figure(figsize=(15, 10))

    plot_avgs_out = plt.scatter(total_files, avgs_out)
    plot_avgs_in = plt.scatter(total_files, avgs_in)
    plt.legend((plot_avgs_out, plot_avgs_in), ('Averages Outside Wedge', 'Averages Inside Wedge'))

    plt.xlim(0, len(total_files))
    plt.xticks(np.arange(0, len(total_files)+4, 2))
    plt.xlabel('Number of Files Used in Averaging')

    plt.ylim(-3.5, 1.5)
    plt.yticks(np.arange(-3, 2))
    plt.ylabel('log10((mK)^2)')

    plt.title(name)
    plt.suptitle('How Well does the Pitchfork Average Down?')

    # Setting range to calculate a fit over.
    # rng = (7, 26)
    x = np.arange(rng[0], rng[1])

    # Outside RMSE and fit calculation.
    m, b = np.polyfit(x, avgs_out[rng[0]:rng[1]], 1)
    fit = np.poly1d([m, b])
    outside = avgs_out[rng[0]:rng[1]]

    RMSE = gu.RMSE(fit(x), outside)
    plt.plot(x, m*x + b, '--', color='k')
    plt.text(8, 0, 'Outside Pitchfork Fit RMSE: {:.3}\nOutside Pitchfork Fit Slope: {:.3}'.format(RMSE, m))

    # Inside RMSE and fit calculation.
    m, b = np.polyfit(x, avgs_in[rng[0]:rng[1]], 1)
    fit = np.poly1d([m, b])
    inside = avgs_in[rng[0]:rng[1]]

    RMSE = gu.RMSE(fit(x), inside)
    plt.plot(x, m*x + b, '--', color='k')
    plt.text(8, -3, 'Inside Pitchfork Fit RMSE: {:.3}\nInside Pitchfork Fit Slope: {:.3}'.format(RMSE, m))

    plt.savefig(name + ".avg_val.png")
    plt.show()
