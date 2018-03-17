"""
This module contains functions to create and plot wedges/pitchforks from HERA data files.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
Paul La Plante <plaplant_at_sas.upenn.edu>
"""
# Python Standard Modules
import os

# Python Community Science Modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
import astropy.cosmology as cos

# HERA Collaboration Modules
from pyuvdata import UVData
import pyuvdata.utils as uvutils
import aipy
import hera_cal

# For Interactive Development
from IPython import embed

class Eris(object):
    def __init__(self, Zeus, pol):
        self.Zeus = Zeus
        self.pol = pol

        self.npz_name = str()

        self.caldata = dict()

        self.uvd = None
        self.uvd2 = None

        self.aa = None

        self.info = dict()
        self.data = dict()
        self.flags = dict()

        self.vis_sq_bl = None
        self.vis_sq_slope = None
        self.vis_sq_antpair = None

        self.bl_slices = list()
        self.sl_slices = list()
        self.ant_slices = list()

    def name_npz(self):
        zen = 'zen'
        HH = 'HH'

        extension = self.Zeus.ext
        pol_naming = self.Zeus.catalog[extension].keys()[0]

        data_file0 = self.Zeus.catalog[extension][pol_naming][0]
        data_filef = self.Zeus.catalog[extension][pol_naming][-1]

        file0 = data_file0['file']
        filef = data_filef['file']

        JD = int(data_file0['JD'])
        JDT0 = data_file0['JD'] - JD
        JDTf = data_filef['JD'] - JD

        JD = str(JD)
        JDT0 = str(JDT0)[2:7]
        JDTf = str(JDTf)[2:7]

        num_files = np.unique(self.Zeus.catalog[extension][pol_naming]['file']).shape[0]
        JDT = '{JDT0}_{num_files}_{JDTf}'.format(JDT0=JDT0, num_files=num_files, JDTf=JDTf)

        pol = self.pol

        if self.Zeus.exants:
            exants = [str(ant) for ant in self.Zeus.exants]
            exants = '_'.join(exants)
        else:
            exants = 'None'

        freqrange = '{start}_{end}'.format(start=self.Zeus.freqrange[0], end=self.Zeus.freqrange[1])

        npz_name = [zen, JD, JDT, pol, exants, freqrange, HH, extension, 'npz']
        npz_name = '.'.join(npz_name)

        self.npz_name = os.path.join(self.Zeus.path, npz_name)
        print(self.npz_name)

    def load_MIRIAD(self):
        """Formats data and flags array for a specific polarization, also apply flags."""
        # Check what type of polarization was specified
        if self.Zeus.pol_type == 'standard':
            self.uvd = UVData()
            self.uvd.read_miriad(self.Zeus.files[self.pol])
        elif self.pol == 'I' or self.pol == 'Q':
            self.uvd = UVData()
            self.uvd.read_miriad(self.Zeus.files['xx'])
            self.uvd2 = UVData()
            self.uvd2.read_miriad(self.Zeus.files['yy'])
        elif self.pol == 'U' or self.pol == 'V':
            self.uvd = UVData()
            self.uvd.read_miriad(self.Zeus.files['xy'])
            self.uvd2 = UVData()
            self.uvd2.read_miriad(self.Zeus.files['yx'])

        self.calculate_caldata()

        # Get metadata and convert freqs array from Hz -> GHz
        self.info['freqs'] = self.uvd.freq_array[0, :] / 1e9
        self.info['times'] = np.unique(self.uvd.time_array)

        # Extract and format data and flags arrays
        for key, d in self.uvd.antpairpol_iter():
            # Ignore auto-correlation
            if key[0] == key[1]:
                continue
            ind1, ind2, ipol = self.uvd._key2inds(key)
            for ind in [ind1, ind2]:
                if len(ind) == 0:
                    continue

                antkey = key[:2]
                if self.Zeus.pol_type == 'standard':
                    f = self.uvd.flag_array[ind, 0, :, ipol]
                    polkey = key[2].lower()
                    self.data[antkey] = {polkey: d}
                    self.flags[antkey] = {polkey: f}
                elif self.Zeus.pol_type == 'stokes':
                    f = self.uvd.flag_array[ind, 0, :, 0]
                    d2 = self.uvd2.data_array[ind, 0, :, 0]
                    f2 = self.uvd2.flag_array[ind, 0, :, 0]

                    if self.pol == 'I':
                        self.data[antkey] = {'I': d + d2}
                        self.flags[antkey] = {'I': f + f2}
                    elif self.pol == 'Q':
                        self.data[antkey] = {'Q': d - d2}
                        self.flags[antkey] = {'Q': f + f2}
                    elif self.pol == 'U':
                        self.data[antkey] = {'U': d + d2}
                        self.flags[antkey] = {'U': f + f2}
                    elif self.pol == 'V':
                        self.data[antkey] = {'V': -1j*d + 1j*d2}
                        self.flags[antkey] = {'V': f + f2}

        del self.uvd2

        # Apply flags to data by looking for where the flags array is true, and zeroing the corresponding data elements.
        for pair in self.caldata['pairs']:
            self.data[pair][self.pol] = np.where(self.flags[pair][self.pol], 0. + 0.*1j, self.data[pair][self.pol])

        for i in self.flags:
            for j in self.flags[i]:
                    self.flags[i][j] = np.logical_not(self.flags[i][j]).astype(int)

    def calculate_caldata(self):
        """Returns a dictionary of baseline lengths and the corresponding pairs.
        The data is based on a calfile. ex_ants is a list of integers that
        specify antennae to be exlcuded from calculation.

        Requires cal file to be in PYTHONPATH."""
        if self.Zeus.calfile:
            try:
                print('Reading calfile: %s...' %self.Zeus.calfile)
                exec('import %s as cal' %self.Zeus.calfile, globals(), locals())
                antennae = cal.prms['antpos_ideal']
            except ImportError:
                raise Exception('Unable to import: %s' %self.Zeus.calfile)
        elif self.uvd:
            # Build antenna positions from data file itself.
            print('Generating calibration information from MIRIAD file.')
            antennae = {}
            lat, lon, alt = self.uvd.telescope_location_lat_lon_alt
            for i, antnum in enumerate(self.uvd.antenna_numbers):
                pos = self.uvd.antenna_positions[i, :] + self.uvd.telescope_location
                xyz = uvutils.ENU_from_ECEF(pos, latitude=lat, longitude=lon, altitude=alt)
                antennae[antnum] = {'top_x': xyz[0], 'top_y': xyz[1], 'top_z': xyz[2]}
        else:
            raise Exception('UVData object does not exist. Try supplying a calfile.')

        # Check that self.Zeus.exants contain only antennae that exist
        if not set(self.Zeus.exants).issubset(set(antennae.keys())):
            raise Exception('You provided invalid antenna(e) to exclude.')

        # Remove all placeholder antennae from consideration
        # Remove all antennae from exants from consideration
        ants = []
        for ant in antennae.keys():
            if (not antennae[ant]['top_z'] == -1) and (ant not in self.Zeus.exants):
                ants.append(ant)
        ants = np.array(ants)

        # Set the bins for binning the baselines
        bins_baseline = np.arange(0, 750, self.Zeus.BIN_WIDTH)
        
        # Initialize baselines, slopes, pairs arrays
        baselines = np.array([], dtype=np.float128)
        slopes = np.array([], dtype=np.float128)
        pairs = {}

        # Cycle through antennae to create pairs
        for ant_i in ants:
            for ant_j in ants:
                if ant_i >= ant_j:
                    continue
                pair = (ant_i, ant_j)

                # Find the baseline length of the pair, bin it, format it, and add it to the baselines array
                dx = antennae[pair[1]]['top_x'] - antennae[pair[0]]['top_x']
                dy = antennae[pair[1]]['top_y'] - antennae[pair[0]]['top_y']
                baseline = np.sqrt(np.power(dx, 2) + np.power(dy, 2))
                baseline = np.round(baseline, decimals=2)
                baseline = np.digitize(baseline, bins_baseline) * self.Zeus.BIN_WIDTH
                baselines = np.append(baselines, baseline) 

                # Find the slope of the pair, format it, and add it to the slopes array
                dy = antennae[pair[1]]['top_y'] - antennae[pair[0]]['top_y']
                dx = antennae[pair[1]]['top_x'] - antennae[pair[0]]['top_x']
                if dx != 0:
                    slope = float(np.round(dy / dx, decimals=2))
                else:
                    slope = np.inf
                slopes = np.append(slopes, slope)

                # Add the pair and its slope and baseline to the pairs dictionary: {(ant_1, ant_2): (baseline, slope), ...}
                pairs[pair] = (np.float(np.round(baseline, 1)), np.float(np.round(slope, 2)))

        # Sort and remove duplicates from baselines and slopes
        baselines = np.unique(np.sort(baselines))
        slopes = np.unique(np.sort(slopes))

        # Sort pairs into antdict and slopedict
        antdict = {}
        slopedict = {}
        for pair, (baseline, slope) in pairs.items():
            if baseline in antdict:
                antdict[baseline].append(pair)
            else:
                antdict[baseline] = [pair]

            if baseline in slopedict:
                if slope in slopedict[baseline]:
                    slopedict[baseline][slope].append(pair)
                else:
                    slopedict[baseline][slope] = [pair]
            else:
                slopedict[baseline] = {slope: []}
                slopedict[baseline][slope].append(pair)

        self.caldata = {'antdict': antdict, 'slopedict': slopedict, 'pairs': pairs, 'baselines': baselines, 'slopes': slopes}
    
    def pitchfork(self):
        self.info['freqs'] = self.info['freqs'][self.Zeus.freqrange[0]: self.Zeus.freqrange[1]]
        for antpair in self.data.keys():
            self.data[antpair][self.pol] = self.data[antpair][self.pol][:, self.Zeus.freqrange[0]: self.Zeus.freqrange[1]]
            self.flags[antpair][self.pol] = self.flags[antpair][self.pol][:, self.Zeus.freqrange[0]: self.Zeus.freqrange[1]]

        ntimes = len(self.info['times'])
        nchan = len(self.info['freqs'])

        if self.Zeus.calfile:
            self.aa = hera_cal.utils.get_aa_from_calfile(self.info['freqs'], self.Zeus.calfile)
        else:
            self.aa = hera_cal.utils.get_aa_from_uv(self.uvd, self.info['freqs'])
        self.aa.set_active_pol(self.pol)

        for baseline in sorted(self.caldata['slopedict']):
            self.vis_sq_bl = np.zeros((ntimes // 2, nchan), dtype=np.complex128)
            for slope in sorted(self.caldata['slopedict'][baseline]):
                self.vis_sq_slope = np.zeros((ntimes // 2, nchan), dtype=np.complex128)
                for antpair in sorted(self.caldata['slopedict'][baseline][slope]):
                    self.vis_sq_antpair = np.zeros((ntimes // 2, nchan), dtype=np.complex128)

                    self.fft_clean_phase(antpair, ntimes, self.Zeus.CLEAN)
                    
                    self.ant_slices.append(self.vis_sq_antpair)

                    self.vis_sq_slope += self.vis_sq_antpair
                    self.vis_sq_bl += self.vis_sq_antpair
            
                self.vis_sq_slope /= len(self.caldata['slopedict'][baseline][slope])
                self.sl_slices.append(self.vis_sq_slope)
            
            self.vis_sq_bl /= len(self.caldata['antdict'][baseline])
            self.bl_slices.append(self.vis_sq_bl)

    def fft_clean_phase(self, antpair, ntimes, CLEAN):
        # FFT
        w = aipy.dsp.gen_window(self.data[antpair][self.pol].shape[-1], window='blackman-harris')
        _dw = np.fft.ifft(self.data[antpair][self.pol] * w)

        # CLEAN
        _ker = np.fft.ifft(self.flags[antpair][self.pol] * w)
        gain = aipy.img.beam_gain(_ker)
        for time in range(_dw.shape[0]):
            _dw[time, :], info = aipy.deconv.clean(_dw[time, :], _ker[time, :], tol=CLEAN)
            _dw[time, :] += info['res'] / gain

        fft_2Ddata = np.ma.array(_dw)

        # Phase
        for i in range(ntimes):
            if i != 0:
                old_zenith = zenith
            time = self.info['times'][i]
            self.aa.set_jultime(time)
            lst = self.aa.sidereal_time()
            zenith = aipy.phs.RadioFixedBody(lst, self.aa.lat)
            zenith.compute(self.aa)

            if i % 2:
                v1 = fft_2Ddata[i - 1, :]
                phase_correction = np.conj(self.aa.gen_phs(zenith, antpair[0], antpair[1])) * self.aa.gen_phs(old_zenith, antpair[0], antpair[1])
                v2 = fft_2Ddata[i, :] * phase_correction
                self.vis_sq_antpair[i // 2, :] = np.conj(v1) * v2

    def save(self):
        channel_width = np.diff(self.info['freqs'])[0]
        num_bins = self.info['freqs'].shape[0]
        delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width))

        np.savez(self.npz_name,
                 Zeus=self.Zeus,
                 bl=self.bl_slices,
                 sl=self.sl_slices,
                 ant=self.ant_slices,
                 caldata=self.caldata,
                 pol=self.pol,
                 freqs=self.info['freqs'],
                 delays=delays)

class Ares(object):
    def __init__(self, Zeus):
        self.Zeus = Zeus
        self.fontsize = 20
        self.xaxis_x, self.xaxis_y = 0.5, 0.03
        self.yaxis_x, self.yaxis_y = 0.07, 0.5
        self.p_axes = [0.9125, 0.25, 0.025, 0.5]
        self.delay_range_multipler = 1. / 2.

    def makePlot(self):
        nplots = len(self.Zeus.inputfiles)
        f, axes = plt.subplots(1, nplots, sharex=True, sharey=True, figsize=(20, 10))
        if nplots == 1:
            axes = [axes]
        f.subplots_adjust(wspace=0)
        self.pols = []
        
        for k, self.ax in enumerate(axes):
            file0 = self.Zeus.inputfiles[k]
            data = np.load(file0) 
            pol = data['pol'].tolist()
            self.fileZeus = data['Zeus'].tolist()
            caldata = data['caldata'].tolist()
            
            self.pols.append(pol)
            self.freqrange = np.array(self.fileZeus.freqrange, dtype=np.int32)
            self.baselines = caldata['baselines']
            self.delays = data['delays']

            self.ax.set_title(pol, fontsize=self.fontsize)

            if self.Zeus.tag_avg == "timeavg":
                if self.Zeus.tag_depth == "bl":
                    self.wedgeslices = np.log10(np.fft.fftshift(np.abs(np.nanmean(data['bl'], axis=1)), axes=1))
                    self.kpr_index()
                    self.unit()
                    self.horizon()
                    if self.Zeus.tag_wedge == "1d":
                        plt.sca(axes[k])
                        j = 0
                        for i, bl in enumerate(self.baselines):
                            plt.plot(
                                self.delays,
                                self.wedgeslices[j],
                                label=self.axis_kpr_b_1d[i])
                            j = self.kpr_indices[i]
                        f.text(self.xaxis_x, self.xaxis_y, self.axis_kpl_tau, fontsize=self.fontsize, ha='center')
                        f.text(self.yaxis_x, self.yaxis_y, self.axis_power, fontsize=self.fontsize, va='center', rotation='vertical')
                        plt.ylim((self.vmin, self.vmax))
                        plt.xlim((-5000 * self.kpl * self.delay_range_multipler, 5000 * self.kpl * self.delay_range_multipler))
                        plt.xticks(self.label_kpl_tau_pf)
                        plt.axvline(x=0, color='#000000', linestyle='--', linewidth=1)
                    
                    elif self.Zeus.tag_wedge == "pf":
                        plt.sca(axes[k])
                        f.text(self.xaxis_x, self.xaxis_y, self.axis_kpl_tau, fontsize=self.fontsize, ha='center')
                        f.text(self.yaxis_x, self.yaxis_y, self.axis_kpr_b, fontsize=self.fontsize, va='center', rotation='vertical')
                        plt.xticks(self.label_kpl_tau_pf)
                        plt.yticks(self.kpr_midindices, [np.round(bl, 3) for bl in self.baselines])
                        plt.xlim((-5000 * self.kpl * self.delay_range_multipler, 5000 * self.kpl * self.delay_range_multipler))
                        plt.axvline(x=0, color='#000000', linestyle='--', linewidth=1)
                        extent = [self.delays[0], self.delays[-1], self.kpr_indices[-1], 0]

                        plot = plt.imshow(
                            self.wedgeslices,
                            aspect='auto',
                            interpolation='nearest',
                            extent=extent,
                            vmax=self.vmax,
                            vmin=self.vmin)

                    elif self.Zeus.tag_wedge == "w":
                        half = self.wedgeslices.shape[1] // 2
                        self.wedgeslices = ((self.wedgeslices[:, None:half-1:-1] + self.wedgeslices[:, None:half]) / 2.).T
                        plt.sca(axes[k])
                        f.text(self.xaxis_x, self.xaxis_y, self.axis_kpr_b, fontsize=self.fontsize, ha='center')
                        f.text(self.yaxis_x, self.yaxis_y, self.axis_kpl_tau, fontsize=self.fontsize, va='center', rotation='vertical')
                        plt.yticks(self.label_kpl_tau_w)
                        plt.xticks(self.kpr_midindices, [np.round(bl, 3) for bl in self.baselines], rotation=45)
                        plt.ylim((0, 5000 * self.kpl * self.delay_range_multipler))
                        extent = [0, self.kpr_indices[-1], 0, self.delays[-1]]

                        plot = plt.imshow(
                            self.wedgeslices,
                            aspect='auto',
                            interpolation='nearest',
                            extent=extent,
                            vmax=self.vmax,
                            vmin=self.vmin)
                    
                elif self.Zeus.tag_depth == "sl":
                    self.wedgeslices = np.log10(np.fft.fftshift(np.abs(np.nanmean(data['sl'], axis=1)), axes=1))
                    if self.Zeus.tag_wedge == "1d": pass
                    elif self.Zeus.tag_wedge == "pf": pass
                    elif self.Zeus.tag_wedge == "w": pass
                
                elif self.Zeus.tag_depth == "ant":
                    self.wedgeslices = np.log10(np.fft.fftshift(np.abs(np.nanmean(data['ant'], axis=1)), axes=1))
                    if self.Zeus.tag_wedge == "1d": pass
                    elif self.Zeus.tag_wedge == "pf": pass
                    elif self.Zeus.tag_wedge == "w": pass
            
            elif self.Zeus.tag == "blavg":
                if self.Zeus.tag_depth == "bl":
                    self.wedgeslices = np.log10(np.fft.fftshift(np.abs(data['bl']), axes=1))
                    if self.Zeus.tag_wedge == "1d": pass
                    elif self.Zeus.tag_wedge == "pf": pass
                    elif self.Zeus.tag_wedge == "w": pass
                
                elif self.Zeus.tag_depth == "sl":
                    self.wedgeslices = np.log10(np.fft.fftshift(np.abs(data['sl']), axes=1))
                    if self.Zeus.tag_wedge == "1d": pass
                    elif self.Zeus.tag_wedge == "pf": pass
                    elif self.Zeus.tag_wedge == "w": pass
                
                elif self.Zeus.tag_depth == "ant":
                    self.wedgeslices = np.log10(np.fft.fftshift(np.abs(data['sl']), axes=1))
                    if self.Zeus.tag_wedge == "1d": pass
                    elif self.Zeus.tag_wedge == "pf": pass
                    elif self.Zeus.tag_wedge == "w": pass

        plt.tick_params(axis='both', direction='inout')
        
        if "1d" in self.Zeus.tag_wedge:
            plt.sca(axes[0])
            plt.legend(loc="upper left")
        else:
            cbar_ax = f.add_axes(self.p_axes)
            cbar = f.colorbar(plot, cax=cbar_ax)
            cbar.set_label(self.axis_power, fontsize=self.fontsize, ha='center')

        self.name_plot()
        plt.savefig(self.plot_name)

    def kpr_index(self):
        # Calculate kpr indices, and stretch wedgeslices array by a factor of 10
        self.kpr_indices = [int(bl * 10) for bl in self.baselines]
        wedgeslices_stretch = np.zeros((self.kpr_indices[-1], self.wedgeslices.shape[-1]), dtype=np.float64)
        j = 0
        for i in range(len(self.kpr_indices)):
            wedgeslices_stretch[j:self.kpr_indices[i]] = self.wedgeslices[i]
            j = self.kpr_indices[i]

        # Find kpr mid-indicies for the tickmarks
        self.kpr_midindices = []
        for i in range(len(self.kpr_indices)):
            if i == 0:
                self.kpr_midindices.append(self.kpr_indices[i] / 2.)
            else:
                self.kpr_midindices.append((self.kpr_indices[i] + self.kpr_indices[i - 1]) / 2.)

        self.wedgeslices = wedgeslices_stretch[...]

    def unit(self):
        # Apply astropy units to necessary quantities
        units_delays = self.delays*u.ns
        units_baselines = self.baselines*u.m

        # Get horizons before cosmo units are applied
        self.horizons = []
        for bl in units_baselines:
            self.horizons.append((bl / c.c).to(u.ns))

        # Convert frequency range from channels to Gigahertz
        units_freqrange = (self.freqrange * .1/1024. + .1)*u.GHz

        # Calculate bandwidth and center frequency
        self.bandwidth = units_freqrange[1] - units_freqrange[0]
        self.center_frequency = self.bandwidth / 2. + units_freqrange[0]

        if self.Zeus.tag_unit == "cosmo":
            """Calculate cosmological units"""
            # Calculate the redshift, z, of the 21cm line given the center frequency
            F21 = 1.42040575177*u.GHz
            z = (F21 / self.center_frequency) - 1.

            # Calculate k_parallel (kpl) and k_perpendicular (kpr) unit conversion factors
            h = cos.Planck15.h
            self.kpl = 2. * np.pi * F21 * cos.Planck15.H(z) * c.c**-1. * (1. + z)**-2.
            self.kpl = self.kpl.to(u.GHz/u.Mpc) / h
            self.kpr = 2. * np.pi * self.center_frequency * c.c**-1. * cos.Planck13.comoving_transverse_distance(z)**-1.
            self.kpr = self.kpr.to(1/(u.m*u.Mpc)) / h

            # Calculate cosmological power conversion 
            X2Y = cos.Planck15.comoving_transverse_distance(z)**2 * cos.Planck15.comoving_distance(z)
            HERA_BEAM_POLY = np.array([
                8.07774113e+08,
                -1.02194430e+09,
                5.59397878e+08,
                -1.72970713e+08,
                3.30317669e+07,
                -3.98798031e+06,
                2.97189690e+05,
                -1.24980700e+04,
                2.27220000e+02])
            BEAM = np.poly1d(HERA_BEAM_POLY)
            Omega = BEAM(self.center_frequency)*u.sr
            self.cosmo = X2Y * Omega**-1 * self.bandwidth**-1 * c.c**4. * 4.**-1. * c.k_B**-2. * self.center_frequency**-4.
            self.cosmo = self.cosmo.to(u.mK**2 * u.Mpc**3 * u.GHz**-1 * u.Jy**-2 * u.sr**-1) * h**3
            self.cosmo = np.log10(self.cosmo.value)

            # Apply unit conversions
            self.delays = ((units_delays * self.kpl).to(u.Mpc**-1)).value
            self.baselines = ((units_baselines * self.kpr).to(u.Mpc**-1)).value
            self.wedgeslices += self.cosmo

            # Declare axis labels
            self.axis_kpl_tau = r'$k_{\parallel}\ [\rm h\ Mpc^{-1}]$'
            self.axis_kpr_b = r'$k_{\perp}\ [\rm h\ Mpc^{-1}]$'
            self.axis_power = r'$\log_{10}({\rm mK^2\ Mpc^3\ h^{-3}})$'
            self.axis_kpr_b_1d = [str(bl) + r' $\rm h\ Mpc^{-1}$' for bl in np.round(self.baselines, 3)]

            # Set kpl/delay axis labels
            self.label_kpl_tau_pf = np.array([-2, -1, 0, 1, 2]) * self.delay_range_multipler
            self.label_kpl_tau_w = np.array([0, 0.5, 1, 1.5, 2, 2.5]) * self.delay_range_multipler

            # Set dynamic range
            self.vmin, self.vmax = 6, 15

            # Remove units from kpl and kpr
            self.kpl, self.kpr = self.kpl.value, self.kpr.value

        elif self.Zeus.tag_unit == "std":
            # Declare axis labels
            self.axis_kpl_tau = r'$\rm \tau\ [ns]$'
            self.axis_kpr_b = r'$\rm {\bf |b|}\ [m]$'
            self.axis_power = r'$\log_{10}({\rm mK}^2)$'
            self.axis_kpr_b_1d = [str(bl) + ' m' for bl in np.round(self.baselines, 1)]

            # Set kpl/delay axis labels
            self.label_kpl_tau_pf = np.array([-4000, -2000, 0, 2000, 4000]) * self.delay_range_multipler
            self.label_kpl_tau_w = np.array([0, 1000, 2000, 3000, 4000]) * self.delay_range_multipler

            # Set dynamic range
            self.vmin, self.vmax = -8, 0

            # Set defualt kpl, kpr values
            self.kpl, self.kpr = 1., 1.

    def horizon(self):
        # Calculate horizon lines
        j = 0
        for i in range(len(self.horizons)):
            x1, y1 = [self.horizons[i] * self.kpl, self.horizons[i] * self.kpl], [j, self.kpr_indices[i]]
            x2, y2 = [-self.horizons[i] * self.kpl, -self.horizons[i] * self.kpl], [j, self.kpr_indices[i]]
            if 'pf' in self.Zeus.tag_wedge:
                self.ax.plot(x1, y1, x2, y2, color='#ffffff', linestyle='--', linewidth=1)
            elif 'w' in self.Zeus.tag_wedge:
                self.ax.plot(y1, x1, y2, x2, color='#ffffff', linestyle='--', linewidth=1)
            j = self.kpr_indices[i]

    def name_plot(self):
        zen = 'zen'
        HH = 'HH'

        extension = self.fileZeus.ext
        pol_naming = self.fileZeus.catalog[extension].keys()[0]

        data_file0 = self.fileZeus.catalog[extension][pol_naming][0]
        data_filef = self.fileZeus.catalog[extension][pol_naming][-1]

        file0 = data_file0['file']
        filef = data_filef['file']

        JD = int(data_file0['JD'])
        JDT0 = data_file0['JD'] - JD
        JDTf = data_filef['JD'] - JD

        JD = str(JD)
        JDT0 = str(JDT0)[2:7]
        JDTf = str(JDTf)[2:7]
        
        LST0 = self.fileZeus.catalog[extension][pol_naming]['LST'][0]
        hours = int(LST0)
        mins = str((LST0 * 60) % 60).split('.')[0]
        if len(mins) == 1:
            mins = '0' + mins
        LST0 = '{HH}.{MM:.2}'.format(HH=hours, MM=mins)

        LSTf = self.fileZeus.catalog[extension][pol_naming]['LST'][-1]
        hours = int(LSTf)
        mins = str((LSTf * 60) % 60).split('.')[0]
        if len(mins) == 1:
            mins = '0' + mins
        LSTf = '{HH}.{MM:.2}'.format(HH=hours, MM=mins)

        if self.Zeus.title:
            plt.suptitle("Start: {}\nStop: {}\nBandwidth: {}\nCentral Frequency: {}".format(
                    LST0,
                    LSTf,
                    np.round(self.bandwidth.to(u.MHz), 2),
                    np.round(self.center_frequency.to(u.MHz), 2)))

        num_files = np.unique(self.fileZeus.catalog[extension][pol_naming]['file']).shape[0]
        # JDT = '{JDT0}_{num_files}_{JDTf}'.format(JDT0=JDT0, num_files=num_files, JDTf=JDTf)
        LST = '{LST0}_{num_files}_{LSTf}'.format(LST0=LST0, num_files=num_files, LSTf=LSTf)

        pol = ''.join(self.pols)

        if self.fileZeus.exants:
            exants = [str(ant) for ant in self.fileZeus.exants]
            exants = '_'.join(exants)
        else:
            exants = 'None'

        freqrange = '{start}_{end}'.format(start=self.fileZeus.freqrange[0], end=self.fileZeus.freqrange[1])
        tag = '_'.join([self.Zeus.tag_avg] + [self.Zeus.tag_depth] + [self.Zeus.tag_wedge] + [self.Zeus.tag_unit])

        # npz_name = [zen, JD, JDT, pol, exants, freqrange, HH, extension, tag, 'png']
        npz_name = [zen, JD, LST, pol, exants, freqrange, HH, extension, tag, 'png']
        npz_name = '.'.join(npz_name)

        self.plot_name = os.path.join(self.Zeus.path, npz_name)
        print(self.plot_name)


