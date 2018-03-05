"""
This module contains functions to create and plot wedges/pitchforks from HERA data files.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
Paul La Plante <plaplant_at_sas.upenn.edu>
"""

import numpy as np
import scipy.constants as sc

import os

from pyuvdata import UVData
import pyuvdata.utils as uvutils
import aipy
import hera_cal
import ephem
import matplotlib.pyplot as plt


# Interactive Development
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
        self.vis_sq_antpair = None

        self.cwedgeslices = list()
        self.wedgeslices = list()
        self.delays = list()

    def name_npz(self):
        file0 = self.Zeus.catalog[self.Zeus.catalog.keys()[0]][0].keys()[0]
        filef = self.Zeus.catalog[self.Zeus.catalog.keys()[0]][-1].keys()[0]
        JDT0 = self.Zeus.catalog[self.Zeus.catalog.keys()[0]][0][file0]['JD'][0]
        JDTf = self.Zeus.catalog[self.Zeus.catalog.keys()[0]][-1][filef]['JD'][0]
        JD = int(JDT0)

        JDT0 = [str(JDT0 - JD)[2:7]]
        JDTf = [str(JDTf - JD)[2:7]]
        JD = [str(JD)]
        JDT = ['_{files}_'.format(files=len(self.Zeus.files[self.Zeus.files.keys()[0]])).join(JDT0 + JDTf)]
        pol = [self.pol]
        freqrange = ['{start}_{end}'.format(start=self.Zeus.freqrange[0], end=self.Zeus.freqrange[1])]
        tag = [self.Zeus.tag]
        zen = ['zen']
        HH = ['HH']

        if self.Zeus.exants:
            exants = [str(ant) for ant in self.Zeus.exants]
            exants = ['_'.join(exants)]
        else:
            exants = ['None']

        ext = [file0.split('.')[-1]]

        npz_name = zen + JD + JDT + pol + exants + freqrange + HH + ext + tag + ['npz']
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
                baselines = np.append(baselines, baseline - 0.1) 

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

        for baseline in sorted(self.caldata['antdict']):
            self.vis_sq_bl = np.zeros((ntimes // 2, nchan), dtype=np.complex128)
            for antpair in sorted(self.caldata['antdict'][baseline]):
                self.vis_sq_antpair = np.zeros((ntimes // 2, nchan), dtype=np.complex128)

                self.fft_clean_phase(antpair, ntimes, self.Zeus.CLEAN)
                
                self.vis_sq_bl += self.vis_sq_antpair
            self.vis_sq_bl /= len(self.caldata['antdict'][baseline])
            self.cwedgeslices.append(self.vis_sq_bl)
            self.wedgeslices = np.log10(np.fft.fftshift(np.abs(np.nanmean(self.cwedgeslices, axis=1)), axes=1))

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
        channel_width = (self.info['freqs'][1] - self.info['freqs'][0]) * (10**3)
        num_bins = len(self.info['freqs'])
        self.delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))

        np.savez(self.npz_name,
                 Zeus=self.Zeus,
                 wslices=self.wedgeslices,
                 cwslices=self.cwedgeslices,
                 caldata=self.caldata,
                 pol=self.pol,
                 freqs=self.info['freqs'],
                 delays=self.delays)

class Ares(object):
    def __init__(self, Zeus):
        self.Zeus = Zeus

    def makePlot(self):
        nplots = len(self.Zeus.inputfiles)
        self.f, self.axes = plt.subplots(1, nplots, sharex=True, sharey=True, figsize=(20, 10))
        if nplots == 1:
            self.axes = [self.axes]
        self.f.subplots_adjust(wspace=0)
        self.pols = []
        
        for k, self.ax in enumerate(self.axes):
            file0 = self.Zeus.inputfiles[k]
            data = np.load(file0)
            delays = data['delays']

            self.pols.append(data['pol'].tolist())
            self.fileZeus = data['Zeus'].tolist()
            self.caldata = data['caldata'].tolist()
            self.wedgeslices = data['wslices']

            self.plot_parameters()

            self.ax.set_title(data['pol'], fontsize=20)
            self.ax.tick_params(axis='both', direction='inout')
            if self.Zeus.tag_wedge == "pf":
                self.ax.set_xlim((-500*self.k_para_factor, 500*self.k_para_factor))
                self.ax.axvline(x=0, color='#000000', linestyle='--', linewidth=1)
                plt.sca(self.axes[k])
                plt.yticks(self.k_perp_midindices, [np.round(self.k_perp_factor * n, 3) for n in self.caldata['baselines']])
                extent = [delays[0]*self.k_para_factor, delays[-1]*self.k_para_factor, self.k_perp_indices[-1], 0]

            elif self.Zeus.tag_wedge == "w":
                self.plotdata = ((self.plotdata[:, None:49:-1] + self.plotdata[:, None:50])/2).T
                self.ax.set_ylim((0*self.k_para_factor, 500*self.k_para_factor))
                plt.sca(self.axes[k])
                plt.xticks(self.k_perp_midindices, [np.round(self.k_perp_factor * n, 3) for n in self.caldata['baselines']], rotation=45)
                extent = [0, self.k_perp_indices[-1], 0, delays[-1]*self.k_para_factor]

            if "1d" not in self.Zeus.tag_wedge:
                plot = self.ax.imshow(
                    self.plotdata,
                    aspect='auto',
                    interpolation='nearest',
                    extent=extent,
                    vmax=self.vmax,
                    vmin=self.vmin)
            else:
                for i, bl in enumerate(self.caldata['baselines']):
                    self.ax.plot(
                        delays*self.k_para_factor,
                        self.wedgeslices[i] + self.cosmo,
                        label=str(np.round(self.caldata['baselines'][i], 1)) + ' m')
                plt.sca(self.axes[0])
                plt.ylim((self.vmin, self.vmax))
                plt.xlim((-500*self.k_para_factor, 500*self.k_para_factor))

        start = self.fileZeus.catalog['xx'][0][self.fileZeus.catalog['xx'][0].keys()[0]]['LST'][0]
        stop = self.fileZeus.catalog['xx'][-1][self.fileZeus.catalog['xx'][-1].keys()[0]]['LST'][-1]
        plt.suptitle("Start: " + start +"\nStop: " + stop, fontsize=20)
        
        if "1d" not in self.Zeus.tag_wedge:
            cbar_ax = self.f.add_axes([0.9125, 0.25, 0.025, 0.5])
            cbar = self.f.colorbar(plot, cax=cbar_ax)
            cbar.set_label(self.power_axis, fontsize=20, ha='center')
        else:
            plt.sca(self.axes[0])
            plt.legend(loc="upper left")

        self.name_plot()
        plt.savefig(self.plot_name)

    def plot_parameters(self):

        # Format cosmo unit factor
        if self.Zeus.args.tag_data == "cosmo":
            if self.fileZeus.freqrange == [580, 680]:
                self.cosmo = 15.55
            elif self.fileZeus.freqrange == [200, 300]:
                self.cosmo = 16.2
        else:
            self.cosmo = 0.
        self.tag_data = self.Zeus.args.tag_data

        # Find k_perp indices
        self.k_perp_indices = [int(i*10) for i in self.caldata['baselines']]
        self.plotdata = np.zeros((self.k_perp_indices[-1], self.wedgeslices.shape[-1]), dtype=np.float64)
        j = 0
        for i in range(len(self.k_perp_indices)):
            self.plotdata[j:self.k_perp_indices[i]] = self.wedgeslices[i]
            j = self.k_perp_indices[i]

        # Find k_perp mid-indicies
        self.k_perp_midindices = []
        for i in range(len(self.k_perp_indices)):
            if i == 0:
                self.k_perp_midindices.append(self.k_perp_indices[i] / 2.)
            else:
                self.k_perp_midindices.append((self.k_perp_indices[i] + self.k_perp_indices[i-1]) / 2.)

        # Declare factors and figure out axis labels
        if self.tag_data == "cosmo":                
            self.plotdata += self.cosmo
            self.power_axis = r'$\log_{10}({\rm mK^2 Mpc^3 h^{-3}})$'
            self.vmin, self.vmax = 7, 17
            if self.fileZeus.freqrange == [580, 680]:
                self.k_para_factor, self.k_perp_factor = 5.474e-4, 6.4e-4
            elif self.fileZeus.freqrange == [200, 300]:
                self.k_para_factor, self.k_perp_factor = 6.239e-4, 4.6e-4
            if "pf" in self.Zeus.tag_wedge:
                if self.Zeus.tag_wedge == "1dpf":
                    self.f.text(0.5, 0.025, r'$k_{\parallel}\ [\rm h\ Mpc^{-1}]$', fontsize=20, ha='center')
                    self.f.text(0.075, 0.5, self.power_axis, fontsize=20, va='center', rotation='vertical')
                    self.ax.set_xticks([-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4])
                else:
                    self.f.text(0.5, 0.025, r'$k_{\parallel}\ [\rm h\ Mpc^{-1}]$', fontsize=20, ha='center')
                    self.f.text(0.06, 0.5, r'$k_{\perp}\ [\rm h\ Mpc^{-1}]$', fontsize=20, va='center', rotation='vertical')
                    self.ax.set_xticks([-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4])
            
            elif "w" in self.Zeus.tag_wedge:
                self.f.text(0.5, 0.025, r'$k_{\perp}\ [\rm h\ Mpc^{-1}]$', fontsize=20, ha='center')
                self.f.text(0.06, 0.5, r'$k_{\parallel}\ [\rm h\ Mpc^{-1}]$', fontsize=20, va='center', rotation='vertical')
                self.ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25])

        elif self.tag_data == "std":
            self.power_axis = r'$\log_{10}({\rm mK}^2)$'
            self.vmin, self.vmax = -8, 1
            self.k_para_factor, self.k_perp_factor = 1., 1.
            if "pf" in self.Zeus.tag_wedge:
                if self.Zeus.tag_wedge == "1dpf":
                    self.f.text(0.5, 0.025, 'Delay [ns]', fontsize=20, ha='center')
                    self.f.text(0.075, 0.5, self.power_axis, fontsize=20, va='center', rotation='vertical')
                    self.ax.set_xticks([-400, -200, 0, 200, 400])
                else:
                    self.f.text(0.5, 0.025, 'Delay [ns]', fontsize=20, ha='center')
                    self.f.text(0.075, 0.5, 'Baseline Length [m]', fontsize=20, va='center', rotation='vertical')
                    self.ax.set_xticks([-400, -200, 0, 200, 400])
            
            elif "w" in self.Zeus.tag_wedge:
                self.f.text(0.5, 0.025, 'Baseline Length [m]', fontsize=20, ha='center')
                self.f.text(0.075, 0.5, 'Delay [ns]', fontsize=20, va='center', rotation='vertical')
                self.ax.set_yticks([0, 100, 200, 300, 400, 500])

        # Calculate horizon lines
        horizons = []
        for length in self.caldata['baselines']:
            horizons.append(length / sc.c * 10**9)
        j = 0
        for i in range(len(horizons)):
            x1, y1 = [horizons[i]*self.k_para_factor, horizons[i]*self.k_para_factor], [j, self.k_perp_indices[i]]
            x2, y2 = [-horizons[i]*self.k_para_factor, -horizons[i]*self.k_para_factor], [j, self.k_perp_indices[i]]
            if self.Zeus.tag_wedge == "pf":
                self.ax.plot(x1, y1, x2, y2, color='#ffffff', linestyle='--', linewidth=1)
            elif self.Zeus.tag_wedge == "w":
                self.ax.plot(y1, x1, y2, x2, color='#ffffff', linestyle='--', linewidth=1)
            j = self.k_perp_indices[i]

    def name_plot(self):
        file0 = self.fileZeus.catalog[self.fileZeus.catalog.keys()[0]][0].keys()[0]
        filef = self.fileZeus.catalog[self.fileZeus.catalog.keys()[0]][-1].keys()[0]
        JDT0 = self.fileZeus.catalog[self.fileZeus.catalog.keys()[0]][0][file0]['JD'][0]
        JDTf = self.fileZeus.catalog[self.fileZeus.catalog.keys()[0]][-1][filef]['JD'][0]
        JD = int(JDT0)

        JDT0 = [str(JDT0 - JD)[2:7]]
        JDTf = [str(JDTf - JD)[2:7]]
        JD = [str(JD)]
        JDT = ['_{files}_'.format(files=len(self.fileZeus.catalog[self.fileZeus.catalog.keys()[0]])).join(JDT0 + JDTf)]
        pol = [''.join(self.pols)]
        freqrange = ['{start}_{end}'.format(start=self.fileZeus.freqrange[0], end=self.fileZeus.freqrange[1])]
        tag = ['_'.join([self.Zeus.tag_wedge] + [self.Zeus.tag] + [self.tag_data])]

        zen = ['zen']
        HH = ['HH']

        if self.fileZeus.exants:
            exants = [str(ant) for ant in self.fileZeus.exants]
            exants = ['_'.join(exants)]
        else:
            exants = ['None']

        ext = [file0.split('.')[-1]]

        self.plot_name = zen + JD + JDT + pol + exants + freqrange + HH + ext + tag + ['pdf']
        self.plot_name = '.'.join(self.plot_name)
        self.plot_name = os.path.join(self.Zeus.path, self.plot_name)
