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
import sys

from pyuvdata import UVData
import pyuvdata.utils as uvutils
import aipy
import hera_cal
import ephem


# Interactive Development
from IPython import embed

class Eris(object):
    def __init__(self, Zeus, pol):
        self.Zeus = Zeus
        self.pol = pol

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
        self.lst = np.array([], dtype=ephem.Angle)
        self.times = list()
        self.delays = list()

    def calculate_caldata(self):
        """Returns a dictionary of baseline lengths and the corresponding pairs.
        The data is based on a calfile. ex_ants is a list of integers that
        specify antennae to be exlcuded from calculation.

        Requires cal file to be in PYTHONPATH."""
        if self.Zeus.calfile:
            try:
                print('Reading calfile: %s ...' %self.Zeus.calfile)
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
            self.lst = np.append(self.lst, lst)

            zenith = aipy.phs.RadioFixedBody(lst, self.aa.lat)
            zenith.compute(self.aa)

            if i % 2:
                v1 = fft_2Ddata[i - 1, :]
                phase_correction = np.conj(self.aa.gen_phs(zenith, antpair[0], antpair[1])) * self.aa.gen_phs(old_zenith, antpair[0], antpair[1])
                v2 = fft_2Ddata[i, :] * phase_correction
                self.vis_sq_antpair[i // 2, :] = np.conj(v1) * v2

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
            self.aa = hera_cal.utils.get_aa_from_uv(uv, self.info['freqs'])
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

    def save(self):

        self.times = np.vstack((self.lst, self.info['times']))
        channel_width = (self.info['freqs'][1] - self.info['freqs'][0]) * (10**3)
        num_bins = len(self.info['freqs'])
        self.delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))

        np.savez(self.npz_name,
                 Wedge=self)


    def name_npz(self):
        embed()
        file_start = self.files[self.files.keys()[0]][0].split('/')[-1].split('.')
        file_end = self.files[self.files.keys()[0]][-1].split('/')[-1].split('.')

        day = [str(self.info['times'][0]).split('.')[0]]
        JDT1 = [str(self.info['times'][0] - np.floor(self.info['times'][0]))[2:7]]
        JDT2 = [str(self.info['times'][-1] - np.floor(self.info['times'][-1]))[2:7]]
        JDT = ['_{files}_'.format(files=len(self.files[self.files.keys()[0]])).join(JDT1 + JDT2)]
        pol = [self.pol]
        tag = [tag]
        freq_range = ['{start}_{end}'.format(start=self.freq_range[0], end=self.freq_range[1])]

        if self.ex_ants:
            ex_ants = [str(ant) for ant in self.ex_ants]
            ex_ants = ["_".join(ex_ants)]
        else:
            ex_ants = ["None"]

        if self.args.sim:
            zen = ['zen']
            HH = ['HH']
            ext = ['SIM']
        else:
            zen = [file_start[0]]
            HH = [file_start[4]]
            ext = [file_start[-1]]

        npz_name = zen + day + JDT + pol + ex_ants + freq_range + HH + ext + tag
        npz_name = '.'.join(npz_name)
        self.npz_name = self.args.path + npz_name

        print(self.npz_name)

