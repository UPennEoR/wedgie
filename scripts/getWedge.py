"""


Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""

import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.constants as sc

import wedgie.wedge_utils as wu

import argparse, os, multiprocessing
from copy import deepcopy

# For Interactive Development
from IPython import embed

parser = argparse.ArgumentParser()

# Arguments that are always necessary:
parser.add_argument('-F',
                    '--filenames',
                    help='Input a list of filenames to be analyzed.',
                    nargs='*',
                    required=True)
parser.add_argument('-C',
                    '--calfile',
                    help='Input the calfile to be used for analysis.',
                    default=None)
parser.add_argument('-P',
                    '--pol',
                    help='Input a comma-delimited list of polatizations to plot.',
                    default='I,Q,U,V')

# Crucial arguments that are frequently needed:
parser.add_argument('-X',
                    '--ex_ants',
                    help='Input a comma-delimited list of antennae to exclude from analysis.')
parser.add_argument('-R',
                    '--freq_range',
                    help='Input a range of frequency channels to use separated by an underscore: "550_650"',
                    default='0_1023')
parser.add_argument('-p',
                    '--path',
                    help='Input the path to where you want your files to be saved.',
                    default='./')

# Types of pitchforks (only one can be used):
parser.add_argument('-t',
                    '--timeavg',
                    help='Toggle time averaging and baseling averging.',
                    action='store_true',
                    default=False)
parser.add_argument('-b',
                    '--blavg',
                    help='Toggle baseline averaging only.',
                    action='store_true',
                    default=False)
parser.add_argument('-f',
                    '--flavors',
                    help='Toggle splitting wedgeslices into a per slope per baseline basis.',
                    action='store_true',
                    default=False)
parser.add_argument('-l',
                    '--bl_type',
                    help='Toggle bltype and input 1 baseline type.',
                    default=False,
                    type=int)

# npz file manipulations:
parser.add_argument('-c',
                    '--combine',
                    help='Combines npz files to represent a longer amount of time in one npz file. Incurs some amount of error.',
                    action='store_true')
parser.add_argument('-d',
                    '--diff',
                    help='Subtract two visibilities from the nearest-neighbor lst.',
                    action='store_true')
parser.add_argument('-D',
                    '--delay_avg',
                    help='Form a wedge from a pitchfork.',
                    default=False)

# Simultaneous analysis arguments:
parser.add_argument('-S',
                    '--step',
                    help='Toggle file stepping.',
                    default=False,
                    type=int)
parser.add_argument('-A',
                    '--stair',
                    help='Compute npz files for 1 file, then 2 files, then 3 files, ...',
                    action='store_true',
                    default=False)
parser.add_argument('-L',
                    '--load',
                    help='How many processes to run at once.',
                    default=1,
                    type=int)

# Specify simulated data:
parser.add_argument('-s',
                    '--sim',
                    help='Specify sim data so that proper naming scheme is upheld.',
                    action='store_true')

args = parser.parse_args()


class Batch(object):
    def __init__(self, args):
        self.args = args
        self.history = vars(args)
        self.calfile = str()
        self.pols = list()
        self.pol_type = str()
        self.file_pols = list()
        self.files = dict()
        self.ex_ants = list()
        self.freq_range = tuple()
        self.step = int()
        self.stair = int()
        self.load = int()

        self.MISSING_TAG_ERR = "You must specify which type of pitchfork to create (--timeavg, --blavg, --flavors, --bl_type=X)."
        self.stokes_pols = ['I', 'Q', 'U', 'V']
        self.standard_pols = ['xx', 'xy', 'yx', 'yy']

        if self.args.path[-1] != '/':
            self.args.path += '/'

    def format_batch(self):
        self.format_pols()
        self.format_files()
        self.format_exants()
        self.format_freqrange()
        self.format_calfile()
        self.create_history()

    def format_pols(self):
        """Format the polarizations, e.g.: translating from IQUV to xx,xy,yx,yy"""
        self.pols = self.args.pol.split(',')

        if any(pol in self.pols for pol in self.stokes_pols):
            self.pol_type = 'stokes'

            if ('I' in self.pols) or ('Q' in self.pols):
                self.file_pols.extend(['xx', 'yy'])
            if ('U' in self.pols) or ('V' in self.pols):
                self.file_pols.extend(['xy', 'yx'])

        elif any(pol in self.pols for pol in self.standard_pols):
            self.pol_type = 'standard'
            self.file_pols = self.pols[:]

    def format_files(self):
        """Generate the filenames from given files and given polarizations"""
        self.files = {pol: [] for pol in self.file_pols}
        for file in self.args.filenames:
            for pol in self.standard_pols:
                if pol in file:
                    self.files[pol].append(file)
                    break

    def format_exants(self):
        """Takes the input ex_ants and creates a list of integers"""
        if self.args.ex_ants:
            self.ex_ants = [int(ant) for ant in self.args.ex_ants.split(',')]

    def format_freqrange(self):
        """Takes the string of frequency ranges as '550_650' and turns it into a tuple of integers as (550, 650)"""
        self.freq_range = (int(self.args.freq_range.split('_')[0]), int(self.args.freq_range.split('_')[1]))

    def format_calfile(self):
        """Remove '.py' from the end of the calfile"""
        self.calfile = self.args.calfile
        if self.calfile is not None:
            self.calfile = self.args.calfile.split('/')[-1].split('.py')[0]

    def create_history(self):
        self.history['filenames'] = self.files
        self.history['calfile'] = self.calfile
        self.history['pol'] = self.file_pols
        self.history['ex_ants'] = self.ex_ants
        self.history['freq_range'] = self.freq_range

    def stepping(self):
        self.step = self.args.step
        self.load = self.args.load

        num_files_unique = len(self.files[self.files.keys()[0]])
        files_copy = deepcopy(self.files)

        for count, index in enumerate(range(0, num_files_unique, self.step)):
            for pol in self.file_pols:
                self.files[pol] = self.files[pol][index:index + self.step]

            self.logic()
            # if (count+1) % self.load:
            #     multiprocessing.Process(target=self.logic).start()
            # else:
            #     self.logic()

            self.files = deepcopy(files_copy)

    def stairing(self):
        self.stair = self.args.stair
        self.load = self.args.load

        num_files_unique = len(self.files[self.files.keys()[0]])
        files_copy = deepcopy(self.files)

        for count, index in enumerate(range(0, num_files_unique, self.stair)):
            for pol in self.file_pols:
                self.files[pol] = self.files[pol][0:index]

            self.logic()
            # if (count+1) % self.load:
            #     multiprocessing.Process(target=self.logic).start()
            # else:
            #     self.logic()

            self.files = deepcopy(files_copy)

    def logic(self):
        if self.pol_type == 'stokes':
            for pol in self.pols:
                wedge = wu.Wedge(self.args, self.files, self.calfile, pol, self.ex_ants, self.freq_range, self.history)
                exec('wedge.form_stokes{}()'.format(pol))
                wedge.apply_flags()

                if self.args.timeavg:
                    wedge.name_npz('timeavg')
                    wedge.starter()
                    wedge.timeavg()

                elif self.args.blavg:
                    wedge.name_npz('blavg')
                    wedge.starter()
                    wedge.blavg()

                elif self.args.flavors:
                    wedge.name_npz('flavors')
                    wedge.flavors()

                elif self.args.bl_type:
                    wedge.name_npz('bl{}'.format(self.args.bl_type))
                    wedge.bltype()

                else:
                    raise Exception(self.MISSING_TAG_ERR)

                wedge.form_times()
                wedge.form_delays()
                wedge.savenpz()

        elif self.pol_type == 'standard':
            for pol in self.pols:
                wedge = wu.Wedge(self.args, self.files, self.calfile, pol, self.ex_ants, self.freq_range, self.history)
                wedge.load_file()
                wedge.apply_flags()

                if self.args.timeavg:
                    wedge.name_npz('timeavg')
                    wedge.starter()
                    wedge.timeavg()

                elif self.args.blavg:
                    wedge.name_npz('blavg')
                    wedge.starter()
                    wedge.blavg()

                elif self.args.flavors:
                    wedge.name_npz('flavors')
                    wedge.flavors()

                elif self.args.bl_type:
                    wedge.name_npz('bl{}'.format(self.args.bl_type))
                    wedge.bltype()

                else:
                    raise Exception(self.MISSING_TAG_ERR)

                wedge.form_times()
                wedge.form_delays()
                wedge.savenpz()
        else:
            raise Exception("Polarization type not understood, be sure you have correctly specified the polarizations you want.")

    def combine(self):
        self.format_pols()

        # This sets up the format of self.files to be {pol: [file1, file2, ...]}
        for pol in self.pols:
            pol_files = []
            for file in self.args.filenames:
                if '.{}.'.format(pol) in file.split('/')[-1]:
                    pol_files.append(file)

            self.files[pol] = pol_files

        for pol in self.pols:

            file_0 = self.files[pol][0]
            with np.load(file_0)  as data_0:
                pol = data_0['pol'].tolist()
                caldata = data_0['caldata'].tolist()
                cwslices = data_0['cwslices']
                delays = data_0['delays']
                times = data_0['times']

            # Cycle through the rest of the files in the polarization
            for npz in self.files[pol][1:]:
                with np.load(npz) as data:
                    cwslices = np.concatenate((cwslices, data['cwslices']), axis=1)
                    times = np.concatenate((times, data['times']), axis=1)

            # Take the log of the fftshift of the abs of the time average of the cwedgeslices (to make wedgeslices)
            wslices = np.log10(np.fft.fftshift(np.abs(np.nanmean(cwslices, axis=1)), axes=1))

            # Naming and Saving
            file1 = self.files[pol][0].split('/')[-1].split('.')
            file2 = self.files[pol][-1].split('/')[-1].split('.')

            zen1 = [file1[0]]
            zen2 = [file2[0]]
            if zen1 == zen2:
                zen = zen1
            else:
                raise Exception("Your files aren't the same (zen).")

            day1 = [file1[1]]
            day2 = [file2[1]]
            day = ['__'.join(day1 + day2)]

            time1 = [file1[2]]
            time2 = [file2[2]]
            time = ['__'.join(time1 + time2)]

            ex_ants1 = [file1[4]]
            ex_ants2 = [file2[4]]
            if ex_ants1 == ex_ants2:
                ex_ants = ex_ants1
            else:
                raise Exception("Your files aren't the same (ex_ants).")

            freq_range1 = [file1[5]]
            freq_range2 = [file2[5]]
            if freq_range1 == freq_range2:
                freq_range = freq_range1
            else:
                raise Exception("Your files aren't the same (freq_range).")

            HH1 = [file1[6]]
            HH2 = [file2[6]]
            if HH1 == HH2:
                HH = HH1
            else:
                raise Exception("Your files aren't the same (HH).")

            ext1 = [file1[7]]
            ext2 = [file2[7]]
            if ext1 == ext2:
                ext = ext1
            else:
                raise Exception("Your files aren't the same (ext).")

            npz_name = '.'.join(zen + day + time + [pol] + ex_ants + freq_range + HH + ext + ['timeavg'])

            np.savez(
                npz_name,
                pol=pol,
                caldata=caldata,
                wslices=wslices,
                cwslices=cwslices,
                delays=delays,
                times=times)

            for npz in self.files[pol]:
                os.remove(npz)
                
    def difference(self):
        self.format_pols()

        # This sets up the format of self.files to be {pol: [file1, file2, ...]}
        for pol in self.pols:
            pol_files = []
            for file in self.args.filenames:
                if '.{}.'.format(pol) in file.split('/')[-1]:
                    pol_files.append(file)

            self.files[pol] = pol_files

        for pol in self.pols:

            if len(self.files[pol]) > 2:
                raise Exception('Can only take the difference of two wedges. Provide only two files.')

            file1 = self.files[pol][0]
            data1 = np.load(file1)
            cwedge1, times1 = data1['cwslices'], data1['times'][:, 1::2]

            file2 = self.files[pol][1]
            data2 = np.load(file2)
            cwedge2, times2 = data2['cwslices'], data2['times'][:, 1::2]

            if np.all(data1['delays'] == data2['delays']) and np.all(data1['caldata'] == data2['caldata']):
                delays = data1['delays']
                caldata = data1['caldata']
            else:
                print 'Warning: delays and/or caldata arrays are not the same.'
                delays = data1['delays']
                caldata = data1['caldata']

            
            # Find overlap in lst
            indices1 = np.where(np.logical_and(times1[2] >= times2[2][0], times1[2] <= times2[2][-1]))[0]
            indices2 = np.where(np.logical_and(times2[2] >= times1[2][0], times2[2] <= times1[2][-1]))[0]

            # Create times array with overlapping lst range.
            times1 = np.take(times1, indices1, axis=1)
            times2 = np.take(times2, indices2, axis=1)
            times = times1[...]

            # Make cwedge arrays represent only data from overlapping lst.
            cwedge1 = np.take(cwedge1, indices1, axis=1)
            cwedge2 = np.take(cwedge2, indices2, axis=1)

            # Create cwedgeslices array to return information before time averaging.
            cwedgeslices = cwedge1 - cwedge2

            # Average the data over time axis.
            wedge1 = np.nanmean(cwedge1, axis=1)
            wedge2 = np.nanmean(cwedge2, axis=1)

            # Create wedgeslices array to be plotted.
            wedgeslices = np.fft.fftshift(np.abs(wedge1 - wedge2), axes=1)
            wedgeslices = np.where(wedgeslices > 0, np.log10(wedgeslices), -100)

            file1 = file1.split('/')[-1].split('.')
            file2 = file2.split('/')[-1].split('.')

            zen1 = [file1[0]]
            zen2 = [file2[0]]
            if zen1 == zen2:
                zen = zen1
            else:
                raise Exception("Your files aren't the same (zen).")

            day1 = [file1[1]]
            day2 = [file2[1]]
            day = ['--'.join(day1 + day2)]

            time1 = [file1[2]]
            time2 = [file2[2]]
            time = ['--'.join(time1 + time2)]

            ex_ants1 = [file1[4]]
            ex_ants2 = [file2[4]]
            if ex_ants1 == ex_ants2:
                ex_ants = ex_ants1
            else:
                ex_ants = ex_ants1
                print("Your files aren't the same (ex_ants).")

            freq_range1 = [file1[5]]
            freq_range2 = [file2[5]]
            if freq_range1 == freq_range2:
                freq_range = freq_range1
            else:
                raise Exception("Your files aren't the same (freq_range).")

            HH1 = [file1[6]]
            HH2 = [file2[6]]
            if HH1 == HH2:
                HH = HH1
            else:
                raise Exception("Your files aren't the same (HH).")

            ext1 = [file1[7]]
            ext2 = [file2[7]]
            if ext1 == ext2:
                ext = ext1
            else:
                ext = ['--'.join(ext1 + ext2)]

            npz_name = '.'.join(zen + day + time + [pol] + ex_ants + freq_range + HH + ext + ['diff'])

            np.savez(
                npz_name,
                pol=pol,
                caldata=caldata,
                wslices=wedgeslices,
                cwslices=cwedgeslices,
                delays=delays,
                times=times)


zen = Batch(args)

if args.diff:
    zen.difference()
    quit()
elif args.combine:
    zen.combine()
    quit()

zen.format_batch()
if args.step:
    zen.stepping()
elif args.stair:
    zen.stairing()
else:
    zen.logic()
