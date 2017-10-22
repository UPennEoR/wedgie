"""


Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.constants as sc

import wedge_utils as wu

import argparse, os, multiprocessing
from copy import deepcopy

# For Interactive Development
from IPython import embed

parser = argparse.ArgumentParser()

# Main Arguemnts
parser.add_argument('-F',
                    '--filenames',
                    help='Input a list of filenames to be analyzed.',
                    nargs='*',
                    required=True)
parser.add_argument('-C',
                    '--calfile',
                    help='Input the calfile to be used for analysis.',
                    default='hsa7458_v001')
parser.add_argument('-P',
                    '--pol',
                    help='Input a comma-delimited list of polatizations to plot.',
                    default='I,Q,U,V')
parser.add_argument('-X',
                    '--ex_ants',
                    help='Input a comma-delimited list of antennae to exclude from analysis.',
                    type=str)
parser.add_argument('-V',
                    '--path',
                    help='Input the path to where you want your files to be saved.',
                    default='')

# Parameters for Changing What is Analyzed
parser.add_argument('-S',
                    '--step',
                    help='Toggle file stepping.',
                    default=False,
                    type=int)
parser.add_argument('-A',
                    '--stair',
                    help='Compute npz files for 1 file, then 2 files, then 3 files, ...',
                    default=False,
                    action='store_true')
parser.add_argument('-L',
                    '--load',
                    help='How many processes to run at once.',
                    type=int,
                    default=1)
parser.add_argument('-R',
                    '--freq_range',
                    help='Input a range of frequency channels to use separated by an underscore: "550_650"',
                    default='0_1023')
parser.add_argument('-D',
                    '--Difference',
                    action='store_true')

# Types of Wedges
# Only One Can Be Used
parser.add_argument('-t',
                    '--timeavg',
                    help='Toggle time averaging.',
                    default=False,
                    action='store_true')
parser.add_argument('-b',
                    '--blavg',
                    help='Toggle blavg for stokes.',
                    default=False,
                    action='store_true')
parser.add_argument('-f',
                    '--flavors',
                    help='Toggle splitting wedgeslices into a per slope per baseline basis.',
                    default=False,
                    action='store_true')
parser.add_argument('-l',
                    '--bl_type',
                    help='Toggle bltype and input 1 baseline type.',
                    default=False,
                    type=int)

# Delay Average (Pitchfork --> Wedge)
parser.add_argument('-d',
                    '--delay_avg',
                    help="sfsdfasdfsf",
                    default=False,
                    action="store_true")

# Combine npz files
parser.add_argument('-c',
                    '--combine',
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

        self.MISSING_TAG_ERR = "You must specify which type of Wedge to create (--timeavg, --blavg, --flavors, --bl_type=X)."

    def format_batch(self):
        self.format_pols()
        self.format_files()
        self.format_exants()
        self.format_freqrange()
        self.format_calfile()
        self.create_history()

    def format_pols(self):
        """Format the polarizations, e.g.: translating from IQUV to xx,xy,yx,yy"""
        stokes_pols = ['I', 'Q', 'U', 'V']
        standard_pols = ['xx', 'xy', 'yx', 'yy']

        self.pols = self.args.pol.split(',')

        if any(pol in self.pols for pol in stokes_pols):
            self.pol_type = 'stokes'

            if ('I' in self.pols) or ('Q' in self.pols):
                self.file_pols.extend(['xx', 'yy'])
            if ('U' in self.pols) or ('V' in self.pols):
                self.file_pols.extend(['xy', 'yx'])

        elif any(pol in self.pols for pol in standard_pols):
            self.pol_type = 'standard'
            self.file_pols = self.pols[:]

    def format_files(self):
        """Generate the filenames from given files and given polarizations"""
        for pol in self.file_pols:
            pol_files = []
            for file in self.args.filenames:
                file_pol = file.split('.')[-3]
                pol_file = file.split(file_pol)[0] + pol + file.split(file_pol)[1]
                pol_files.append(pol_file)

            self.files[pol] = sorted(list(set(pol_files)))

    def format_exants(self):
        """Takes the input ex_ants and creates a list of integers"""
        if self.args.ex_ants:
            self.ex_ants = [int(ant) for ant in self.args.ex_ants.split(',')]

    def format_freqrange(self):
        """Takes the string of frequency ranges as '550_650' and turns it into a tuple of integers as (550, 650)"""
        self.freq_range = (int(self.args.freq_range.split('_')[0]), int(self.args.freq_range.split('_')[1]))

    def format_calfile(self):
        """Remove '.py' from the end of the calfile"""
        self.calfile = self.args.calfile.split('.py')[0]

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
                wedge.format_flags()

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
                wedge.format_flags()

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
        self.pols = self.args.pol.split(',')

        # This sets up the format of self.files to be {pol: [file1, file2, ...]}
        for pol in self.pols:
            pol_files = []
            for file in self.args.filenames:
                if '.{}.'.format(pol) in file.split('/')[-1]:
                    pol_files.append(file)

            self.files[pol] = pol_files

        for pol in self.pols:

            # Load data from the first file to start out with
            # Grab the delays and caldata which will be the same for every file
            # Grab the times which will eventually be combined with every other file's times
            file_0 = self.files[pol][0]
            data_0 = np.load(file_0)
            caldata = data_0['caldata']
            delays = data_0['delays']
            cwedgeslices = data_0['cwslices']
            times = data_0['times']

            # Cycle through the rest of the files in the polarization
            # Add to rolling sum of cwedgeslice data
            # Combine times
            for i, npz in enumerate(self.files[pol][1:]):
                data = np.load(npz)
                cwedgeslices = np.concatenate((cwedgeslices, data['cwslices']), axis=1)
                times = np.concatenate((times, data['times']), axis=1)

            # Average together cwedgeslices
            # cwedgeslices /= len(self.files[pol])

            # Take the log of the fftshift of the abs of the time average of the cwedgeslices (to make wedgeslices)
            wedgeslices = np.log10(np.fft.fftshift(np.abs(np.mean(cwedgeslices, axis=1)), axes=1))

            # Naming and Saving
            start = file_0.split('/')[-1].split('.')[2].split('_')[0]
            end = file.split('/')[-1].split('.')[2].split('_')[1]
            time_rng = ['_'.join([start, end])]
            npz_name = '.'.join(file_0.split('/')[-1].split('.')[:2] + time_rng + file.split('/')[-1].split('.')[3:])
            np.savez(npz_name, cwslices=cwedgeslices, wslices=wedgeslices, delays=delays, caldata=caldata, pol=pol, times=times)

            # Remove npz files that have been combined
            for npz in self.files[pol]:
                os.remove(npz)

    def difference(self):
        if len(self.args.filenames) > 2:
            raise Exception('Can only take the difference of two wedges. Provide only two files.')

        file1 = self.args.filenames[0]
        data1 = np.load(file1)
        cwedge1, times1 = data1['cwslices'], data1['times'][:, 1::2]

        file2 = self.args.filenames[1]
        data2 = np.load(file2)
        cwedge2, times2 = data2['cwslices'], data2['times'][:, 1::2]

        if np.all(data1['delays'] == data2['delays']) and np.all(data1['caldata'] == data2['caldata']):
            delays = data1['delays']
            caldata = data1['caldata']
        else:
            raise Exception('delays and caldata arrays are not the same.')

        for i, time in enumerate(times1[0]):
            if time >= times2[0, 0]:
                start = i
                break

        for i, time in enumerate(times2[0]):
            if time >= times1[0, -1]:
                end = i
                break

        times1 = times1[:, start:]
        times2 = times2[:, :end]

        cwedge1 = np.mean(cwedge1[:, start:, :], axis=1)
        cwedge2 = np.mean(cwedge2[:, :end, :], axis=1)

        wedgeslices = np.fft.fftshift(np.abs(cwedge1 - cwedge2), axes=0)

        # Plot difference data.
        plotindeces = [int(round(i*10)) for i in caldata[3]]
        plotdata = np.zeros((plotindeces[-1], wedgeslices.shape[-1]), dtype=np.float64)
        j = 0
        for i in range(len(plotindeces)):
            plotdata[j:plotindeces[i]] = wedgeslices[i]
            j = plotindeces[i]

        plt.imshow(
            plotdata,
            aspect='auto',
            interpolation='nearest',
            extent=[delays[0], delays[-1], plotindeces[-1], 0],
            vmin=0,
            vmax=0.025)

        plt.colorbar().set_label('Linear Difference')
        plt.xlabel("Delay [ns]")
        plt.xlim((-450, 450))
        plt.ylabel("Baseline Length [m]")
        plt.yticks(plotindeces, [round(n, 1) for n in caldata[3]])

        plt.suptitle("JD: {JD1} - {JD2}; LST {start} to {end}".format(
            JD1=file1.split('/')[-1].split('.')[1],
            JD2=file2.split('/')[-1].split('.')[1],
            start=times1[1][0][:-6],
            end=times2[1][-1][:-6]))
        plt.title(file1.split('/')[-1].split('.')[3])

        plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

        horizons = []
        for length in caldata[3]:
            horizons.append(length / sc.c * 10**9)
        j = 0
        for i in range(len(horizons)):
            x1, y1 = [horizons[i], horizons[i]], [j, plotindeces[i]]
            x2, y2 = [-horizons[i], -horizons[i]], [j, plotindeces[i]]
            plt.plot(x1, y1, x2, y2, color='white', linestyle='--', linewidth=.75)
            j = plotindeces[i]

        file1 = file1.split('/')[-1].split('.')
        file2 = file2.split('/')[-1].split('.')

        time1 = ['_'.join(file1[1:3])]
        time2 = ['_'.join(file2[1:3])]
        time = ['-'.join(time1 + time2)]

        end1 = file1[3:-2]
        end2 = file2[3:-2]

        zen1 = [file1[0]]
        zen2 = [file2[0]]

        if end1 != end2 or zen1 != zen2:
            raise Exception('You did not supply the same type of file!')

        file = '.'.join(zen1 + time + end1)

        plt.savefig(self.args.path + file + '.diff.png')
        plt.show()
        plt.close()
        plt.clf()


zen = Batch(args)

if args.Difference:
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
