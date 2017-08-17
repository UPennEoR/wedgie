"""


Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""
import argparse
import wedge_utils as wu
import multiprocessing
from copy import deepcopy

# For Interactive Decelopment
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
                    '--bl_num',
                    help='Toggle bltype and input 1 baseline type.',
                    default=False,
                    type=int)

# Delay Average (Pitchfork --> Wedge)
parser.add_argument('-d',
                    '--delay_avg',
                    help="sfsdfasdfsf",
                    default=False,
                    action="store_true")

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

        self.format_pols()
        self.format_files()
        self.format_exants()
        self.format_freqrange()
        self.format_calfile()
        self.create_history()

    def format_pols(self):
        """Format the polarizations, e.g.: translating from IQUV to xx,xy,yx,yy"""
        stokes_pols = ['I', 'Q', 'U', 'V']
        standard_pols = ['xx','xy','yx','yy']

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
                self.files[pol] = self.files[pol][index : index + self.step]

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
                self.files[pol] = self.files[pol][0 : index]

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
                
                if self.args.timeavg:
                    wedge.name_npz('timeavg')
                    wedge.timeavg()

                elif self.args.blavg:
                    wedge.name_npz('blavg')
                    wedge.blavg()

                elif self.args.flavors:
                    wedge.name_npz('flavors')
                    wedge.flavors()

                elif self.args.bl_num:
                    wedge.name_npz('bl{}'.format(self.args.bl_num))
                    wedge.bltype()

                else:
                    raise Exception("You must specify which type of Wedge to create (--timeavg, --blavg, --flavors, --bl_num=X).")

        elif self.pol_type == 'standard':
            for pol in self.pols:
                wedge = wu.Wedge(self.args, self.files, self.calfile, pol, self.ex_ants, self.freq_range, self.history)
                wedge.load_file()

                if self.args.timeavg:
                    wedge.name_npz('timeavg')
                    wedge.timeavg()

                elif self.args.blavg:
                    wedge.name_npz('blavg')
                    wedge.blavg()

                elif self.args.flavors:
                    wedge.name_npz('flavors')
                    wedge.flavors()

                elif self.args.bl_num:
                    wedge.name_npz('bl{}'.format(self.args.bl_num))
                    wedge.bltype()

                else:
                    raise Exception("You must specify which type of Wedge to create (--timeavg, --blavg, --flavors, --bl_num=X).")

        else:
            raise Exception("Polarization type not understood, be sure you have correctly specified the polarizations you want.")


zen = Batch(args)
if args.step:
    zen.stepping()
elif args.stair:
    zen.stairing()
else:
    zen.logic()