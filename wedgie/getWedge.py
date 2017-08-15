"""


Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""

import argparse
import wedge_utils as wu
from IPython import embed
import multiprocessing

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

class Batch:
    def __init__(self, args):
        """
        Runs preliminary formatting of ex_ants list, polarization, and filenames
        based on given arguments.
        """
        self.args = args
        self.history = vars(args)
        
        self.files = None
        self.ex_ants = []
        self.pols = []
        self.pol_type = []

        self.arg_pols = [pol for pol in self.args.pol.split(',')]
        self.calfile = args.calfile.split('.')[0]
        self.freq_range = (int(args.freq_range.split('_')[0]), int(args.freq_range.split('_')[1]))


        # Generate ex_ants list from args.ex_ants.
        if args.ex_ants:
            self.ex_ants = map(int, args.ex_ants.split(','))

        # Format the polarizations to be used from args.pol.
        stokes_pols = ['I', 'Q', 'U', 'V']
        standard_pols = ['xx','xy','yx','yy']

        if any(pol in self.arg_pols for pol in stokes_pols):
            self.pol_type = 'stokes'

            if ('I' in self.arg_pols) or ('Q' in self.arg_pols):
                self.pols.extend(['xx', 'yy'])
            if ('U' in self.arg_pols) or ('V' in self.arg_pols):
                self.pols.extend(['xy', 'yx'])

        elif any(pol in self.arg_pols for pol in standard_pols):
            self.pol_type = 'standard'
            self.pols = self.arg_pols

        self.pols.sort()

        # Generates correct file names depending on polarization chosen and files given.
        self.files = {}
        for pol in self.pols:

            pol_files = []
            for file in self.args.filenames:
                file_pol = file.split('.')[-3]
                new_file = file.split(file_pol)[0] + pol + file.split(file_pol)[1]

                if not new_file in pol_files:
                    pol_files.append(new_file)

            self.files[pol] = pol_files

    def __repr__(self):
        """
        Returns self.history, the arguments passed when this instance of the Batch class 
        was initialized.
        """
        return str(self.history)

    def logic(self):
        """
        Based on the arguments passed when this instance of the Batch class was initialized 
        this function decides which functions from wedge_utils.py to execute.
        """
        if self.args.delay_avg:
            for file in self.files[0]:
                 wu.wedge_delayavg(file)

        if self.pol_type == 'stokes':
            for i in self.arg_pols:
                if args.timeavg:
                    wu.wedge_timeavg(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                elif args.blavg:
                    wu.wedge_blavg(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                elif args.flavors:
                    wu.wedge_flavors(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                elif args.bl_num:
                    wu.wedge_bltype(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                else:
                    raise Exception('You must choose a wedge type (timeavg, blavg, flavors, bl_num)')

        elif self.pol_type == 'standard':
            for i in self.arg_pols:
                if args.timeavg:
                    wu.wedge_timeavg(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                elif args.blavg:
                    wu.wedge_blavg(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                elif args.flavors:
                    wu.wedge_flavors(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                elif args.bl_num:
                    wu.wedge_bltype(self.args, self.files, i, self.calfile, self.history, self.freq_range, self.ex_ants)
                else:
                    raise Exception('You must choose a wedge type (timeavg, blavg, flavors, bl_num)')

if args.step:
    step, files = args.step, args.filenames
    del args.step, args.filenames

    files_xx = [file for file in files if 'xx' in file]
    num_files_xx = len(files_xx)

    count = 1
    for index in range(0, num_files_xx, step):
        args.filenames = files_xx[index : index + step]
        zen = Batch(args)

        if count % args.load:
            multiprocessing.Process(target=zen.logic).start()
        else:
            zen.logic()

        count += 1

    print 'Step program complete.'

elif args.stair:
    files = args.filenames
    del args.filenames

    files_xx = [file for file in files if 'xx' in file]
    num_files_xx = len(files_xx)

    count = 1
    for index in range(1, num_files_xx):
        args.filenames = files_xx[0: index]
        zen = Batch(args)
        
        if count % args.load:
            multiprocessing.Process(target=zen.logic).start()
        else:
            zen.logic()

        count += 1
    print 'Stair program complete.'

else:
    zen = Batch(args)
    zen.logic()