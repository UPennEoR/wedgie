import argparse
import wedge_utils as wu
from IPython import embed
import multiprocessing

parser = argparse.ArgumentParser()
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
                    default='stokes')
parser.add_argument('-f',
                    '--flavors',
                    help='Toggle splitting wedgeslices into a per slope per baseline basis.',
                    action='store_true')
parser.add_argument('-t',
                    '--time_avg',
                    help='Toggle off time averaging.',
                    default=True,
                    action='store_false')
parser.add_argument('-x',
                    '--ex_ants',
                    help='Input a comma-delimited list of antennae to exclude from analysis.',
                    type=str)
parser.add_argument('-s',
                    '--step',
                    help='Toggle file stepping.',
                    type=int)
parser.add_argument('-L',
                    '--load',
                    help='How many processes to run at once.',
                    type=int,
                    default=1)
parser.add_argument('-r',
                    '--freq_range',
                    help='Input a range of frequency channels to use separated by an underscore: "550_650"',
                    default='0_1023')
parser.add_argument('-a',
                    '--stair',
                    help='Compute npz files for 1 file, then 2 files, then 3 files, ...',
                    action='store_true')
parser.add_argument('-d',
                    '--delay_avg',
                    help="sfsdfasdfsf",
                    action="store_true")
parser.add_argument('-b',
                    '--blavg',
                    help='Toggle blavg for stokes.',
                    action='store_true')
parser.add_argument('-l',
                    '--bl_num',
                    help='Toggle bltype and input 1 baseline type.',
                    type=int)
args = parser.parse_args()

class Batch:
    def __init__(self, args):
        self.args = args
        self.history = vars(args)
        
        self.files = None
        self.pols = None
        self.pol_type = None
        self.calfile = args.calfile.split('.')[0]
        self.freq_range = (int(args.freq_range.split('_')[0]), int(args.freq_range.split('_')[1]))
        self.ex_ants = []

        # Generate ex_ants list from args.ex_ants.
        if args.ex_ants is not None:
            self.ex_ants = map(int, args.ex_ants.split(','))

        # Format the polarizations to be used from args.pol.
        self.pols = [pol.lower() for pol in self.args.pol.split(',')]
        num_pols = len(self.pols)

        if self.pols == ['stokes']:
            self.pols = ['xx','xy','yx','yy']
            self.pol_type = 'stokes'
        elif num_pols == 1:
            self.pol_type = 'single'
        elif num_pols > 1:
            self.pol_type = 'multi'

        # Generates correct file names depending on polarization chosen and files given.
        self.files = []
        for pol in self.pols:

            pol_files = []
            for file in self.args.filenames:
                file_pol = file.split('.')[-3]
                new_file = file.split(file_pol)[0] + pol + file.split(file_pol)[1]

                if not new_file in pol_files:
                    pol_files.append(new_file)

            self.files.append(pol_files)

    def __repr__(self):
        return str(self.history)

    def logic(self):
        if self.pol_type == 'stokes':
                wu.wedge_stokes(self.args, self.files, self.calfile, self.history, self.freq_range, self.ex_ants)

        elif self.pol_type == 'multi':
            for i in range(len(pols)):
                if self.args.flavors:
                    wu.wedge_flavors(self.args, self.files[i], self.pols[i], self.calfile, self.history, self.freq_range, self.ex_ants)
                elif self.args.time_avg:
                    wu.wedge_timeavg(self.args, self.files[i], self.pols[i], self.calfile, self.history, self.freq_range, self.ex_ants)

        elif self.pol_type == 'single':
            if self.args.delay_avg:
                for file in self.files[0]:
                    wu.wedge_delayavg(file)
            elif self.args.flavors:
                wu.wedge_flavors(self.args, self.files[0], self.pols[0], self.calfile, self.history, self.freq_range, self.ex_ants)
            elif self.args.time_avg:
                wu.wedge_timeavg(self.args, self.files[0], self.pols[0], self.calfile, self.history, self.freq_range, self.ex_ants)
            else:
                wu.wedge_blavg(self.args, self.files[0], self.pols[0], self.calfile, self.history, self.freq_range, self.ex_ants)

if args.step is not None:
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