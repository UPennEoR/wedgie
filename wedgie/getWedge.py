"""
This code contains an importable function getWedge intended to be used on HERA data.  
It takes the delay transform across all visibilities and creates a wedge plot 
comprised of subplots averaged over antennae of same baselengths.

It requires a HERA data filename string such as:
    "path/to/zen.2457700.47314.xx.HH.uvcRR"

wedge_utils.py, and the calfile should be in the PYTHONPATH

Co-Author: Paul Chichura <pchich@sas.upenn.edu>
Co-Author: Amy Igarashi <igarashiamy@gmail.com>
Co-Author: Austin Fox Fortino <fortino@sas.upenn.edu>
Date Created: 6/21/2017
"""
import argparse, wedge_utils
# from pprint import pprint
# from IPython import embed

def callers(args):
    # Sets up the frequency range to analyze over.
    freq_range = (int(args.freq_range.split('_')[0]), int(args.freq_range.split('_')[1]))

    # Set up the list of pols for later use
    if args.pol == 'stokes':
        pols = ['xx','xy','yx','yy']
    else:
        pols = args.pol.split(",")

    # Generate the list of filenames depending on what polarizations were chosen.
    filenames = []
    for pol in pols:
        #make a list of all filenames for each polarization
        pol_filenames = []
        for filename in args.filenames:
            #replace polarization in the filename with pol we want to see
            filepol = filename.split('.')[-3]
            new_filename = filename.split(filepol)[0]+pol+filename.split(filepol)[1]
            #append it if it's not already there
            if not any(new_filename in s for s in pol_filenames):
                pol_filenames.append(new_filename)
        filenames.append(pol_filenames)

    # Generates the history dictionary that gets attached to each npz file.
    history = vars(args)

    # Checks if any antennae were specified by --ex_ants.
    if not args.ex_ants is None:
        ex_ants_list = map(int, args.ex_ants.split(','))
    else:
        ex_ants_list = []


    #calculate and get the names of the npz files for
    if args.pol == 'stokes' and args.blavg:
        print freq_range
        print 'that fre;a;kfdasf'
        print freq_range[0]
        print freq_range[1]
        npz_names = wedge_utils.wedge_stokes(filenames, args.calfile.split('.')[0], history, freq_range, ex_ants_list, blavg=True)

    if args.pol == 'stokes' and args.bl_num:
        npz_names = wedge_utils.wedge_stokes(filenames, args.calfile.split('.')[0], args.bl_num, history, freq_range, ex_ants_list, bltype=True)

    elif args.pol == 'stokes':
        #calculate and get the names of the npz files
        npz_names = wedge_utils.wedge_stokes(filenames, args.calfile.split('.')[0], history, freq_range, args.flavors, ex_ants_list)

    elif args.delay_avg and (len(pols) == 1):
        for filename in args.filenames:
            wedge_utils.wedge_delayavg(filename)

    elif len(pols) == 1:
        if args.flavors:
            npz_name = wedge_utils.wedge_flavors(args.filenames, args.pol, args.calfile.split('.')[0], history, freq_range, ex_ants_list)
        elif args.time_avg:
            npz_name = wedge_utils.wedge_timeavg(args.filenames, args.pol, args.calfile.split('.')[0], history, freq_range, ex_ants_list)
        else:
            npz_name = wedge_utils.wedge_blavg(args.filenames, args.pol, args.calfile.split('.')[0], history, ex_ants_list)

    elif len(pols) > 1:
        npz_names = []

        for i in range(len(pols)):
            if args.flavors:
                npz_names.append(wedge_utils.wedge_flavors(filenames[i], pols[i], args.calfile.split('.')[0], history, freq_range, ex_ants_list))
            elif args.time_avg:
                npz_names.append(wedge_utils.wedge_timeavg(filenames[i], pols[i], args.calfile.split('.')[0], history, freq_range, ex_ants_list))

parser = argparse.ArgumentParser()
parser.add_argument('-F', '--filenames', nargs='*', required=True, help='Input a list of filenames to be analyzed.')
parser.add_argument('-C', '--calfile', default='hsa7458_v001', help='Input the calfile to be used for analysis.')
parser.add_argument('-P', '--pol', default='stokes', help='Input a comma-delimited list of polatizations to plot.')
parser.add_argument('-f', '--flavors', default=False, action='store_true', help='Toggle splitting wedgeslices into a per slope per baseline basis.')
parser.add_argument('-t', '--time_avg', default=True, action='store_false', help='Toggle off time averaging.')
parser.add_argument('-x', '--ex_ants', type=str, help='Input a comma-delimited list of antennae to exclude from analysis.')
parser.add_argument('-s', '--step', type=int, help='Toggle file stepping.')
parser.add_argument('-r', '--freq_range', default='0_1023', help='Input a range of frequency channels to use separated by an underscore: "550_650"')
parser.add_argument('-a', '--stair', help='Compute npz files for 1 file, then 2 files, then 3 files, ...  to see how wedges change over multiple files averaged together.', action='store_true')
# parser.add_argument('-l', '--load', default=3, type=int, help='How many operations do you want to run at once. Used only in conjuction with --step.') XXX Deprecated
parser.add_argument('-d','--delay_avg', help="sfsdfasdfsf", action="store_true")
parser.add_argument('-b', '--blavg', action='store_true', default=False, help='Toggle blavg for stokes.')
parser.add_argument('-l','--bl_num', type=int, help='Toggle bltype and input 1 baseline type.')
args = parser.parse_args()

# If --step or --stair is used, they are executed here.
if not args.step is None:
    step, files = args.step, args.filenames
    del args.step, args.filenames

    files_xx = [file for file in files if 'xx' in file]
    num_files_xx = len(files_xx)

    for index in range(0, num_files_xx, step):
        args.filenames = files_xx[index : index + step]
        callers(args)
    print 'Step program complete.'

elif not args.stair is None:
    files = args.filenames
    del args.filenames

    files_xx = [file for file in files if 'xx' in file]
    num_files_xx = len(files_xx)

    for index in range(1, num_files_xx):
        args.filenames = files_xx[0: index]
        callers(args)
    print 'Stair program complete.'

else:
    callers(args)