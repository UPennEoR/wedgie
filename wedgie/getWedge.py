"""
This code contains an importable function getWedge intended to be used on HERA data.  
It takes the delay transform across all visibilities and creates a wedge plot 
comprised of subplots averaged over antennae of same baselengths.

Arguments:
-f path/to/FILENAME [path/to/FILENAME [...]]
-c=CALFILE
--pol=[stokes], [xx] [xy] [yx] [yy]
-t
-x=antenna,antenna,...
-s=STEP

It requires a HERA data filename string such as:
    "path/to/zen.2457700.47314.xx.HH.uvcRR"

wedge_utils.py, and the calfile should be in the PYTHONPATH

Co-Author: Paul Chichura <pchich@sas.upenn.edu>
Co-Author: Amy Igarashi <igarashiamy@gmail.com>
Co-Author: Austin Fox Fortino <fortino@sas.upenn.edu>
Date Created: 6/21/2017
"""
import argparse, wedge_utils, os, pprint, threading, sys

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filenames', nargs='*', required=True, help='Input a list of filenames to be analyzed.')
parser.add_argument('-c', '--calfile', help='Input the calfile to be used for analysis.')
parser.add_argument('-p', '--pol', help='Input a comma-delimited list of polatizations to plot.')
parser.add_argument('-t', '--time_avg', action='store_true', help='Toggle time averaging.')
parser.add_argument('-x', '--ex_ants', type=str, help='Input a comma-delimited list of antennae to exclude from analysis.')
parser.add_argument('-s', '--step', type=int, help='Toggle file stepping.')
parser.add_argument('-r', '--range', help='Supply a range of times throughout a day to process.')
parser.add_argument('-F', '--freq', default='0_1024')
parser.add_argument('-d','--delay_avg', help="sfsdfasdfsf", action="store_true")
parser.add_argument('-b', '--blavg', action='store_true', default=False, help='Toggle blavg for stokes.')
parser.add_argument('-l','--bl_num', type=int, help='Toggle bltype and input 1 baseline type.')

args = parser.parse_args()

freq_range = (int(args.freq.split('_')[0]), int(args.freq.split('_')[1]))


files = args.filenames[:] 

if not args.range == None:
    range_start = args.range.split('_')[0]
    range_end = args.range.split('_')[1]

    for file in files[:]:
        time = file.split('.')[-4]
        if time < range_start:
            files.remove(file)
        elif time > range_end:
            files.remove(file)

if not args.step is None:
    wedge_utils.step(sys.argv, args.step, files, args.filenames)
    quit()

if args.pol == 'stokes':
    pols = ['xx','xy','yx','yy']
else:
    pols = args.pol.split(",")

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

history = {
    'filenames': filenames, 
    'calfile': args.calfile, 
    'pol': args.pol, 
    'time_avg': args.time_avg, 
    'ex_ants': args.ex_ants, 
    'step': args.step,
    'range': args.range,
    'freq_range': args.freq
    }

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
    npz_names = wedge_utils.wedge_stokes(filenames, args.calfile.split('.')[0], history, freq_range, ex_ants_list)

elif args.delay_avg and (len(pols) == 1):
    for filename in args.filenames:
        wedge_utils.wedge_delayavg(filename)

elif len(pols) == 1:
    if args.time_avg:
        npz_name = wedge_utils.wedge_timeavg(args.filenames, args.pol, args.calfile.split('.')[0], history, freq_range, ex_ants_list)
    else:
        npz_name = wedge_utils.wedge_blavg(args.filenames, args.pol, args.calfile.split('.')[0], history, freq_range, ex_ants_list)

elif len(pols) > 1:
    npz_names = []
    if args.time_avg:
        for i in range(len(pols)):
            npz_names.append(wedge_utils.wedge_timeavg(filenames[i], pols[i], args.calfile.split('.')[0], history, freq_range,
ex_ants_list))
    else:
            npz_names.append(wedge_utils.wedge_blavg(filenames[i], pols[i], args.calfile.split('.')[0], history, freq_range, ex_ants_list))

