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
import argparse, wedge_utils, os, pprint

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filenames', help='Input a list of filenames to be analyzed.', nargs='*', required=True)
parser.add_argument('-c', '--calfile', help='Input the calfile to be used for analysis.', required=True)
parser.add_argument('--pol', help='Input a comma-delimited list of polatizations to plot.', required=True)
parser.add_argument('-t', '--time_avg', help='Toggle time averaging.', action='store_true')
parser.add_argument('-x', '--ex_ants', help='Input a comma-delimited list of antennae to exclude from analysis.', type=str)
parser.add_argument('-s', '--step', help='Toggle file stepping.', action='store')
args = parser.parse_args()

if not args.step is None:
    opts = ["-c " + args.calfile, "--pol " + args.pol]
    if args.time_avg:
        opts.append("-t")
    if args.plot:
        opts.append("-p")
    elif args.only_plot:
        opts.append("-o")
    if not args.ex_ants is None:
        opts.append("-x={}".format(args.ex_ants))

    files_all = [file for file in args.filenames if 'xx' in file]

    step = int(args.step)
    for file_index in range(0, len(files_all), step):
        cmd = opts + ["-f"] + files_all[file_index : file_index + step]
        
        print "I just executed the following arguments:"
        pprint.pprint(cmd)
        
        os.system("python2.7 getWedge.py {}".format(" ".join(cmd)))
        
        print
    quit()

pols = args.pol.split(",")

if not args.ex_ants is None:
    ex_ants_list = map(int, args.ex_ants.split(','))
else:
    ex_ants_list = []

if pols == ['stokes']:
    filenames = []
    for pol in ['xx','xy','yx','yy']:
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
    #calculate and get the names of the npz files
    npz_names = wedge_utils.wedge_stokes(filenames, args.calfile.split('.')[0], ex_ants_list)

elif len(pols) == 1:
    if args.time_avg:
        npz_name = wedge_utils.wedge_timeavg(args.filenames, args.pol, args.calfile.split('.')[0], ex_ants_list)
    else:
        npz_name = wedge_utils.wedge_blavg(args.filenames, args.pol, args.calfile.split('.')[0], ex_ants_list)

elif len(pols) > 1:
    npz_names = []
    for pol in pols:
        #make a list of all filenames for each polarization
        filenames = []
        for filename in args.filenames:
            #replace polarization in the filename with pol we want to see
            filepol = filename.split('.')[-3]
            new_filename = filename.split(filepol)[0]+pol+filename.split(filepol)[1]
            #append it if it's not already there
            if not any(new_filename in s for s in filenames):
                filenames.append(new_filename)
        npz_names.append(wedge_utils.wedge_timeavg(filenames, pol, args.calfile.split('.')[0], ex_ants_list))