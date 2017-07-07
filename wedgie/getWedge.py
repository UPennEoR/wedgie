"""
This code contains an importable function getWedge intended to be used on HERA
telescope data.  It takes the delay transform across all visibilities and
creates a wedge plot comprised of subplots averaged over antennae of same
baselengths.

It requires a HERA data filename string such as:
    "zen.2457700.47314.xx.HH.uvcRR"
It is reliant on Austin Fox Fortino's baseline_lengths.py
Can be imported as a function or used in the terminal with argument filename.

Co-Author: Paul Chichura <pchich@sas.upenn.edu>
Co-Author: Amy Igarashi <igarashiamy@gmail.com>
Date: 6/21/2017
"""

import argparse
import wedge_utils
import glob

#get filename from command line argument
parser = argparse.ArgumentParser()
parser.add_argument("filenames", help="your HERA data file(s)", nargs="*")
parser.add_argument("calfile", help="your calfile")
parser.add_argument("pol", help="comma-delimited list of pols to plot for filenames")
parser.add_argument("--time_avg", help="Toggle time averaging", action="store_true")
parser.add_argument("--ex_ants", type=str, help='comma-delimited list of antennae to exclude.')
parser.add_argument("--plot", help="toggle plotting the data in addition to saving as a .npz", action="store_true")
parser.add_argument("--only_plot", help="call just the plot functions for filenames=npz name", action="store_true")
args=parser.parse_args()

pols = args.pol.split(",")

# format ex_ants argument for intake
if not args.ex_ants is None:
    ex_ants_list = map(int, args.ex_ants.split(','))

if args.only_plot and (len(pols) == 1 ):
    for filename in args.filenames:
        if filename.split('.')[-2] == 'timeavg':
            wedge_utils.plot_timeavg(filename)
        elif filename.split('.')[-2] == 'blavg':
            wedge_utils.plot_blavg(filename)

elif (not args.only_plot) and (len(pols) == 1):
    #make wedge
    if not args.ex_ants is None:
        if args.time_avg:
            npz_name = wedge_utils.wedge_timeavg(args.filenames, args.pol, args.calfile.split('.')[0], ex_ants_list)
            if args.plot: 
                wedge_utils.plot_timeavg(npz_name)
        else:
            npz_name = wedge_utils.wedge_blavg(args.filenames, args.pol, args.calfile.split('.')[0], ex_ants_list)
            if args.plot:
                wedge_utils.plot_blavg(npz_name)

    else:
        if args.time_avg:
            npz_name = wedge_utils.wedge_timeavg(args.filenames, args.pol, args.calfile.split('.')[0])
            if args.plot: 
                wedge_utils.plot_timeavg(npz_name)
        else:
            npz_name = wedge_utils.wedge_blavg(args.filenames, args.pol, args.calfile.split('.')[0])
            if args.plot:
                wedge_utils.plot_blavg(npz_name)


elif (len(pols) > 1):

    npz_names = []
    
    if args.only_plot:
        npz_names = args.filenames
    else:
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

    if args.plot or args.only_plot:
        wedge_utils.plot_multi_timeavg(npz_names)

    


