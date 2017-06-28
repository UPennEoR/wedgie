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
#get filename from command line argument
parser = argparse.ArgumentParser()
parser.add_argument("filenames", help="your HERA data file(s)", nargs="*")
parser.add_argument("calfile", help="your calfile")
parser.add_argument("pol", help="polarization of data")
parser.add_argument("--time_avg", help="Toggle time averaging", action="store_true")
parser.add_argument("--ex_ants", type=str, help='comma-delimited list of antennae to exclude.')
args=parser.parse_args()

# format ex_ants argument for intake
ex_ants_list = map(int, args.ex_ants.split(','))

#make wedge
if args.time_avg:
    wedge_utils.plot_wedge_timeavg(args.filenames, args.pol, args.calfile.split(".")[0], ex_ants_list) 
else:
    wedge_utils.plot_wedge_blavg(args.filenames, args.pol, args.calfile.split(".")[0], ex_ants_list) 

