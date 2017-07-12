"""
This program compliments getWedge.py by being able to plot the npz files generated from getWedge.py.

Author: Austin Fox Fortino <fortino@sas.upenn.edu>
Created: July 11, 2017
Last Updated: July 11, 2017
"""
import argparse, wedge_utils, os, pprint

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filenames', help='Input a list of filenames to be analyzed.', nargs='*', required=True)
parser.add_argument('-p', '--path', help='Path to save destination for png files.', default='./')
parser.add_argument('-s', '--single_plot', help='Plot a single plot from supplied npz files.', action='store_true')
parser.add_argument('-m', '--multi_plot', help='Plot 4 plots at once from supplied npz files.', action='store_true')
parser.add_argument('-d', '--delay_plot', help='Plot a single plot from supplied delayavg npz file', action='store_true')
args = parser.parse_args()

if args.delay_plot:
    for filename in args.filenames:
        wedge_utils.plot_delayavg(filename)

if args.single_plot:
    for filename in args.filenames:
        if 'timeavg' in filename:
            wedge_utils.plot_timeavg(filename, args.path)
        elif 'blavg' in filename:
            wedge_utils.plot_blavg(filename, args.path)

elif args.multi_plot:
    wedge_utils.plot_multi_timeavg(args.filenames, args.path)
