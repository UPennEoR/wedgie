"""
This program compliments getWedge.py by being able to plot the npz files generated from getWedge.py.

Author: Austin Fox Fortino ,fortino@sas.upenn.edu
"""
import argparse, wedge_utils

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filenames', help='Input a list of filenames to be analyzed.', nargs='*', required=True)
parser.add_argument('-s', '--single_plot', help='Plot a single plot from supplied npz files.', action='store_true')
parser.add_argument('-m', '--multi_plot', help='Plot 4 plots at once from supplied npz files.', action='store_true')
parser.add_argument('-a', '--avg_plot', help='Plots average value inside and outside wedge per files analyzed.',action='store_true')
parser.add_argument('-d', '--delay_plot', help='Plot a single plot from supplied delayavg npz file', action='store_true')
parser.add_argument('-M','--multi_delayplot', help ='Plot multiple delay average plots', action='store_true')
parser.add_argument('-o', '--plot_1D', help="Plot (optional: specified as comma delimited list) baselines' wedges on a 1D plot from supplied npz file", default=None, const='all', nargs='?', action='store')

args = parser.parse_args()

if args.plot_1D is not None:
    if args.plot_1D == 'all':
        baselines = []
    else:
        baselines = [int(x) for x in args.plot_1D.split(',')]
    for filename in args.filenames:
        wedge_utils.plot_1D(filename, baselines)

if args.delay_plot:
    for filename in args.filenames:
        wedge_utils.plot_delayavg(filename)

"""if args.multi_delayplot:
    for filename in args.filenames:
        wedge_utils.plot_multi_delayavg(filename, args.save_path)"""

if args.single_plot:
    for filename in args.filenames:
        if filename.split('.')[-2] == 'timeavg':
            wedge_utils.plot_timeavg(filename)
        elif filename.split('.')[-2] == 'blavg':
            wedge_utils.plot_blavg(filename)

elif args.multi_plot:
    wedge_utils.plot_multi_timeavg(args.filenames)

elif args.avg_plot:
    wedge_utils.plot_avgs(args.filenames)
