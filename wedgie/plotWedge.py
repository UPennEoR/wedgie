"""
This program compliments getWedge.py by being able to plot the npz files generated from getWedge.py.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""
import argparse, wedge_utils

parser = argparse.ArgumentParser()
parser.add_argument('-F', '--filenames', help='Input a list of filenames to be analyzed.', nargs='*', required=True)
parser.add_argument('-s', '--single_plot', help='Plot a single plot from supplied npz files.', action='store_true')
parser.add_argument('-m', '--multi_plot', help='Plot 4 plots at once from supplied npz files.', action='store_true')
parser.add_argument('-f', '--flavors_plot', action='store_true')
parser.add_argument('-b', '--multi_bl_plot', help='Plot 4 plots for blavg.', action='store_true')
parser.add_argument('-a', '--avg_plot', help='Plots average value inside and outside wedge per files analyzed.', action='store_true')
parser.add_argument('-d', '--delay_plot', help='Plot a single plot from supplied delayavg npz file', action='store_true')
parser.add_argument('-l', '--plot_bltype', help='Plot non-averaged plots for given bltype file.', action='store_true')
parser.add_argument('-o', '--plot_1D', help="Plot (optional: specified as comma delimited list) baselines' wedges on a 1D plot from supplied npz file", default=None, const='all', nargs='?', action='store')
args = parser.parse_args()

if (args.plot_1D is not None) and not args.multi_plot:
    if args.plot_1D == 'all':
        baselines = []
    else:
        baselines = [int(x) for x in args.plot_1D.split(',')]
    for filename in args.filenames:
        wedge_utils.plot_1D(filename, baselines)

elif args.delay_plot:
    for filename in args.filenames:
        wedge_utils.plot_delayavg(filename)

elif args.single_plot:
    for filename in args.filenames:
        if filename.split('.')[-2] == 'timeavg':
            wedge_utils.plot_timeavg(filename)
        elif filename.split('.')[-2] == 'blavg':
            wedge_utils.plot_blavg(filename)

elif args.plot_bltype:
    for filename in args.filenames:
        wedge_utils.plot_bltype(filename)

elif args.multi_plot:
    if args.plot_1D is not None:
        if args.plot_1D == 'all':
            baselines = []
        else:
            baselines = [int(x) for x in args.plot_1D.split(',')]
        wedge_utils.plot_multi_1D(args.filenames, baselines)
    else:
        if len(args.filenames) > 1:
            for index in range(0, len(args.filenames), 4):
                wedge_utils.plot_multi_timeavg(args.filenames[index:index+4])
        else:
            wedge_utils.plot_timeavg(args.filenames)

elif args.multi_bl_plot:
    wedge_utils.plot_multi_blavg(args.filenames)

elif args.avg_plot:
    wedge_utils.plot_avgs(args.filenames, rng=(0,len(args.filenames)))

elif args.flavors_plot:
    if len(args.filenames) > 1:
        for index in range(0, len(args.filenames), 4):
            wedge_utils.plot_multi_flavors(args.filenames[index:index+4])
    else:
        wedge_utils.plot_flavors(args.filenames)