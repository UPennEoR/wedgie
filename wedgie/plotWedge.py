"""
This program compliments getWedge.py by being able to plot the npz files generated from getWedge.py.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""
import argparse
import wedge_utils as wu

parser = argparse.ArgumentParser()
parser.add_argument('-F', '--filenames', help='Input a list of filenames to be analyzed.', nargs='*', required=True)
parser.add_argument('-t', '--timeavg', action='store_true')
parser.add_argument('-b', '--blavg', action='store_true')
parser.add_argument('-f', '--flavors', action='store_true')
parser.add_argument('-s', '--single', help='Plot a single plot from supplied npz files.', action='store_true')
parser.add_argument('-m', '--multi', help='Plot 4 plots at once from supplied npz files.', action='store_true')

parser.add_argument('-a', '--avg', help='Plots average value inside and outside wedge per files analyzed.', action='store_true')
parser.add_argument('-d', '--delay', help='Plot a single plot from supplied delayavg npz file', action='store_true')
parser.add_argument('-l', '--bl_type', help='Plot non-averaged plots for given bltype file.', action='store_true')
parser.add_argument('-o', '--one_D', help="Plot (optional: specified as comma delimited list) baselines' wedges on a 1D plot from supplied npz file", default=None, const='all', nargs='?', action='store')
parser.add_argument('-P', '--path', default='./')
args = parser.parse_args()


class Graph(object):
    def __init__(self, args):
        self.args = args
        self.files = args.filenames
        self.path = args.path

        self.MISSING_TAG_ERR = "You must specify which type of Wedge to plot."

    def logic(self):
        if self.args.single:
            for file in self.files:
                if self.args.timeavg:
                    wu.plot_timeavg(file, self.path)
                elif self.args.blavg:
                    wu.plot_blavg(file, self.path)
                elif self.args.flavors:
                    wu.plot_flavors(file, self.path)
                elif self.args.bl_type:
                    wu.plot_bltype(file)
                elif self.args.delay:
                    wu.plot_delayavg(file)
                elif self.args.one_D:
                    if self.args.one_D == 'all':
                        baselines = []
                    else:
                        baselines = [int(x) for x in self.args.one_D.split(',')]
                    wu.plot_1D(file, baselines)
                else:
                    raise Exception(self.MISSING_TAG_ERR)

        elif self.args.multi:
            for index in range(0, len(self.files), 4):
                files = self.files[index:index+4]
                if self.args.timeavg:
                    wu.plot_timeavg_multi(files, self.path)
                elif self.args.flavors:
                    wu.plot_flavors_multi(files, self.path)
                elif self.args.one_D:
                    wu.plot_multi_1D
                else:
                    raise Exception(self.MISSING_TAG_ERR)

        elif self.args.avg:
            wu.plot_avgs(self.files)


graph = Graph(args)
graph.logic()

# if (args.plot_1D is not None) and not args.multi_plot:
#     if args.plot_1D == 'all':
#         baselines = []
#     else:
#         baselines = [int(x) for x in args.plot_1D.split(',')]
#     for filename in args.filenames:
#         wedge_utils.plot_1D(filename, baselines)

# elif args.delay_plot:
#     for filename in args.filenames:
#         wedge_utils.plot_delayavg(filename)

# elif args.single_plot:
#     for filename in args.filenames:
#         if filename.split('.')[-2] == 'timeavg':
#             wedge_utils.plot_timeavg(filename)
#         elif filename.split('.')[-2] == 'blavg':
#             wedge_utils.plot_blavg(filename)

# elif args.plot_bltype:
#     for filename in args.filenames:
#         wedge_utils.plot_bltype(filename)

# elif args.multi_plot:
#     if args.plot_1D is not None:
#         if args.plot_1D == 'all':
#             baselines = []
#         else:
#             baselines = [int(x) for x in args.plot_1D.split(',')]
#         wedge_utils.plot_multi_1D(args.filenames, baselines)
#     else:
#         if len(args.filenames) > 1:
#             for index in range(0, len(args.filenames), 4):
#                 wedge_utils.plot_multi_timeavg(args.filenames[index:index+4])
#         else:
#             wedge_utils.plot_timeavg(args.filenames)

# elif args.multi_bl_plot:
#     wedge_utils.plot_multi_blavg(args.filenames)

# elif args.avg_plot:
#     wedge_utils.plot_avgs(args.filenames, rng=(0,len(args.filenames)))

# elif args.flavors_plot:
#     if len(args.filenames) > 1:
#         for index in range(0, len(args.filenames), 4):
#             wedge_utils.plot_multi_flavors(args.filenames[index:index+4])
#     else:
#         wedge_utils.plot_flavors(args.filenames)
