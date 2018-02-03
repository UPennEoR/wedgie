"""
This script contains all of the logic necessary to call wedge_utils.py to plot or generate polarized pitchforks or wedges.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
Paul La Plante <plaplant_at_sas.upenn.edu>
"""
# Python Standard Modules
import argparse
import os

# Community Made Modules
import numpy as np
import ephem

# UPenn-HERA Moduls
import makeWedge as mW

# HERA Community Modules
from pyuvdata import UVData

# Interactive Development
from IPython import embed


parser = argparse.ArgumentParser()

parser.add_argument("-T",
                    "--filetype",
                    help="Designate whether you will be inputing raw MIRIAD files or npz files, or if you want to use the catalog.",
                    choices=["MIRIAD", "npz", "catalog"],
                    default="MIRIAD")

parser.add_argument("-f",
                    "--filepath",
                    help="Designate the relative path to the folder that contains the files you want to be analyzed.",
                    required=True)

parser.add_argument("-F",
                    "--inputfiles",
                    help="Designate specific files to be analyzed.",
                    nargs='*')
parser.add_argument("-J",
                    "--JDRange",
                    help="Designate an inclusive range of JD to analyze. (e.g.: '2457548.3_2457549.5')")
parser.add_argument("-L",
                    "--LSTRange",
                    help="Designate an inclusive range of LST in hours to analyze. (e.g.: '13.0_15.5' will analyze from 1pm to 3:30pm)")
parser.add_argument("-r",
                    "--LSTrRange",
                    help="Designate an inclusive range of LST in hours to analyze. (e.g.: '2.7_3.7')")

parser.add_argument("-P",
                    "--inputpols",
                    help="Input a comma-delimited list of polarizations to be used in analysis.",
                    default="I,Q,U,V")

parser.add_argument("-X",
                    "--exants",
                    help="Input a comma-delimited list of antennae to exclude from analysis.")

parser.add_argument("-C",
                    "--calfile",
                    help="Enter the calfile to be used for analysis. H0C data only.")

parser.add_argument("-R",
                    "--freqrange",
                    help="Designate with frequency band to analyze.",
                    default='high')

parser.add_argument("-k",
                    "--keyword",
                    help="Designate which file type (by any keyword in the files) you wish to catalog. (e.g.: '.uvcRK', '2457548')")
parser.add_argument("-D",
                    "--datatype",
                    help="Designate abscal or sim data",
                    choices=["abscal", "sim"])

parser.add_argument("--CLEAN",
                    help="Designate the CLEAN tolerance.",
                    type=float,
                    default=1e-3)

args = parser.parse_args()

class Zeus(object):
    def __init__(self, args):
        """Class that handles all the calling and stuff."""
        # Attributes from command line arguments
        self.args = args
        self.filetype = args.filetype
        self.filepath = os.path.abspath(args.filepath)
        self.inputpols = args.inputpols
        self.JDRange = args.JDRange
        self.LSTrRange = args.LSTrRange
        self.LSTRange = args.LSTRange
        self.inputfiles = args.inputfiles
        self.keyword = args.keyword
        self.CLEAN = args.CLEAN
        
        if args.datatype:
            if args.freqrange == 'high':
                self.cosmo = 15.55
            elif args.freqrange == 'low':
                self.cosmo = 16.2
            else:
                self.cosmo = 1.
            self.datatype = args.datatype

        if args.calfile:
            self.calfile = os.path.splitext(os.path.basename(args.calfile))[0]
        else:
            self.calfile = None

        if args.exants:
            self.exants = sorted([int(ant) for ant in args.exants.split(',')])
        else:
            self.exants = list()

        if args.freqrange:
            if args.freqrange == 'high':
                args.freqrange = '580_680'
            elif args.freqrange == 'low':
                args.freqrange = '200_300'
            self.freqrange = [int(freq) for freq in args.freqrange.split('_')]

        # Copy of the initial working directory
        self.cwd = os.getcwd()

        # Polarization specific class attributes
        self.pol_dipole = str()
        self.pol_type = str()

        # File specific class attributes
        self.files_basename = list()
        self.files_filepath = list()
        self.files_keyword = list()
        self.files = dict()
        self.catalog = list()

        # Class constants:
        self.BIN_WIDTH = 0.3
        self.STOKES_POLS = ['I', 'Q', 'U', 'V']
        self.STANDARD_POLS = ['xx', 'xy', 'yx', 'yy']

        # Internal method calls initiated upon instance creation
        self.logic()

    def logic(self):
        """Function to determine which process (catalog, wedge creation, plotting, etc.) to initiate."""
        if self.filetype == 'MIRIAD':
            self.format_pols()
            self.find_files()
            for pol in self.inputpols:
                eris = mW.Eris(self, pol)
                eris.name_npz()
                eris.load_MIRIAD()
                eris.pitchfork()

        elif self.filetype == 'catalog':
            self.catalog_directory()

    def catalog_directory(self):
        # """Catalogs the MIRIAD files in a directory and looks for each file's polarization, JD, and LST,
        # and saves that information in an npz file called 'catalog.npz'."""
        catalog = {}

        array_jd = np.array([], dtype=float)
        array_lstr = np.array([], dtype=float)
        array_lst = np.array([], dtype=float)

        files_filepath = os.listdir(self.filepath)
        files_keyword = sorted([file for file in files_filepath if self.keyword in file])
        if len(files_keyword) == 0:
            raise Exception("There are no files with keyword '{}' in file path '{}'.".format(self.keyword, self.filepath))

        # Cycle through the MIRIAD files in the directory and grab their info
        for file in files_keyword:
            uv = UVData()
            uv.read_miriad(os.path.join(self.filepath, file))
            array_jd = np.array(np.unique(uv.time_array), dtype=float)
            array_lstr = np.array(np.unique(uv.lst_array), dtype=float)
            array_lst = np.array([str(ephem.hours(lstr)) for lstr in array_lstr], dtype=str)

            pol_file = uv.get_pols()[0].lower()

            if pol_file in catalog:
                catalog[pol_file].append({file: {'JD': array_jd, 'LSTr': array_lstr, 'LST': array_lst}})
            else:
                catalog[pol_file] = [{file: {'JD': array_jd, 'LSTr': array_lstr, 'LST': array_lst}}]

        np.savez(os.path.join(self.filepath, 'catalog.npz'), cat=catalog)

        # """Catalogs the MIRIAD files in a directory and looks for each file's polarization, JD, and LST,
        # and saves that information in an npz file called 'catalog.npz'."""
        # catalog = np.array([[], [], [], [], []])

        # files_filepath = os.listdir(self.filepath)
        # files_keyword = sorted([file for file in files_filepath if self.keyword in file])

        # """Maybe include a check here with os.path.splitext()"""

        # # Cycle through the MIRIAD files in the directory and grab their info
        # uv = UVData()
        # for file in files_keyword:
        #     uv.read_miriad(os.path.join(self.filepath, file))
        #     Ntimes = uv.Ntimes
        #     pol_file = uv.get_pols()[0].lower()
        #     array_lstr = np.unique(uv.lst_array)
        #     array_lst = np.array([str(ephem.hours(lstr)) for lstr in array_lstr])
        #     array_jd = np.unique(uv.time_array)
        #     array_file = np.array([file] * Ntimes)
        #     array_pol = np.array([pol_file] * Ntimes)

        #     file_array = np.array((array_file, array_pol, array_jd.astype(float), array_lstr.astype(float), array_lst.astype(str)))
        #     catalog = np.concatenate((catalog, file_array), axis=1)
        # np.savez(os.path.join(self.filepath, 'catalog.npz'), cat=catalog.T)

    def find_files(self):
        """Properly format the given files specified directly or by JD, LST (hours (coming soon)), or LST (radiams)."""
        # Change directory to the directory with catalog.npz
        os.chdir(self.filepath)

        # Try to open catalog.npz
        try:
            self.catalog = np.load('catalog.npz')['cat'].tolist()
        except IOError:
            raise Exception("There is no catalog (catalog.npz) in the specified path: %s" %self.filepath)

        if self.JDRange:
            JDRange_start, JDRange_stop = [float(x) for x in self.JDRange.split('_')]
        elif self.LSTrRange:
            LSTrRange_start, LSTrRange_stop = [float(x) for x in self.LSTrRange.split('_')]
        elif self.LSTRange:
            LSTRange_start, LSTRange_stop = [ephem.hours(x) for x in self.LSTRange.split('_')]

        for pol in self.STANDARD_POLS:
            if (pol not in self.pol_dipole) and (pol in self.catalog):
                del self.catalog[pol]
                continue

            for index, file in enumerate(self.catalog[pol]):
                key = self.catalog[pol][index].keys()[0]
                if self.JDRange:
                    indices = np.where(np.logical_and(self.catalog[pol][index][key]['JD'] > JDRange_start,
                                                      self.catalog[pol][index][key]['JD'] < JDRange_stop))[0]
                elif self.LSTrRange:
                    indices = np.where(np.logical_and(self.catalog[pol][index][key]['LSTr'] > LSTrRange_start,
                                                      self.catalog[pol][index][key]['LSTr'] < LSTrRange_stop))[0]
                elif self.LSTRange:
                    indices = np.where(np.logical_and(self.catalog[pol][index][key]['LSTr'] > LSTRange_start,
                                                      self.catalog[pol][index][key]['LSTr'] < LSTRange_stop))[0]
                else:
                    raise Exception("You must indicate which type of range you would like to use (JDRange, LSTrRange, LSTRange).")

                if len(indices):
                    if pol not in self.files:
                        self.files[pol] = [key]
                    else:
                        self.files[pol].append(key)
                else:
                    del self.catalog[pol][index][key]

            self.catalog[pol] = [x for x in self.catalog[pol] if len(x) != 0]

        embed()

        # # Remove files from catalog that don't have a pol specified by args.pol
        # # Necessary for all types of file finding
        # row_index = 0
        # for pol in self.catalog[:, 1]:
        #     if not np.any(pol in self.pol_dipole):
        #         self.catalog = np.delete(self.catalog, row_index, axis=0)
        #         row_index -= 1
        #     row_index += 1

        # # Check if the user wants to find files based on JD, LST (hours), or LST (radians)
        # if self.JDRange or self.LSTRange or self.LSTrRange:
        #     # Find indices of files in catalog that are within the specified range.
        #     if self.JDRange:
        #         JDRange_start, JDRange_stop = [float(x) for x in self.JDRange.split('_')]
        #         indices = np.where(np.logical_and(self.catalog[:, 2].astype(float) >= JDRange_start,
        #                                           self.catalog[:, 2].astype(float) <= JDRange_stop))[0]
        #     elif self.LSTrRange:
        #         LSTrRange_start, LSTrRange_stop = [float(x) for x in self.LSTrRange.split('_')]
        #         indices = np.where(np.logical_and(self.catalog[:, 3].astype(float) >= LSTrRange_start,
        #                                           self.catalog[:, 3].astype(float) <= LSTrRange_stop))[0]
        #     self.catalog = np.take(self.catalog, indices, axis=0)

        # # Check if the user wants to input their own files directly
        # elif self.inputfiles:
        #     # Take just the basename of the given files
        #     # (This necessitates not including the ending '/' when inputting files)
        #     self.files_basename = [os.path.basename(file) for file in self.inputfiles]

        #     # Remove files from the catalog that are not in args.inputfiles
        #     row_index = 0
        #     for file in self.catalog[:, 0]:
        #         if not np.any(file in self.files_basename):
        #             self.catalog = np.delete(self.catalog, row_index, axis=0)
        #             row_index -= 1
        #         row_index += 1

        # # Organizes the good files into self.files, a dictionary: {pol: [file1, file2, ...], ...}
        # files = self.catalog[:, :2]
        # self.files = {pol: [] for pol in self.pol_dipole}
        # for row in files:
        #     if row[0] not in self.files[row[1]]:
        #         self.files[row[1]].append(row[0])

    def format_pols(self):
        """Format the polarizations, e.g.: translating from IQUV to xx,xy,yx,yy"""
        self.inputpols = self.inputpols.split(',')
        self.pol_dipole = []

        # Check if the given pols from args.input pols are stokes or standard
        if any(pol in self.inputpols for pol in self.STOKES_POLS):
            self.pol_type = 'stokes'
            if ('I' in self.inputpols) or ('Q' in self.inputpols):
                self.pol_dipole.extend(['xx', 'yy'])
            if ('U' in self.inputpols) or ('V' in self.inputpols):
                self.pol_dipole.extend(['xy', 'yx'])
        elif any(pol in self.inputpols for pol in self.STANDARD_POLS):
            self.pol_type = 'standard'
            self.pol_dipole = self.inputpols[:]
        else:
            raise Exception("You provided nonsensical polarization types: %s" %self.inputpols)

zeus = Zeus(args)
# Name class in new wedge_utils 'Eris'

# parser.add_argument("-A",
#                     "--AltAnalysis",
#                     help="Indicate that you want to use an alternative wedge creation method.",
#                     action="store_true",
#                     default=False)
# parser.add_argument("-b",
#                     "--blavg",
#                     help="Toggle baseline averaging only.")
# parser.add_argument("-f",
#                     "--flavors",
#                     help="Toggle splitting wedgeslices into a per slope, per baseline basis.",
#                     action="store_true")
# parser.add_argument("-l",
#                     "--BaselineNum",
#                     help="Input one baseline type (1-8 for H0C).",
#                     type=int)

# parser.add_argument("-N",
#                     "--npz_operations",
#                     help="Various functions that work on npz files.",
#                     choices=["pspec_diff", "wedge_diff", "combine", "delayavg"])

# parser.add_argument("-p",
#                     "--path",
#                     help="Enter the path where you want save the files.")