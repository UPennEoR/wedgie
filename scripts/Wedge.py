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
import glob

# Python Community Science Modules
import numpy as np
import astropy.units as u
import astropy.coordinates as aco

# UPenn-HERA Modules
import makeWedge as mW

# HERA Collaboration Modules
from pyuvdata import UVData

# For Interactive Development
from IPython import embed

parser = argparse.ArgumentParser()

# Used for every operation
parser.add_argument("-f",
    "--filepath",
    help="Designate the relative path to the folder that contains the files you want to be analyzed.",
    required=True)

# File designating arguments
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

# Designate polarization types to be analyzed
parser.add_argument("-P",
    "--inputpols",
    help="Input a comma-delimited list of polarizations to be used in analysis.",
    default="I,Q,U,V")

# MIRIAD analysis options
parser.add_argument("-X",
    "--exants",
    help="Input a comma-delimited list of antennae to exclude from analysis.",
    choices=[None, "H0C", "H1C"],
    default=None)
parser.add_argument("-C",
    "--calfile",
    help="Enter the calfile to be used for analysis. H0C data only.",
    choices=[None, "H0C"],
    default=None)
parser.add_argument("-R",
    "--freqrange",
    help="Designate with frequency band to analyze.",
    default=None)

# Time-Averaging Tag
parser.add_argument("-A",
    "--tag_avg",
    help="Time average or non-time averaged plots",
    choices=["timeavg", "blavg"],
    default="timeavg")
# Depth of analysis tag
parser.add_argument("-D",
    "--tag_depth",
    choices=["bl", "sl", "ant"],
    default="bl")
# Pitchfork or Wedge tag
parser.add_argument("-W",
    "--tag_wedge",
    choices=["pf", "w", "1d", "1d"],
    default="pf")
# Unit type tag
parser.add_argument("-U",
    "--tag_unit",
    help="Designate the units to be used",
    choices=["std", "cosmo"],
    default="std")

# Save path
parser.add_argument("-p",
    "--path",
    help="Enter the path where you want save the files.",
    default=".")
# Catalog/Extension
parser.add_argument("-e",
    "--extension",
    help="Designate which file type (by its extension) you wish to catalog. (e.g.: 'uvcRK')",
    action="append")
# Print Catalog
parser.add_argument("--printc",
    help="Print the directory catalog.",
    action="store_true")
# Analysis constants
parser.add_argument("--CLEAN",
    help="Designate the CLEAN tolerance.",
    type=float,
    default=1e-3)
parser.add_argument("-T",
    "--title",
    type=int,
    default=1)

# Stepping through files
parser.add_argument("-S",
    "--step",
    type=int,
    default=0)

args = parser.parse_args()

class Zeus(object):
    def __init__(self, args):
        """Class that handles all the calling and stuff."""
        # Attributes from command line arguments
        self.args = args
        self.filepath = os.path.abspath(args.filepath)

        # Excluded Antennae Formatting
        if args.exants == None:
            self.exants = list()
        elif args.exants == "H0C":
            self.exants = [22, 43, 80, 81]
        elif args.exants == "H1C":
            self.exants = [50, 98]
        else:
            self.exants = sorted([int(ant) for ant in args.exants.split(',')])

        # Calibration File Formatting
        if args.calfile == None:
            self.calfile = None
        elif args.calfile == "H0C":
            self.calfile = "hsa7458_v001"
        else:
            self.calfile = os.path.splitext(os.path.basename(args.calfile))[0]

        # Frequency Range Formatting
        if args.freqrange == None:
            args.freqrange = "0_1024"
        elif args.freqrange == "high":
            args.freqrange = '580_680'
        elif args.freqrange == 'low':
            args.freqrange = '200_300'
        self.freqrange = [int(freq) for freq in args.freqrange.split('_')]

        # Plotting Tags
        self.tag_avg = args.tag_avg
        self.tag_depth = args.tag_depth
        self.tag_unit = args.tag_unit
        self.tag_wedge = args.tag_wedge

        # Turn plot titles on/off
        self.title = bool(args.title)

        # Extension used for cataloging directories
        if args.extension:
            self.extension = [ext.split('.')[-1] for ext in args.extension]

        # Print catalog boolean attribute
        self.printc = args.printc

        # Stepping attribute
        self.step = args.step
        self.i = int()

        # Copy of the initial working directory
        self.cwd = os.getcwd()

        # Save path
        self.path = os.path.abspath(args.path)

        # Polarization specific class attributes
        self.inputpols = args.inputpols
        self.pol_dipole = str()
        self.pol_type = str()

        # File specific class attributes
        self.inputfiles = args.inputfiles
        self.JDRange = args.JDRange
        self.LSTRange = args.LSTRange
        self.LSTrRange = args.LSTrRange
        self.files_basename = list()
        self.files_filepath = list()
        self.files_extension = list()
        self.files = dict()
        self.files_full = dict()
        self.catalog = dict()
        self.catalog_full = dict()

        # Constants:
        self.BIN_WIDTH = 0.3 
        self.STOKES_POLS = ['I', 'Q', 'U', 'V']
        self.STANDARD_POLS = ['xx', 'xy', 'yx', 'yy']
        self.CLEAN = args.CLEAN

        # Internal method calls initiated upon instance creation
        self.logic()

    def logic(self):
        """Function to determine which process (catalog, wedge creation, plotting, etc.) to initiate."""
        if self.inputfiles:
            ares = mW.Ares(self)
            ares.makePlot()

        elif self.LSTRange or self.LSTrRange or self.JDRange or self.printc:
            self.format_pols()
            self.find_files()
            if self.step:
                for i in range(0, len(self.files_full[self.files_full.keys()[0]]), self.step):
                    self.catalog[self.ext] = dict()                    
                    for pol in self.files_full:
                        self.files[pol] = self.files_full[pol][i: i + self.step]
                        indices_files = np.array([], dtype=np.int64)
                        for j in range(len(self.files[pol])):
                            index = np.where(self.catalog_full[self.ext][pol]['file'] == self.files[pol][j])[0]
                            indices_files = np.append(indices_files, index)
                        self.catalog[self.ext][pol] = self.catalog_full[self.ext][pol][indices_files]

                    for pol in self.inputpols:
                        eris = mW.Eris(self, pol)
                        eris.name_npz()
                        eris.load_MIRIAD()
                        eris.pitchfork()
                        eris.save()

                    self.files = dict()
                    self.catalog = dict()
            else:
                self.files = self.files_full
                self.catalog = self.catalog_full
                for pol in self.inputpols:
                    eris = mW.Eris(self, pol)
                    eris.name_npz()
                    eris.load_MIRIAD()
                    eris.pitchfork()
                    eris.save()

        elif self.extension:
            self.catalog_directory()

    def catalog_directory(self):
        """Catalogs the MIRIAD files in a directory and looks for each file's polarization, JD, and LST,
        and saves that information in an npz file called 'catalog.npz'."""

        catalog = {ext: {} for ext in self.extension}

        for ext in self.extension:
            files_extension = glob.glob(self.filepath + '/zen.*' + ext)
            if len(files_extension) == 0:
                raise Exception("There are no files with extension '{}' in file path '{}'.".format(ext, self.filepath))

            # Cycle through the MIRIAD files in the directory and grab their info
            for file in files_extension:
                uv = UVData()
                uv.read_miriad(file)

                array_jd = np.array(np.unique(uv.time_array), dtype=np.float128)
                array_lstr = np.array(np.unique(uv.lst_array), dtype=np.float128)
                array_lst = aco.Angle(array_lstr*u.rad).hour
                pol_file = uv.get_pols()[0].lower()
                array_file = np.array([os.path.basename(file)] * array_jd.shape[0])

                array_cat = np.array([array_file, array_jd, array_lstr, array_lst])
                struct = np.dtype([('file', '|S32'),('JD', np.double),('LSTr', np.double),('LST', np.double)])
                array_cat = np.core.records.fromarrays(array_cat, dtype=struct)

                if pol_file in catalog[ext]:
                    catalog[ext][pol_file] = np.sort(np.append(catalog[ext][pol_file], array_cat), order='JD')
                else:
                    catalog[ext][pol_file] = np.array(array_cat)

        np.savez(os.path.join(self.filepath, 'catalog.npz'), cat=catalog)

    def find_files(self):
        """Properly format the given files specified directly or by JD, LST (hours (coming soon)), or LST (radiams)."""
        # Change directory to the directory with catalog.npz
        os.chdir(self.filepath)

        # Try to open catalog.npz
        try:
            self.catalog_full = np.load('catalog.npz')['cat'].tolist()
            if self.printc:
                print(self.catalog)
                quit()
        except IOError:
            raise Exception("There is no catalog (catalog.npz) in the specified path: %s" %self.filepath)

        self.ext = self.extension[0]

        if self.JDRange:
            JDRange_start, JDRange_stop = [float(x) for x in self.JDRange.split('_')]
        elif self.LSTrRange:
            LSTrRange_start, LSTrRange_stop = [float(x) for x in self.LSTrRange.split('_')]
        elif self.LSTRange:
            LSTRange_start, LSTRange_stop = [float(x) for x in self.LSTRange.split('_')]
        else:
            raise Exception("You must indicate which type of range you would like to use (JDRange, LSTrRange, LSTRange).")

        for ext in self.catalog_full.keys():
            if ext != self.ext:
                del self.catalog_full[ext]

        for pol in self.STANDARD_POLS:
            if (pol not in self.pol_dipole) and (pol in self.catalog_full[self.ext]):
                del self.catalog_full[self.ext][pol]
                continue

            if self.JDRange:
                indices = np.where(np.logical_and(self.catalog_full[self.ext][pol]['JD'] > JDRange_start,
                                                  self.catalog_full[self.ext][pol]['JD'] < JDRange_stop))[0]
            elif self.LSTrRange:
                indices = np.where(np.logical_and(self.catalog_full[self.ext][pol]['LSTr'] > LSTrRange_start,
                                                  self.catalog_full[self.ext][pol]['LSTr'] < LSTrRange_stop))[0]
            elif self.LSTRange:
                indices = np.where(np.logical_and(self.catalog_full[self.ext][pol]['LST'] > LSTRange_start,
                                                  self.catalog_full[self.ext][pol]['LST'] < LSTRange_stop))[0]
            else:
                raise Exception("You must indicate which type of range you would like to use (JDRange, LSTrRange, LSTRange).")

            files = np.unique(self.catalog_full[self.ext][pol][indices]['file']).tolist()
            indices_files = np.array([], dtype=np.int64)
            for i in range(len(files)):
                index = np.where(self.catalog_full[self.ext][pol]['file'] == files[i])[0]
                indices_files = np.append(indices_files, index)

            if len(indices):
                self.files_full[pol] = files
            else:
                raise Exception("There are no files in the specified time range. Try printing the catalog with --printc.")

            self.catalog_full[self.ext][pol] = self.catalog_full[self.ext][pol][indices_files]

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
