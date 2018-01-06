"""
This script contains all of the logic necessary to call wedge_utils.py to plot or generate polarized pitchforks or wedges.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
Paul La Plante <plaplant_at_sas.upenn.edu>
"""

import argparse, os
import numpy as np
from pyuvdata import UVData

from IPython import embed


parser = argparse.ArgumentParser()

parser.add_argument("-T",
                    "--filetype",
                    help="Designate whether you will be inputing raw MIRIAD files or npz files, or if you want to use the catalog.",
                    choices=["MIRIAD", "npz", "catalog"],
                    default="MIRIAD")

# File Designating Args
parser.add_argument("-f",
                    "--filepath",
                    help="Designate the relative path to the folder that contains the files you want to be analyzed.",
                    required=True)

parser.add_argument("-P",
                    "--inputpols",
                    help="Input a comma-delimited list of polarizations to be used in analysis.",
                    default="I,Q,U,V")

parser.add_argument("-J",
                    "--JDRange",
                    help="Designate an inclusive range of JD to analyze. (e.g.: '2457548.3_2457549.5')")
parser.add_argument("-L",
                    "--LSTRange",
                    help="Designate an inclusive range of LST in hours to analyze. (e.g.: '13.0_15.5' will analyze from 1pm to 3:30pm)")
parser.add_argument("-r",
                    "--LSTrRange",
                    help="Designate an inclusive range of LST in hours to analyze. (e.g.: '2.7_3.7')")
parser.add_argument("-F",
                    "--inputfiles",
                    help="Designate specific files to be analyzed.",
                    nargs='*')


# Catalog Args
parser.add_argument("-k",
                    "--keyword",
                    help="Designate which file type (by any keyword in the files) you wish to catalog. (e.g.: '.uvcRK', '2457548')")


args = parser.parse_args()

class Zeus(object):
    def __init__(self, args):
        self.args = args
        self.filetype = args.filetype
        self.filepath = os.path.abspath(args.filepath)
        self.inputpols = args.inputpols
        self.JDRange = args.JDRange
        self.LSTrRange = args.LSTrRange
        self.LSTRange = args.LSTRange
        self.inputfiles = args.inputfiles
        self.keyword = args.keyword

        self.cwd = os.getcwd()

        self.pol_dipole = str()
        self.pol_type = str()

        self.files_basename = list()
        self.files_filepath = list()
        self.files_keyword = list()
        self.files = dict()

        self.STOKES_POLS = ['I', 'Q', 'U', 'V']
        self.STANDARD_POLS = ['xx', 'xy', 'yx', 'yy']

        self.logic()

    def logic(self):
        """Function to determine which process (catalog, wedge creation, plotting, etc.) to initiate."""
        if self.filetype == 'MIRIAD':
            self.format_pols()
            self.find_files()
        elif self.filetype == 'catalog':
            self.catalog()

    def catalog(self):
        """Catalogs the MIRIAD files in a directory and looks for each file's polarization, JD, and LST,
        and saves that information in an npz file called 'catalog.npz'."""
        catalog = np.array([[], [], [], []])

        files_filepath = os.listdir(self.filepath)
        files_keyword = sorted([file for file in files_filepath if self.keyword in file])

        """Maybe include a check here with os.path.splitext()"""

        uv = UVData()
        for file in files_keyword:
            uv.read_miriad(os.path.join(self.filepath, file))
            Ntimes = uv.Ntimes
            pol_file = uv.get_pols()[0].lower()
            array_lst = np.unique(uv.lst_array)
            array_jd = np.unique(uv.time_array)
            array_file = np.array([file] * Ntimes)
            array_pol = np.array([file_pol] * Ntimes)

            file_array = np.array((array_file, array_pol, array_jd.astype(float), array_lst.astype(float)))
            catalog = np.concatenate((catalog, file_array), axis=1)

        np.savez(os.path.join(self.filepath, 'catalog.npz'), cat=catalog.T)

    def find_files(self):
        """Properly format the given files specified directly or by JD, LST (hours (coming soon)), or LST (radiams)."""
        os.chdir(self.filepath)

        try:
            catalog = np.load('catalog.npz')['cat']
        except IOError:
            raise Exception("There is no catalog (catalog.npz) in the specified path: %s" %self.filepath)

        row_index = 0
        for pol in catalog[:, 1]:
            if not np.any(pol in self.pol_dipole):
                catalog = np.delete(catalog, row_index, axis=0)
                row_index -= 1
            row_index += 1

        if (self.JDRange is not None) or (self.LSTRange is not None) or (self.LSTrRange is not None):
            if self.JDRange is not None:
                JDRange_start, JDRange_stop = [float(x) for x in self.JDRange.split('_')]
                indices = np.where(np.logical_and(catalog[:, 2].astype(float) >= JDRange_start,
                                                  catalog[:, 2].astype(float) <= JDRange_stop))[0]
            elif self.LSTrRange is not None:
                LSTrRange_start, LSTrRange_stop = [float(x) for x in self.LSTrRange.split('_')]
                indices = np.where(np.logical_and(catalog[:, 3].astype(float) >= LSTrRange_start,
                                                  catalog[:, 3].astype(float) <= LSTrRange_stop))[0]
            catalog = np.take(catalog, indices, axis=0)

        elif self.inputfiles is not None:
            self.files_basename = [os.path.basename(file) for file in self.inputfiles]

            row_index = 0
            for file in catalog[:, 0]:
                if not np.any(file in self.files_basename):
                    catalog = np.delete(catalog, row_index, axis=0)
                    row_index -= 1
                row_index += 1

        files = catalog[:, :2]
        self.files = {pol: [] for pol in self.pol_dipole}
        for row in files:
            if row[0] not in self.files[row[1]]:
                self.files[row[1]].append(row[0])

    def format_pols(self):
        """Format the polarizations, e.g.: translating from IQUV to xx,xy,yx,yy"""
        self.inputpols = self.inputpols.split(',')
        self.pol_dipole = []

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
            raise Exception("You provided nonsensical polarization types: %s" %self.args.pol)

temp = Zeus(args)
embed()


# Wedge Creation Arguments
# parser.add_argument("-X",
#                     "--ex_ants",
#                     help="Input a comma-delimited list of antennae to exclude from analysis.")
# parser.add_argument("-C",
#                     "--calfile",
#                     help="Enter the calfile to be used for analysis. H0C data only.")


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

# parser.add_argument("-D",
#                     "--DataType",
#                     help="Designate abscal or sim data",
#                     choices=["regular", "abscal", "sim"],
#                     default="regular")

# parser.add_argument("-R",
#                     "--FreqRange",
#                     help="Designate with frequency band to analyze.")

# parser.add_argument("-p",
#                     "--path",
#                     help="Enter the path where you want save the files.")
