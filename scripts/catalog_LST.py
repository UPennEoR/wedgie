# -*- coding: utf-8 -*-
"""
Created on Sun Nov 05 16:54:58 2017

@author: Paul Chichura <pchich@sas.upenn.edu>

Contains two functions:
catalog_directory saves all LSTs in a given directory for later use
find_LST searches the saved LST catalog for files within a given LST range
"""

import argparse
import os
import numpy as np
from pyuvdata import UVData
uv = UVData()


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path',
                    help='path to the directory you want to catalogue by LST.',
                    default='.')
parser.add_argument('-e','--extension',
                    help='file extension you want listed in the catalogue.')
parser.add_argument('-c','--catalog',default=False,action='store_true',
                    help='toggle on catalog maker')
parser.add_argument('-r','--range',
                    help='LST1_LST2 you want to find appropriate files for')
parser.add_argument('-f','--find',default=False,action='store_true',
                    help='toggle on time finder')
args = parser.parse_args()

#function that catalogues a given directory
def catalog_directory(path, extension):
    #variable that stores data in the form:
    #   [['LST','filename','index_in_file','JD']]
    catalog = []
    
    #get the list of files to catalogue
    files = [] #list of all files ending with extension, only 1 file per time
    all_files = os.listdir(path)
    most_recent_file_time = ''
    for filename in all_files:
        if (filename[-1*len(extension):] == extension and 
                     filename.split('.')[2] != most_recent_file_time):
            files.append(filename)
            most_recent_file_time = filename.split('.')[2]
    
    #read through each file, get times
    for filename in files:
        uv.read_miriad(os.path.join(path, filename))
        LSTs = np.unique(uv.lst_array)
        indices = np.arange(0,uv.Ntimes)
        JDs = np.unique(uv.get_times(np.unique(uv.ant_1_array)[0],
                                     np.unique(uv.ant_2_array)[-1]))
        #save the times
        for index in indices:
            catalog.append([LSTs[index],
                            filename,
                            indices[index],
                            JDs[index]])
    
    #save the catalog
    np.savez(os.path.join(path, 'LST_catalog.npz'),cat=catalog)
if args.catalog: catalog_directory(args.path,args.extension)

#function that returns [['LST','filename','index_in_file','JD']] for files in a
#given path over a given LST range, in radians
def find_LST(LST_range, path='.'):
    LST_info = []
    LST_start = LST_range.split('_')[0]
    LST_stop = LST_range.split('_')[1]
    LST_catalog = np.load(os.path.join(path, 'LST_catalog.npz'))['cat']
    
    for i in range(len(LST_catalog)):
        if LST_catalog[i][0] >= LST_start and LST_catalog[i][0] <= LST_stop:
            LST_info.append(LST_catalog[i])

    return LST_info
if args.find: print(find_LST(args.range, args.path))
