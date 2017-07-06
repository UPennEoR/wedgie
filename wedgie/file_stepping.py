"""
This script analyzes HERA data files.  It analyzes the first file, then the 
first two files, then the first three files, and so on, until there are
as many plots as there are files.

There is a lot of hardcoded information in this file.  That might change.  The
directory of data, calfile, polarization, time averaged toggle and plot toggle are
all hardcoded.

Author: Austin Fox Fortino <fortino@sas.upenn.edu>
Created: July 5, 2017
Last Updated: July 5, 2017
"""

from os import system, listdir
from glob import glob
from pprint import pprint
from time import time

t0 = time()

files = []
for file in glob("/data4/paper/HERA2015/2457746/PennData/fourpol/*.OR"):
    files.append("/data4/paper/HERA2015/2457746/PennData/fourpol/" + file)
sorted(files)

calfile = "hsa7458_v001"
pol = "xx,xy,yx,yy"
time_avg = "--time_avg"
ex_ants = "--ex_ants=81"
args = [calfile, pol, time_avg, ex_ants]

troubleshoot_times = []

for file_index in range(len(files)):
    
    t1 = time()
    
    cmd = files[:file_index+1]
    cmd.extend(args)
    cmd = " ".join(cmd)
    system("python2.7 getWedge.py %s" %(cmd))
    
    t2 = time()
    t = t2-t1
    print t
    troubleshoot_times.append(t)

tf = time()

troubleshoot_times.append(sum(troubleshoot_times))
pprint("Troubleshot times in order:", troubleshoot_times)
print "Total time:", t0-tf