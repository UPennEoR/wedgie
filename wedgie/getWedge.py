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

def getWedge(filenames):
    
    import capo
    import matplotlib.pyplot as plt
    import numpy as np
    import baseline_lengths
    
    #get data from file
    t,d,f = capo.miriad.read_files(filenames,antstr='cross',polstr='xx')
    
    #create variable to store the wedge subplots in
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = baseline_lengths.get_baselines()
    
    #for each baselength in the dictionary
    for length in antdict.keys():
    
        totald = None #variable to store cumulative fft data
        
        #cycle through every ant pair of the given baselength
        for antpair in antdict[length]:
            
            #fft the data wrt freq and sum it with the previous iterations
            ftd_2D_data = np.fft.ifft(d[antpair]['xx'],axis=1)
            if totald == None:
                totald = ftd_2D_data
            else:
                totald += ftd_2D_data 
        
        #get average of all values for this baselength, store in wedgeslices
        totald /= len(antdict[length])
        wedgeslices.append(totald)
    
    #create subplot to plot data
    f,axarr = plt.subplots(len(wedgeslices),1,sharex=True,sharey=True)
    
    #add axes labels
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', 
                    right='off')
    plt.xlabel("Delay (bins)")
    plt.ylabel("Time (bins)")
    
    #plot individual wedge slices
    for i in range(len(wedgeslices)):
        axarr[i].imshow(np.log10(np.fft.fftshift(np.abs(wedgeslices[i]),
             axes=1)),aspect='auto',interpolation='nearest',vmin=-6,vmax=-1)
    
    #scale x axis to the significant information
    axarr[0].set_xlim(475,550)
    
    plt.show()

import argparse

#get filename from command line argument
parser = argparse.ArgumentParser()
parser.add_argument("filenames", help="your HERA data file(s)", nargs="*")
args=parser.parse_args()

#make wedge
getWedge(args.filenames)